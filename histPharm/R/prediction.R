#' Compare state, clock, and dose predictors with held-out biological models
#'
#' Fits simple linear models within leave-one-biological-model-out folds. This
#' deliberately favors transparent prediction and strict leakage control over
#' complex machine learning for small experimental datasets.
#'
#' @param data Data frame containing the outcome and predictors.
#' @param outcome Numeric outcome column, for example history_effect.
#' @param model_col Biological model column defining held-out folds.
#' @param state_vars Measured biological-state predictors.
#' @param clock_var Elapsed-time predictor.
#' @param dose_vars Nominal-dose predictors.
#' @return A list containing held-out predictions and summary metrics.
#' @export
compare_state_clock_models <- function(
  data,
  outcome,
  model_col = "model",
  state_vars,
  clock_var = "time",
  dose_vars = intersect(c("dose_A", "dose_B"), names(data))
) {
  .assert_columns(data, unique(c(outcome, model_col, state_vars, clock_var, dose_vars)))
  predictor_sets <- list(
    clock = clock_var,
    dose = dose_vars,
    state = state_vars,
    state_clock = unique(c(state_vars, clock_var))
  )
  predictor_sets <- predictor_sets[lengths(predictor_sets) > 0]
  held_out <- unique(as.character(data[[model_col]]))
  pred_rows <- list()

  for (m in held_out) {
    train <- data[as.character(data[[model_col]]) != m, , drop = FALSE]
    test <- data[as.character(data[[model_col]]) == m, , drop = FALSE]
    for (nm in names(predictor_sets)) {
      vars <- predictor_sets[[nm]]
      f <- stats::reformulate(vars, response = outcome)
      fit <- try(stats::lm(f, data = train), silent = TRUE)
      if (inherits(fit, "try-error")) next
      pr <- try(stats::predict(fit, newdata = test), silent = TRUE)
      if (inherits(pr, "try-error")) next
      pred_rows[[length(pred_rows) + 1L]] <- data.frame(
        held_out_model = m,
        predictor = nm,
        row_id = rownames(test),
        observed = test[[outcome]],
        predicted = as.numeric(pr),
        stringsAsFactors = FALSE
      )
    }
  }

  predictions <- if (length(pred_rows)) do.call(rbind, pred_rows) else data.frame()
  if (!nrow(predictions)) return(list(predictions = predictions, summary = data.frame()))
  summary <- .group_apply(predictions, "predictor", function(z) {
    ok <- is.finite(z$observed) & is.finite(z$predicted)
    z <- z[ok, , drop = FALSE]
    if (!nrow(z)) return(data.frame())
    residual <- z$observed - z$predicted
    sst <- sum((z$observed - mean(z$observed))^2)
    r2 <- if (sst > 0) 1 - sum(residual^2) / sst else NA_real_
    data.frame(
      predictor = z$predictor[1],
      n = nrow(z),
      rmse = sqrt(mean(residual^2)),
      mae = mean(abs(residual)),
      held_out_r2 = r2,
      stringsAsFactors = FALSE
    )
  })
  list(predictions = predictions, summary = summary)
}

#' Evaluate delayed post-washout regenerative fate
#'
#' Acute burden and delayed regrowth are returned as separate quantities rather
#' than collapsed into an arbitrary weighted score.
#'
#' @param data Long-form time-course response data.
#' @param acute_time Acute endpoint time.
#' @param late_time Delayed endpoint time after washout.
#' @param value Response column.
#' @param time_col Time column, interpreted in hours for the reported log growth rate.
#' @param floor Positive numerical floor.
#' @return Matched acute and delayed responses with fold regrowth and log growth rate.
#' @export
evaluate_durable_fate <- function(data, acute_time, late_time, value = "response", time_col = "time", floor = 1e-8) {
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", time_col, value))
  key <- c("model", "replicate", "sequence", "dose_A", "dose_B")
  acute <- data[data[[time_col]] == acute_time, c(key, value), drop = FALSE]
  late <- data[data[[time_col]] == late_time, c(key, value), drop = FALSE]
  names(acute)[ncol(acute)] <- "acute_response"
  names(late)[ncol(late)] <- "late_response"
  out <- merge(acute, late, by = key, all = FALSE)
  out$fold_regrowth <- pmax(out$late_response, floor) / pmax(out$acute_response, floor)
  out$log_growth_per_hour <- log(out$fold_regrowth) / (late_time - acute_time)
  out
}
