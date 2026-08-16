#' Compare state, clock, and dose predictors with held-out biological models
#'
#' Fits transparent linear models within leave-one-biological-model-out folds.
#' Input must already be at the biological-observation level; technical wells
#' are not permitted.
#'
#' @param data Data frame containing outcome and predictors.
#' @param outcome Numeric outcome column, for example history_effect.
#' @param model_col Biological model column defining held-out folds.
#' @param replicate_col Biological replicate identifier.
#' @param state_vars Measured biological-state predictors.
#' @param clock_var Elapsed-time predictor.
#' @param dose_vars Nominal-dose predictors.
#' @return A list containing held-out predictions and summary metrics.
#' @export
compare_state_clock_models <- function(
  data,
  outcome,
  model_col = "model",
  replicate_col = "replicate",
  state_vars,
  clock_var = "time",
  dose_vars = intersect(c("dose_A", "dose_B"), names(data))
) {
  .assert_columns(data, unique(c(outcome, model_col, replicate_col, state_vars, clock_var, dose_vars)))
  if ("technical" %in% names(data)) {
    stop("State-vs-clock modelling requires biological-observation-level data; collapse technical wells first.", call. = FALSE)
  }
  if (length(unique(as.character(data[[model_col]]))) < 3L) {
    warning("Fewer than three biological models: leave-one-model-out comparison is highly exploratory.", call. = FALSE)
  }

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
      complete_train <- stats::complete.cases(train[, unique(c(outcome, vars)), drop = FALSE])
      complete_test <- stats::complete.cases(test[, vars, drop = FALSE]) & is.finite(test[[outcome]])
      train2 <- train[complete_train, , drop = FALSE]
      test2 <- test[complete_test, , drop = FALSE]
      if (nrow(train2) <= length(vars) || !nrow(test2)) next
      f <- stats::reformulate(vars, response = outcome)
      fit <- try(stats::lm(f, data = train2), silent = TRUE)
      if (inherits(fit, "try-error")) next
      pr <- try(stats::predict(fit, newdata = test2), silent = TRUE)
      if (inherits(pr, "try-error")) next
      pred_rows[[length(pred_rows) + 1L]] <- data.frame(
        held_out_model = m,
        held_out_replicate = test2[[replicate_col]],
        predictor = nm,
        row_id = rownames(test2),
        observed = test2[[outcome]],
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
      n_observations = nrow(z),
      n_held_out_models = length(unique(z$held_out_model)),
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
#' Acute burden and delayed response are returned as separate quantities rather
#' than collapsed into an arbitrary weighted score. The per-hour quantity is
#' deliberately named a log response change; it is a biological growth rate
#' only when the assay is validated as proportional to cell number.
#'
#' @param data Long-form time-course response data.
#' @param acute_time Acute endpoint time.
#' @param late_time Delayed endpoint time after washout.
#' @param value Response column.
#' @param time_col Time column, interpreted in hours.
#' @return Matched acute and delayed responses with fold change and log response change per hour.
#' @export
evaluate_durable_fate <- function(data, acute_time, late_time, value = "response", time_col = "time") {
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", time_col, value))
  if (!is.numeric(data[[time_col]]) || any(!is.finite(data[[time_col]]))) stop("time must be finite numeric values.", call. = FALSE)
  if (!is.finite(acute_time) || !is.finite(late_time) || late_time <= acute_time) stop("late_time must be strictly greater than acute_time.", call. = FALSE)
  if ("technical" %in% names(data)) stop("Durable-fate analysis requires technical replicates to be collapsed first.", call. = FALSE)

  key <- c("model", "replicate", "sequence", "dose_A", "dose_B")
  target <- data[data[[time_col]] %in% c(acute_time, late_time), , drop = FALSE]
  dup_key <- c(key, time_col)
  if (any(duplicated(target[dup_key]))) {
    stop("Duplicate biological condition/time rows detected in durable-fate input.", call. = FALSE)
  }

  acute <- data[data[[time_col]] == acute_time, c(key, value), drop = FALSE]
  late <- data[data[[time_col]] == late_time, c(key, value), drop = FALSE]
  names(acute)[ncol(acute)] <- "acute_response"
  names(late)[ncol(late)] <- "late_response"
  out <- merge(acute, late, by = key, all = FALSE)
  exact <- is.finite(out$acute_response) & out$acute_response > 0 & is.finite(out$late_response) & out$late_response > 0
  out$fold_response_change <- NA_real_
  out$log_response_change_per_hour <- NA_real_
  out$fold_response_change[exact] <- out$late_response[exact] / out$acute_response[exact]
  out$log_response_change_per_hour[exact] <- log(out$fold_response_change[exact]) / (late_time - acute_time)
  out$durable_fate_status <- ifelse(exact, "exact", "censored_or_missing")
  out$interpretation_note <- "log_response_change_per_hour is a growth rate only if the assay is validated as proportional to cell number"
  out
}
