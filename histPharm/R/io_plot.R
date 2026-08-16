#' Export a fixed-first sequence slice in the public d-chain example schema
#'
#' This adapter creates AB, A0, and A rows compatible in column structure with
#' the public d-chain example dataset. It does not claim equivalence between
#' histPharm's reciprocal estimands and d-chain's Bayesian model assumptions.
#'
#' @param data histPharm response data.
#' @param sequence Sequence to export.
#' @param first_dose Fixed first-drug dose.
#' @param drug_A Drug A display name.
#' @param drug_B Drug B display name.
#' @param value Response column.
#' @param plate_col Optional plate identifier column.
#' @return A data frame with d-chain example column names.
#' @export
as_dchain_viability_data <- function(
  data,
  sequence = c("A_to_B", "B_to_A"),
  first_dose,
  drug_A = "DrugA",
  drug_B = "DrugB",
  value = "response",
  plate_col = NULL
) {
  sequence <- match.arg(sequence)
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", value))
  z <- data[data$sequence == sequence, , drop = FALSE]
  if (!nrow(z)) stop("Requested sequence is absent.", call. = FALSE)
  plate <- if (!is.null(plate_col) && plate_col %in% names(z)) z[[plate_col]] else z$replicate

  mk <- function(idx, experiment, pretreatment, compound, concentration) {
    if (!any(idx)) return(data.frame())
    data.frame(
      Experiment = experiment,
      CellLine = z$model[idx],
      Run = z$replicate[idx],
      Plate = plate[idx],
      Pretreatment = pretreatment,
      Compound = compound,
      Concentration = concentration(idx),
      RelCount = z[[value]][idx],
      stringsAsFactors = FALSE
    )
  }

  if (sequence == "A_to_B") {
    ab_idx <- z$dose_A == first_dose & z$dose_B > 0
    a0_idx <- z$dose_A == first_dose & z$dose_B == 0
    a_idx <- z$dose_A == 0 & z$dose_B > 0
    parts <- list(
      mk(ab_idx, "AB", drug_A, drug_B, function(i) z$dose_B[i]),
      mk(a0_idx, "A0", "None", drug_A, function(i) z$dose_A[i]),
      mk(a_idx, "A", "None", drug_B, function(i) z$dose_B[i])
    )
  } else {
    ab_idx <- z$dose_B == first_dose & z$dose_A > 0
    a0_idx <- z$dose_B == first_dose & z$dose_A == 0
    a_idx <- z$dose_B == 0 & z$dose_A > 0
    parts <- list(
      mk(ab_idx, "AB", drug_B, drug_A, function(i) z$dose_A[i]),
      mk(a0_idx, "A0", "None", drug_B, function(i) z$dose_B[i]),
      mk(a_idx, "A", "None", drug_A, function(i) z$dose_A[i])
    )
  }
  parts <- parts[vapply(parts, nrow, integer(1)) > 0]
  if (!length(parts)) return(data.frame())
  do.call(rbind, parts)
}

#' Plot a history-effect dose surface
#'
#' Uses measured grid values only. No smoothing or interpolation is performed.
#'
#' @param metrics Output of estimate_history_effect or a compatible table.
#' @param metric Numeric column to display.
#' @param main Optional title.
#' @return Aggregated plotting data invisibly.
#' @export
plot_history_surface <- function(metrics, metric = "history_effect", main = NULL) {
  .assert_columns(metrics, c("dose_A", "dose_B", metric))
  agg <- stats::aggregate(metrics[[metric]], metrics[c("dose_A", "dose_B")], mean, na.rm = TRUE)
  names(agg)[3] <- "value"
  mat <- stats::xtabs(value ~ dose_A + dose_B, data = agg)
  x <- as.numeric(rownames(mat)); y <- as.numeric(colnames(mat))
  graphics::image(x, y, unclass(mat), xlab = "Dose A", ylab = "Dose B", main = main)
  graphics::contour(x, y, unclass(mat), add = TRUE, drawlabels = FALSE)
  invisible(agg)
}

#' Plot interaction asymmetry against absolute efficacy
#'
#' The two axes deliberately remain separate. Quadrants reveal concordant and
#' discordant sequence regimes.
#'
#' @param reciprocal Output of estimate_reciprocal_asymmetry.
#' @param absolute Output of estimate_absolute_sequence_effect.
#' @param main Plot title.
#' @return Joined plotting data invisibly.
#' @export
plot_sequence_effect_map <- function(reciprocal, absolute, main = "Sequence effect map") {
  z <- classify_sequence_regime(reciprocal, absolute)
  graphics::plot(z$delta_history, z$absolute_log_ratio,
    xlab = "Reciprocal history asymmetry (DeltaH)",
    ylab = "Absolute sequence contrast",
    main = main, pch = 19)
  graphics::abline(h = 0, v = 0, lty = 2)
  invisible(z)
}
