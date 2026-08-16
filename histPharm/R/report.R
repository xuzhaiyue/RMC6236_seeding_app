# Publication-oriented reporting helpers for reciprocal sequence experiments.
# All plots are raw-grid / raw-point first. No interpolation, smoothing, or
# imputation is performed by these functions.

.hp_finite_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

.hp_slug <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

.hp_surface_table <- function(data, value) {
  .assert_columns(data, c("dose_A", "dose_B", value))
  stats::aggregate(data[[value]], data[c("dose_A", "dose_B")], .hp_finite_mean)
}

.hp_draw_surface <- function(data, value, main, xlab = "Dose A", ylab = "Dose B") {
  z <- .hp_surface_table(data, value)
  names(z)[3] <- "value"
  xs <- sort(unique(z$dose_A))
  ys <- sort(unique(z$dose_B))
  mat <- matrix(NA_real_, nrow = length(xs), ncol = length(ys),
                dimnames = list(as.character(xs), as.character(ys)))
  for (i in seq_len(nrow(z))) {
    mat[match(z$dose_A[i], xs), match(z$dose_B[i], ys)] <- z$value[i]
  }
  graphics::image(xs, ys, mat, xlab = xlab, ylab = ylab, main = main)
  ok <- which(is.finite(z$value))
  if (length(ok)) graphics::points(z$dose_A[ok], z$dose_B[ok], pch = 21, bg = "white", cex = 0.7)
  attr(z, "surface_matrix") <- mat
  invisible(z)
}

#' Plot a measured response surface without interpolation
#'
#' @param data Normalized histPharm data.
#' @param value Response column.
#' @param model Optional model selector.
#' @param sequence Optional sequence selector.
#' @param main Optional title.
#' @return Aggregated measured-grid data invisibly.
#' @export
plot_response_surface <- function(data, value = "response", model = NULL,
                                  sequence = NULL, main = NULL) {
  .assert_columns(data, c("model", "sequence", "dose_A", "dose_B", value))
  z <- data
  if (!is.null(model)) z <- z[z$model == model, , drop = FALSE]
  if (!is.null(sequence)) z <- z[z$sequence == sequence, , drop = FALSE]
  if (!nrow(z)) stop("No rows remain after model/sequence selection.", call. = FALSE)
  if (is.null(main)) main <- paste(c(model, sequence, value), collapse = " | ")
  .hp_draw_surface(z, value, main)
}

#' Plot any replicate-level sequence metric on the measured dose grid
#'
#' Missing/censored coordinates remain missing; no smoothing is applied.
#'
#' @param metrics Metric table containing model, dose_A, dose_B.
#' @param metric Numeric metric column.
#' @param model Optional model selector.
#' @param sequence Optional sequence selector when present.
#' @param main Optional title.
#' @return Aggregated measured-grid data invisibly.
#' @export
plot_metric_surface <- function(metrics, metric, model = NULL,
                                sequence = NULL, main = NULL) {
  .assert_columns(metrics, c("model", "dose_A", "dose_B", metric))
  z <- metrics
  if (!is.null(model)) z <- z[z$model == model, , drop = FALSE]
  if (!is.null(sequence)) {
    .assert_columns(z, "sequence")
    z <- z[z$sequence == sequence, , drop = FALSE]
  }
  if (!nrow(z)) stop("No rows remain after model/sequence selection.", call. = FALSE)
  if (is.null(main)) main <- paste(c(model, sequence, metric), collapse = " | ")
  .hp_draw_surface(z, metric, main)
}

#' Build reciprocal mean-marginal effect-space coordinates
#'
#' Coordinates are visualization-only averages of the observed single-agent
#' residual fractions from the two reciprocal schedules. This does not alter
#' the estimand and must not be used as a substitute for replicate-level
#' inference.
#'
#' @param data Normalized reciprocal sequence data.
#' @param reciprocal Output of estimate_reciprocal_asymmetry().
#' @param value Response column.
#' @return Replicate-level effect-space table containing DeltaH.
#' @export
build_effect_space_projection <- function(data, reciprocal = NULL, value = "response") {
  if (is.null(reciprocal)) reciprocal <- estimate_reciprocal_asymmetry(data, value = value)
  es <- transform_effect_space(data, value = value)
  key <- c("model", "replicate", "dose_A", "dose_B")
  pieces <- split(seq_len(nrow(es)), .key_string(es, key))
  coords <- lapply(pieces, function(i) {
    z <- es[i, , drop = FALSE]
    data.frame(
      model = z$model[1], replicate = z$replicate[1],
      dose_A = z$dose_A[1], dose_B = z$dose_B[1],
      marginal_A_residual = .hp_finite_mean(z$marginal_A_residual),
      marginal_B_residual = .hp_finite_mean(z$marginal_B_residual),
      coordinate_method = "mean_reciprocal_marginal",
      stringsAsFactors = FALSE
    )
  })
  coords <- do.call(rbind, coords)
  rownames(coords) <- NULL
  merge(coords, reciprocal, by = key, all.x = TRUE)
}

#' Plot an effect-space projection from observed marginal responses
#'
#' Point location is the observed effect-space coordinate. Point size reflects
#' absolute metric magnitude and symbol direction reflects the sign. There is
#' no interpolation and no predefined "golden" region.
#'
#' @param projection Output of build_effect_space_projection().
#' @param metric Metric to encode, default DeltaH.
#' @param model Optional model selector.
#' @param main Optional title.
#' @return Plotted data invisibly.
#' @export
plot_effect_space_projection <- function(projection, metric = "delta_history",
                                         model = NULL, main = NULL) {
  .assert_columns(projection, c("model", "marginal_A_residual", "marginal_B_residual", metric))
  z <- projection
  if (!is.null(model)) z <- z[z$model == model, , drop = FALSE]
  ok <- is.finite(z$marginal_A_residual) & is.finite(z$marginal_B_residual) & is.finite(z[[metric]])
  z <- z[ok, , drop = FALSE]
  if (!nrow(z)) stop("No finite effect-space points are available.", call. = FALSE)
  mag <- abs(z[[metric]])
  denom <- max(mag, na.rm = TRUE)
  cex <- if (is.finite(denom) && denom > 0) 0.7 + 1.8 * mag / denom else rep(1, nrow(z))
  pch <- ifelse(z[[metric]] > 0, 24, ifelse(z[[metric]] < 0, 25, 21))
  if (is.null(main)) main <- paste(c(model, "Effect-space", metric), collapse = " | ")
  graphics::plot(z$marginal_A_residual, z$marginal_B_residual,
                 xlim = c(0, 1.05), ylim = c(0, 1.05),
                 xlab = "Drug A single-agent residual fraction",
                 ylab = "Drug B single-agent residual fraction",
                 main = main, pch = pch, cex = cex)
  graphics::abline(h = c(0.5, 0.7), v = c(0.5, 0.7), lty = 3)
  graphics::legend("topright", legend = c("metric > 0", "metric < 0", "metric = 0"),
                   pch = c(24, 25, 21), bty = "n", cex = 0.8)
  invisible(z)
}

#' Plot cross-model effect-space overlay
#'
#' Models are overlaid at their observed effect-space coordinates only. This
#' plot is descriptive and does not assert that high-effect regions overlap.
#'
#' @param projection Output of build_effect_space_projection().
#' @param metric Metric controlling point size.
#' @param models Optional models to include.
#' @param main Plot title.
#' @return Plotted data invisibly.
#' @export
plot_cross_model_effect_space <- function(projection, metric = "delta_history",
                                          models = NULL,
                                          main = "Cross-model effect-space overlay") {
  .assert_columns(projection, c("model", "marginal_A_residual", "marginal_B_residual", metric))
  z <- projection
  if (!is.null(models)) z <- z[z$model %in% models, , drop = FALSE]
  ok <- is.finite(z$marginal_A_residual) & is.finite(z$marginal_B_residual) & is.finite(z[[metric]])
  z <- z[ok, , drop = FALSE]
  if (!nrow(z)) stop("No finite cross-model effect-space points are available.", call. = FALSE)
  mods <- unique(as.character(z$model))
  symbols <- rep(c(21, 22, 23, 24, 25), length.out = length(mods))
  names(symbols) <- mods
  mag <- abs(z[[metric]])
  denom <- max(mag, na.rm = TRUE)
  cex <- if (is.finite(denom) && denom > 0) 0.65 + 1.5 * mag / denom else rep(1, nrow(z))
  graphics::plot(NA, xlim = c(0, 1.05), ylim = c(0, 1.05),
                 xlab = "Drug A single-agent residual fraction",
                 ylab = "Drug B single-agent residual fraction", main = main)
  for (m in mods) {
    i <- which(z$model == m)
    graphics::points(z$marginal_A_residual[i], z$marginal_B_residual[i],
                     pch = symbols[m], cex = cex[i])
  }
  graphics::abline(h = c(0.5, 0.7), v = c(0.5, 0.7), lty = 3)
  graphics::legend("topright", legend = mods, pch = symbols[mods], bty = "n", cex = 0.8)
  invisible(z)
}

.hp_plot_single_agent <- function(data, model, axis = c("A", "B"), value = "response") {
  axis <- match.arg(axis)
  z <- data[data$model == model, , drop = FALSE]
  if (axis == "A") {
    z <- z[z$dose_B == 0, , drop = FALSE]
    dose <- z$dose_A
    xlab <- "Drug A dose"
  } else {
    z <- z[z$dose_A == 0, , drop = FALSE]
    dose <- z$dose_B
    xlab <- "Drug B dose"
  }
  y <- z[[value]]
  ok <- is.finite(dose) & is.finite(y)
  dose <- dose[ok]; y <- y[ok]
  graphics::plot(dose, y, xlab = xlab, ylab = "Normalized residual response",
                 main = paste(model, axis, "single-agent response"), pch = 19)
  if (length(dose)) {
    agg <- stats::aggregate(y, list(dose = dose), .hp_finite_mean)
    graphics::lines(agg$dose, agg$x, lwd = 2)
    graphics::abline(h = 0.5, lty = 3)
  }
}

#' Generate the reciprocal discovery report package
#'
#' This is the one-command reporting entry point for a reciprocal dose matrix.
#' It validates raw data, collapses technical replicates, applies explicit LLOQ
#' masking, computes H_A->B, H_B->A, DeltaH and (when calibrated) the absolute
#' sequence contrast, writes machine-readable tables, and exports raw-first PDF
#' figures including fixed cross-model effect-space overlays.
#'
#' @param data Raw long-form histPharm data.
#' @param outdir Output directory.
#' @param value Raw readout column.
#' @param blank_col Optional row-wise blank column.
#' @param lloq Explicit LLOQ after blank correction.
#' @param lloq_col Optional row-wise LLOQ column.
#' @param baseline Normalization baseline.
#' @param require_absolute If TRUE, failure of calibrated absolute comparison is recorded, not silently ignored.
#' @param calibration_col Calibration identifier.
#' @return A list of analysis tables and manifest, invisibly.
#' @export
run_reciprocal_report <- function(data, outdir, value = "readout", blank_col = NULL,
                                  lloq = NULL, lloq_col = NULL,
                                  baseline = "within_sequence",
                                  require_absolute = TRUE,
                                  calibration_col = "calibration_group") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dirs <- file.path(outdir, c("01_QC", "02_Response", "03_History", "04_Sequence", "05_EffectSpace"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  validate_sequence_data(data)
  collapsed <- collapse_technical_replicates(data, value = value)
  norm <- normalize_sequence_data(collapsed, value = value, blank_col = blank_col,
                                  baseline = baseline, lloq = lloq, lloq_col = lloq_col)
  h <- estimate_history_effect(norm)
  reciprocal <- estimate_reciprocal_asymmetry(norm)
  absolute_error <- NULL
  absolute <- tryCatch(
    estimate_absolute_sequence_effect(norm, value = "response",
                                      calibration_col = calibration_col,
                                      require_calibration = require_absolute),
    error = function(e) { absolute_error <<- conditionMessage(e); data.frame() }
  )
  projection <- build_effect_space_projection(norm, reciprocal = reciprocal)

  utils::write.csv(collapsed, file.path(outdir, "01_QC", "collapsed_biological_conditions.csv"), row.names = FALSE)
  utils::write.csv(norm, file.path(outdir, "01_QC", "normalized_response_with_LLOQ_flags.csv"), row.names = FALSE)
  utils::write.csv(h, file.path(outdir, "03_History", "history_effects.csv"), row.names = FALSE)
  utils::write.csv(reciprocal, file.path(outdir, "04_Sequence", "reciprocal_deltaH.csv"), row.names = FALSE)
  if (nrow(absolute)) utils::write.csv(absolute, file.path(outdir, "04_Sequence", "absolute_sequence_contrast.csv"), row.names = FALSE)
  utils::write.csv(projection, file.path(outdir, "05_EffectSpace", "effect_space_projection.csv"), row.names = FALSE)

  manifest <- data.frame(file = character(), class = character(), status = character(), note = character(), stringsAsFactors = FALSE)
  add_manifest <- function(file, class, status = "generated", note = "") {
    manifest <<- rbind(manifest, data.frame(file = file, class = class, status = status, note = note, stringsAsFactors = FALSE))
  }

  mods <- unique(as.character(norm$model))
  seqs <- c("A_to_B", "B_to_A")
  for (m in mods) {
    f <- file.path(outdir, "02_Response", paste0(.hp_slug(m), "_single_agent_curves.pdf"))
    grDevices::pdf(f, width = 9, height = 4.5); old <- graphics::par(mfrow = c(1, 2));
    .hp_plot_single_agent(norm, m, "A"); .hp_plot_single_agent(norm, m, "B"); graphics::par(old); grDevices::dev.off()
    add_manifest(f, "single_agent")

    for (s in seqs) {
      f <- file.path(outdir, "02_Response", paste0(.hp_slug(m), "_", s, "_response_surface.pdf"))
      grDevices::pdf(f); plot_response_surface(norm, model = m, sequence = s, main = paste(m, s, "response")); grDevices::dev.off()
      add_manifest(f, "response_surface")

      f <- file.path(outdir, "03_History", paste0(.hp_slug(m), "_", s, "_H_surface.pdf"))
      grDevices::pdf(f); plot_metric_surface(h, "history_effect", model = m, sequence = s, main = paste(m, s, "history effect H")); grDevices::dev.off()
      add_manifest(f, "history_surface")
    }

    f <- file.path(outdir, "04_Sequence", paste0(.hp_slug(m), "_DeltaH_surface.pdf"))
    grDevices::pdf(f); plot_metric_surface(reciprocal, "delta_history", model = m, main = paste(m, "reciprocal DeltaH")); grDevices::dev.off()
    add_manifest(f, "deltaH_surface")

    if (nrow(absolute)) {
      f <- file.path(outdir, "04_Sequence", paste0(.hp_slug(m), "_absolute_sequence_surface.pdf"))
      grDevices::pdf(f); plot_metric_surface(absolute, "absolute_log_ratio", model = m, main = paste(m, "absolute sequence contrast")); grDevices::dev.off()
      add_manifest(f, "absolute_surface")

      f <- file.path(outdir, "04_Sequence", paste0(.hp_slug(m), "_interaction_vs_absolute.pdf"))
      grDevices::pdf(f); plot_sequence_effect_map(reciprocal[reciprocal$model == m, , drop = FALSE],
                                                  absolute[absolute$model == m, , drop = FALSE],
                                                  main = paste(m, "interaction vs absolute")); grDevices::dev.off()
      add_manifest(f, "interaction_vs_absolute")
    }

    f <- file.path(outdir, "05_EffectSpace", paste0(.hp_slug(m), "_effect_space_DeltaH.pdf"))
    grDevices::pdf(f); plot_effect_space_projection(projection, model = m, main = paste(m, "effect-space DeltaH")); grDevices::dev.off()
    add_manifest(f, "effect_space")
  }

  f <- file.path(outdir, "05_EffectSpace", "CROSS_MODEL_effect_space_overlay_DeltaH.pdf")
  grDevices::pdf(f, width = 7, height = 7); plot_cross_model_effect_space(projection); grDevices::dev.off()
  add_manifest(f, "cross_model_effect_space", note = "Raw observed effect-space coordinates; point size encodes |DeltaH|; no interpolation.")

  if (!nrow(absolute) && !is.null(absolute_error)) {
    add_manifest("absolute_sequence_contrast", "absolute_surface", "not_generated", absolute_error)
  }
  utils::write.csv(manifest, file.path(outdir, "report_manifest.csv"), row.names = FALSE)
  result <- list(collapsed = collapsed, normalized = norm, history = h,
                 reciprocal = reciprocal, absolute = absolute,
                 effect_space = projection, manifest = manifest,
                 absolute_error = absolute_error)
  invisible(result)
}
