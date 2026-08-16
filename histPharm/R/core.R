.histpharm_required <- c("model", "replicate", "sequence", "dose_A", "dose_B", "readout")
.histpharm_timing <- c("phase1_start", "phase1_end", "phase2_start", "phase2_end", "endpoint_time")

.assert_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

.safe_log <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- log(x[ok])
  out
}

.key_string <- function(data, cols) {
  do.call(paste, c(data[cols], sep = "|||"))
}

.group_apply <- function(data, group_cols, FUN) {
  if (!length(group_cols)) return(FUN(data))
  keys <- .key_string(data, group_cols)
  out <- lapply(split(seq_len(nrow(data)), keys), function(i) FUN(data[i, , drop = FALSE]))
  out <- out[lengths(out) > 0]
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

.lookup_value <- function(data, dose_A, dose_B, value) {
  hit <- data$dose_A == dose_A & data$dose_B == dose_B
  if (sum(hit) != 1L) return(NA_real_)
  as.numeric(data[[value]][hit])
}

.lookup_flag <- function(data, dose_A, dose_B, value) {
  if (!value %in% names(data)) return(FALSE)
  hit <- data$dose_A == dose_A & data$dose_B == dose_B
  if (sum(hit) != 1L) return(TRUE)
  isTRUE(as.logical(data[[value]][hit]))
}

.assert_collapsed_unique <- function(data, extra_cols = character()) {
  key <- unique(c("model", "replicate", "sequence", "dose_A", "dose_B", extra_cols))
  .assert_columns(data, key)
  dup <- duplicated(data[key]) | duplicated(data[key], fromLast = TRUE)
  if (any(dup)) {
    stop(
      "Duplicate biological conditions detected after technical-replicate collapse. ",
      "Each model/replicate/sequence/dose pair/endpoint must be unique before estimation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate reciprocal sequence-response data
#'
#' Publication-grade validation checks the complete reciprocal grid, raw
#' technical-replicate identifiers, timing metadata, plate identifiers, and
#' matched schedule timing. Technical replicates are allowed only when an
#' explicit technical identifier is present.
#'
#' @param data A data frame.
#' @param strict If TRUE, validation problems are errors; otherwise warnings.
#' @param require_complete_grid Require every dose_A x dose_B coordinate.
#' @param require_timing Require phase timing metadata and matched reciprocal timing.
#' @param require_plate Require a non-missing plate identifier.
#' @param technical_col Technical-replicate identifier column.
#' @return The input data invisibly.
#' @export
validate_sequence_data <- function(
  data,
  strict = TRUE,
  require_complete_grid = TRUE,
  require_timing = TRUE,
  require_plate = TRUE,
  technical_col = "technical"
) {
  .assert_columns(data, .histpharm_required)
  if (!is.numeric(data$dose_A) || !is.numeric(data$dose_B)) stop("dose_A and dose_B must be numeric.", call. = FALSE)
  if (!is.numeric(data$readout)) stop("readout must be numeric.", call. = FALSE)
  if (any(!is.finite(data$dose_A)) || any(!is.finite(data$dose_B))) stop("Doses must be finite.", call. = FALSE)
  if (any(data$dose_A < 0 | data$dose_B < 0)) stop("Doses cannot be negative.", call. = FALSE)
  if (any(!is.finite(data$readout)) || any(data$readout < 0)) stop("readout must be finite and non-negative before blank correction.", call. = FALSE)
  bad_seq <- setdiff(unique(as.character(data$sequence)), c("A_to_B", "B_to_A"))
  if (length(bad_seq)) stop("sequence must contain only A_to_B or B_to_A.", call. = FALSE)

  problems <- character()

  if (require_plate) {
    if (!"plate" %in% names(data)) {
      problems <- c(problems, "missing required plate metadata")
    } else if (any(is.na(data$plate) | !nzchar(as.character(data$plate)))) {
      problems <- c(problems, "plate identifiers must be non-missing")
    }
  }

  if (require_timing) {
    missing_timing <- setdiff(.histpharm_timing, names(data))
    if (length(missing_timing)) {
      problems <- c(problems, paste0("missing timing metadata: ", paste(missing_timing, collapse = ", ")))
    } else {
      for (nm in .histpharm_timing) {
        if (!is.numeric(data[[nm]]) || any(!is.finite(data[[nm]]))) {
          problems <- c(problems, paste0(nm, " must be finite numeric time metadata"))
        }
      }
      if (!length(problems)) {
        bad_order <- !(data$phase1_start < data$phase1_end &
          data$phase1_end <= data$phase2_start &
          data$phase2_start < data$phase2_end &
          data$phase2_end <= data$endpoint_time)
        if (any(bad_order)) problems <- c(problems, "invalid phase timing order")
      }
    }
  }

  condition_cols <- c("model", "replicate", "sequence", "dose_A", "dose_B")
  if ("endpoint" %in% names(data)) condition_cols <- c(condition_cols, "endpoint")
  keys <- .key_string(data, condition_cols)
  for (idx in split(seq_len(nrow(data)), keys)) {
    if (length(idx) > 1L) {
      if (!technical_col %in% names(data)) {
        problems <- c(problems, "duplicate biological condition without explicit technical-replicate identifiers")
      } else {
        tech <- data[[technical_col]][idx]
        if (any(is.na(tech)) || anyDuplicated(tech)) {
          problems <- c(problems, "technical replicate identifiers must be present and unique within each biological condition")
        }
      }
    }
  }

  blocks <- interaction(data$model, data$replicate, data$sequence, drop = TRUE)
  for (idx in split(seq_len(nrow(data)), blocks)) {
    z <- data[idx, , drop = FALSE]
    tag <- paste(z$model[1], z$replicate[1], z$sequence[1], sep = "/")
    if (!any(z$dose_A == 0 & z$dose_B == 0)) {
      problems <- c(problems, paste0(tag, ": missing (0,0) vehicle"))
      next
    }
    if (!0 %in% z$dose_A || !0 %in% z$dose_B) {
      problems <- c(problems, paste0(tag, ": both dose grids must include zero"))
    }
    for (a in unique(z$dose_A[z$dose_A > 0])) {
      if (!any(z$dose_A == a & z$dose_B == 0)) problems <- c(problems, paste0(tag, ": missing A-only margin for dose_A=", a))
    }
    for (b in unique(z$dose_B[z$dose_B > 0])) {
      if (!any(z$dose_A == 0 & z$dose_B == b)) problems <- c(problems, paste0(tag, ": missing B-only margin for dose_B=", b))
    }
    if (require_complete_grid) {
      combos <- unique(z[c("dose_A", "dose_B")])
      expected <- length(unique(z$dose_A)) * length(unique(z$dose_B))
      if (nrow(combos) != expected) {
        problems <- c(problems, paste0(tag, ": incomplete dose grid; expected ", expected, " unique coordinates but found ", nrow(combos)))
      }
    }
    if (require_timing && all(.histpharm_timing %in% names(z))) {
      for (nm in .histpharm_timing) {
        if (length(unique(z[[nm]])) != 1L) problems <- c(problems, paste0(tag, ": ", nm, " varies within a sequence block"))
      }
    }
  }

  reciprocal_blocks <- interaction(data$model, data$replicate, drop = TRUE)
  for (idx in split(seq_len(nrow(data)), reciprocal_blocks)) {
    z <- data[idx, , drop = FALSE]
    seqs <- unique(as.character(z$sequence))
    tag <- paste(z$model[1], z$replicate[1], sep = "/")
    if (!all(c("A_to_B", "B_to_A") %in% seqs)) {
      problems <- c(problems, paste0(tag, ": reciprocal sequence pair is incomplete"))
      next
    }
    za <- z[z$sequence == "A_to_B", , drop = FALSE]
    zb <- z[z$sequence == "B_to_A", , drop = FALSE]
    if (!setequal(unique(za$dose_A), unique(zb$dose_A)) || !setequal(unique(za$dose_B), unique(zb$dose_B))) {
      problems <- c(problems, paste0(tag, ": reciprocal sequences use different dose grids"))
    }
    if (require_timing && all(.histpharm_timing %in% names(z))) {
      ta <- as.numeric(za[1, .histpharm_timing, drop = TRUE])
      tb <- as.numeric(zb[1, .histpharm_timing, drop = TRUE])
      if (!isTRUE(all.equal(ta, tb, tolerance = 0))) {
        problems <- c(problems, paste0(tag, ": reciprocal schedules are not time matched"))
      }
    }
  }

  if (length(problems)) {
    msg <- paste(unique(problems), collapse = "\n")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(data)
}

#' Collapse technical replicates explicitly
#'
#' @param data A histPharm long-form data frame.
#' @param value Name of the response column.
#' @param id_cols Columns defining a biological condition.
#' @param fun Aggregation function.
#' @return One row per biological condition.
#' @export
collapse_technical_replicates <- function(
  data,
  value = "readout",
  id_cols = intersect(
    c("model", "replicate", "sequence", "dose_A", "dose_B", "endpoint", "time",
      "plate", "calibration_group", .histpharm_timing),
    names(data)
  ),
  fun = mean
) {
  .assert_columns(data, c(id_cols, value))
  key <- .key_string(data, id_cols)
  pieces <- lapply(split(seq_len(nrow(data)), key), function(i) {
    z <- data[i, , drop = FALSE]
    out <- z[1, id_cols, drop = FALSE]
    out[[value]] <- fun(z[[value]], na.rm = TRUE)
    out$technical_n <- sum(is.finite(z[[value]]))
    out$technical_sd <- if (out$technical_n > 1) stats::sd(z[[value]], na.rm = TRUE) else NA_real_
    out
  })
  ans <- do.call(rbind, pieces)
  rownames(ans) <- NULL
  .assert_collapsed_unique(ans, extra_cols = intersect(c("endpoint", "time"), names(ans)))
  ans
}

#' Normalize sequence-response data with explicit LLOQ handling
#'
#' Values at or below the LLOQ are censored by default and are never replaced
#' by an arbitrary tiny positive constant. Exact history effects are therefore
#' unavailable when a required component is below quantification.
#'
#' @param data A histPharm data frame.
#' @param value Response column.
#' @param blank_col Optional row-wise blank column.
#' @param baseline within_sequence, shared_replicate, or none.
#' @param output Output column name.
#' @param lloq Optional scalar lower limit of quantification after blank correction.
#' @param lloq_col Optional row-wise LLOQ column.
#' @param lloq_action mask or error.
#' @return Data with normalized response and censoring metadata.
#' @export
normalize_sequence_data <- function(
  data,
  value = "readout",
  blank_col = NULL,
  baseline = c("within_sequence", "shared_replicate", "none"),
  output = "response",
  lloq = NULL,
  lloq_col = NULL,
  lloq_action = c("mask", "error")
) {
  baseline <- match.arg(baseline)
  lloq_action <- match.arg(lloq_action)
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", value))
  .assert_collapsed_unique(data, extra_cols = intersect(c("endpoint", "time"), names(data)))
  if (!is.null(lloq) && !is.null(lloq_col)) stop("Specify either lloq or lloq_col, not both.", call. = FALSE)

  x <- as.numeric(data[[value]])
  if (!is.null(blank_col)) {
    .assert_columns(data, blank_col)
    x <- x - as.numeric(data[[blank_col]])
  }

  threshold <- rep(NA_real_, length(x))
  if (!is.null(lloq_col)) {
    .assert_columns(data, lloq_col)
    threshold <- as.numeric(data[[lloq_col]])
    if (any(!is.finite(threshold) | threshold < 0)) stop("LLOQ values must be finite and non-negative.", call. = FALSE)
  } else if (!is.null(lloq)) {
    if (length(lloq) != 1L || !is.finite(lloq) || lloq < 0) stop("lloq must be one finite non-negative number.", call. = FALSE)
    threshold[] <- as.numeric(lloq)
  }

  if (any(x <= 0, na.rm = TRUE) && all(is.na(threshold))) {
    stop(
      "Non-positive values remain after blank correction. Provide an explicit LLOQ; histPharm no longer floors such values to a tiny positive constant.",
      call. = FALSE
    )
  }

  below <- ifelse(is.na(threshold), FALSE, x <= threshold)
  below <- below | !is.finite(x) | x <= 0
  if (lloq_action == "error" && any(below)) stop("One or more observations are at/below LLOQ.", call. = FALSE)

  used <- x
  used[below] <- NA_real_
  data$blank_corrected_readout <- x
  data$lloq_value <- threshold
  data$below_lloq <- below

  if (baseline == "none") {
    data[[output]] <- used
  } else if (baseline == "within_sequence") {
    grp <- interaction(data$model, data$replicate, data$sequence, drop = TRUE)
    base <- ave(seq_len(nrow(data)), grp, FUN = function(i) {
      z <- data[i, , drop = FALSE]
      v <- used[i][z$dose_A == 0 & z$dose_B == 0]
      v <- v[is.finite(v)]
      if (!length(v)) return(rep(NA_real_, length(i)))
      rep(mean(v), length(i))
    })
    data[[output]] <- used / base
  } else {
    grp <- interaction(data$model, data$replicate, drop = TRUE)
    base <- ave(seq_len(nrow(data)), grp, FUN = function(i) {
      z <- data[i, , drop = FALSE]
      v <- used[i][z$dose_A == 0 & z$dose_B == 0]
      v <- v[is.finite(v)]
      if (!length(v)) return(rep(NA_real_, length(i)))
      rep(mean(v), length(i))
    })
    data[[output]] <- used / base
  }
  data$baseline_method <- baseline
  data
}

#' Estimate history-conditioned incremental drug effects
#'
#' @param data Matched, collapsed sequence data.
#' @param value Response column.
#' @return Replicate-level exact history effects; censored/missing contrasts are NA.
#' @export
estimate_history_effect <- function(data, value = "response") {
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", value))
  .assert_collapsed_unique(data, extra_cols = intersect(c("endpoint", "time"), names(data)))
  .group_apply(data, c("model", "replicate", "sequence"), function(z) {
    pairs <- unique(z[z$dose_A > 0 & z$dose_B > 0, c("dose_A", "dose_B"), drop = FALSE])
    if (!nrow(pairs)) return(data.frame())
    out <- lapply(seq_len(nrow(pairs)), function(i) {
      a <- pairs$dose_A[i]
      b <- pairs$dose_B[i]
      vals <- c(
        .lookup_value(z, 0, 0, value),
        .lookup_value(z, a, 0, value),
        .lookup_value(z, 0, b, value),
        .lookup_value(z, a, b, value)
      )
      censored <- c(
        .lookup_flag(z, 0, 0, "below_lloq"),
        .lookup_flag(z, a, 0, "below_lloq"),
        .lookup_flag(z, 0, b, "below_lloq"),
        .lookup_flag(z, a, b, "below_lloq")
      )
      exact <- all(is.finite(vals) & vals > 0) && !any(censored)
      if (exact && z$sequence[1] == "A_to_B") {
        with_history <- log(vals[2] / vals[4])
        no_history <- log(vals[1] / vals[3])
        first_drug <- "A"
        second_drug <- "B"
      } else if (exact) {
        with_history <- log(vals[3] / vals[4])
        no_history <- log(vals[1] / vals[2])
        first_drug <- "B"
        second_drug <- "A"
      } else {
        with_history <- NA_real_
        no_history <- NA_real_
        first_drug <- if (z$sequence[1] == "A_to_B") "A" else "B"
        second_drug <- if (z$sequence[1] == "A_to_B") "B" else "A"
      }
      data.frame(
        model = z$model[1], replicate = z$replicate[1], sequence = z$sequence[1],
        dose_A = a, dose_B = b, first_drug = first_drug, second_drug = second_drug,
        second_effect_with_history = with_history,
        second_effect_without_history = no_history,
        history_effect = if (exact) with_history - no_history else NA_real_,
        history_effect_status = if (exact) "exact" else "censored_or_missing",
        n_censored_components = sum(censored | !is.finite(vals) | vals <= 0),
        V00 = vals[1], Va0 = vals[2], V0b = vals[3], Vab = vals[4],
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, out)
  })
}

#' Estimate reciprocal history asymmetry
#' @param data Reciprocal sequence data.
#' @param value Response column.
#' @return Replicate-level reciprocal contrasts.
#' @export
estimate_reciprocal_asymmetry <- function(data, value = "response") {
  h <- estimate_history_effect(data, value = value)
  a <- h[h$sequence == "A_to_B", c("model", "replicate", "dose_A", "dose_B", "history_effect", "history_effect_status")]
  b <- h[h$sequence == "B_to_A", c("model", "replicate", "dose_A", "dose_B", "history_effect", "history_effect_status")]
  names(a)[5:6] <- c("H_A_to_B", "status_A_to_B")
  names(b)[5:6] <- c("H_B_to_A", "status_B_to_A")
  out <- merge(a, b, by = c("model", "replicate", "dose_A", "dose_B"), all = FALSE)
  out$delta_history <- out$H_A_to_B - out$H_B_to_A
  out$delta_status <- ifelse(out$status_A_to_B == "exact" & out$status_B_to_A == "exact", "exact", "censored_or_missing")
  out
}

#' Estimate absolute sequence effect
#'
#' Absolute cross-sequence comparisons require an explicit shared calibration
#' group linking reciprocal plates/schedules.
#'
#' @param data Reciprocal sequence data.
#' @param value Cross-sequence-comparable response column.
#' @param calibration_col Shared calibration identifier.
#' @param require_calibration Require explicit cross-plate calibration metadata.
#' @return Replicate-level absolute sequence contrasts.
#' @export
estimate_absolute_sequence_effect <- function(
  data,
  value = "response",
  calibration_col = "calibration_group",
  require_calibration = TRUE
) {
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", value))
  .assert_collapsed_unique(data, extra_cols = intersect(c("endpoint", "time"), names(data)))
  if (require_calibration && !calibration_col %in% names(data)) {
    stop("Absolute cross-sequence comparison requires explicit calibration_group metadata.", call. = FALSE)
  }
  z <- data[data$dose_A > 0 & data$dose_B > 0, , drop = FALSE]
  keep <- c("model", "replicate", "dose_A", "dose_B", value)
  if (calibration_col %in% names(z)) keep <- c(keep, calibration_col)
  a <- z[z$sequence == "A_to_B", keep, drop = FALSE]
  b <- z[z$sequence == "B_to_A", keep, drop = FALSE]
  names(a)[names(a) == value] <- "V_A_to_B"
  names(b)[names(b) == value] <- "V_B_to_A"
  if (calibration_col %in% names(z)) {
    names(a)[names(a) == calibration_col] <- "calibration_A_to_B"
    names(b)[names(b) == calibration_col] <- "calibration_B_to_A"
  }
  out <- merge(a, b, by = c("model", "replicate", "dose_A", "dose_B"), all = FALSE)
  if (require_calibration) {
    bad <- is.na(out$calibration_A_to_B) | is.na(out$calibration_B_to_A) |
      as.character(out$calibration_A_to_B) != as.character(out$calibration_B_to_A)
    if (any(bad)) stop("Reciprocal observations are not linked by the same calibration group.", call. = FALSE)
  }
  exact <- is.finite(out$V_A_to_B) & out$V_A_to_B > 0 & is.finite(out$V_B_to_A) & out$V_B_to_A > 0
  out$absolute_log_ratio <- NA_real_
  out$absolute_log_ratio[exact] <- log(out$V_B_to_A[exact]) - log(out$V_A_to_B[exact])
  out$absolute_status <- ifelse(exact, "exact", "censored_or_missing")
  out$preferred_absolute_sequence <- ifelse(
    !exact, NA_character_,
    ifelse(out$absolute_log_ratio > 0, "A_to_B", ifelse(out$absolute_log_ratio < 0, "B_to_A", "tie"))
  )
  out
}

#' Classify joint interaction-efficacy sequence regimes
#' @param reciprocal Reciprocal history table.
#' @param absolute Absolute sequence table.
#' @param tolerance Neutrality threshold.
#' @return Joint sequence-regime table.
#' @export
classify_sequence_regime <- function(reciprocal, absolute, tolerance = 0) {
  key <- c("model", "replicate", "dose_A", "dose_B")
  .assert_columns(reciprocal, c(key, "delta_history"))
  .assert_columns(absolute, c(key, "absolute_log_ratio"))
  z <- merge(reciprocal, absolute, by = key, all = FALSE)
  valid <- is.finite(z$delta_history) & is.finite(z$absolute_log_ratio)
  sh <- ifelse(z$delta_history > tolerance, 1L, ifelse(z$delta_history < -tolerance, -1L, 0L))
  sa <- ifelse(z$absolute_log_ratio > tolerance, 1L, ifelse(z$absolute_log_ratio < -tolerance, -1L, 0L))
  z$sequence_regime <- NA_character_
  z$sequence_regime[valid] <- ifelse(
    sh[valid] == 1L & sa[valid] == 1L, "A_to_B_concordant",
    ifelse(sh[valid] == -1L & sa[valid] == -1L, "B_to_A_concordant",
      ifelse(sh[valid] == 1L & sa[valid] == -1L, "A_to_B_interaction_B_to_A_absolute",
        ifelse(sh[valid] == -1L & sa[valid] == 1L, "B_to_A_interaction_A_to_B_absolute",
          ifelse(sh[valid] == 0L & sa[valid] == 0L, "neutral", "partially_neutral")
        )
      )
    )
  )
  z
}

#' Transform nominal doses to observed single-agent effect space
#' @param data Sequence-response data.
#' @param value Response column.
#' @return Effect-space coordinates.
#' @export
transform_effect_space <- function(data, value = "response") {
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", value))
  .assert_collapsed_unique(data, extra_cols = intersect(c("endpoint", "time"), names(data)))
  .group_apply(data, c("model", "replicate", "sequence"), function(z) {
    comb <- z[z$dose_A > 0 & z$dose_B > 0, , drop = FALSE]
    if (!nrow(comb)) return(data.frame())
    v00 <- .lookup_value(z, 0, 0, value)
    ans <- lapply(seq_len(nrow(comb)), function(i) {
      va0 <- .lookup_value(z, comb$dose_A[i], 0, value)
      v0b <- .lookup_value(z, 0, comb$dose_B[i], value)
      exact <- all(is.finite(c(v00, va0, v0b))) && v00 > 0
      ar <- if (exact) va0 / v00 else NA_real_
      br <- if (exact) v0b / v00 else NA_real_
      data.frame(
        model = comb$model[i], replicate = comb$replicate[i], sequence = comb$sequence[i],
        dose_A = comb$dose_A[i], dose_B = comb$dose_B[i],
        marginal_A_residual = ar, marginal_B_residual = br,
        marginal_A_suppression = if (is.finite(ar)) 1 - ar else NA_real_,
        marginal_B_suppression = if (is.finite(br)) 1 - br else NA_real_,
        combination_response = comb[[value]][i],
        effect_space_status = if (exact) "exact" else "censored_or_missing",
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, ans)
  })
}

#' Bootstrap uncertainty for sequence metrics
#'
#' Resampling is performed only over biological replicate-level metrics. For
#' fewer than min_replicates, output is explicitly labelled exploratory.
#'
#' @param metrics Replicate-level metric table.
#' @param metric Numeric metric column.
#' @param strata Surface strata.
#' @param B Bootstrap draws.
#' @param conf Confidence level.
#' @param seed Optional random seed.
#' @param min_replicates Minimum n for non-exploratory labelling.
#' @return Bootstrap estimates, intervals, and stability labels.
#' @export
bootstrap_sequence_metrics <- function(
  metrics,
  metric,
  strata = intersect(c("model", "dose_A", "dose_B"), names(metrics)),
  B = 2000,
  conf = 0.95,
  seed = NULL,
  min_replicates = 4
) {
  .assert_columns(metrics, c(strata, "replicate", metric))
  if ("technical" %in% names(metrics)) stop("Bootstrap input must be biological-replicate-level metrics, not technical wells.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  alpha <- (1 - conf) / 2
  .group_apply(metrics, strata, function(z) {
    ok <- is.finite(z[[metric]])
    z <- z[ok, , drop = FALSE]
    if (!nrow(z)) return(data.frame())
    reps <- unique(z$replicate)
    if (nrow(z) != length(reps)) stop("Each bootstrap stratum must contain one metric per biological replicate.", call. = FALSE)
    x <- z[[metric]]
    boots <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
    out <- z[1, strata, drop = FALSE]
    out$metric <- metric
    out$estimate <- mean(x)
    out$ci_low <- as.numeric(stats::quantile(boots, alpha, names = FALSE))
    out$ci_high <- as.numeric(stats::quantile(boots, 1 - alpha, names = FALSE))
    out$prob_gt_zero <- mean(boots > 0)
    out$n_biological_replicates <- length(x)
    out$inference_grade <- if (length(x) >= min_replicates) "bootstrap_descriptive" else "exploratory_small_n"
    out
  })
}
