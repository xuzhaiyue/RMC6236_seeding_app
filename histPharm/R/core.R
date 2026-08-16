.histpharm_required <- c("model", "replicate", "sequence", "dose_A", "dose_B", "readout")

.assert_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

.safe_log <- function(x, floor = 1e-8) {
  log(pmax(as.numeric(x), floor))
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

#' Validate reciprocal sequence-response data
#'
#' @param data A data frame.
#' @param strict If TRUE, missing matched margins are errors; otherwise warnings.
#' @return The input data invisibly.
#' @export
validate_sequence_data <- function(data, strict = TRUE) {
  .assert_columns(data, .histpharm_required)
  if (!is.numeric(data$dose_A) || !is.numeric(data$dose_B)) stop("dose_A and dose_B must be numeric.", call. = FALSE)
  if (!is.numeric(data$readout)) stop("readout must be numeric.", call. = FALSE)
  if (any(!is.finite(data$dose_A)) || any(!is.finite(data$dose_B))) stop("Doses must be finite.", call. = FALSE)
  if (any(data$dose_A < 0 | data$dose_B < 0)) stop("Doses cannot be negative.", call. = FALSE)
  if (any(!is.finite(data$readout)) || any(data$readout < 0)) stop("readout must be finite and non-negative.", call. = FALSE)
  bad_seq <- setdiff(unique(as.character(data$sequence)), c("A_to_B", "B_to_A"))
  if (length(bad_seq)) stop("sequence must contain only A_to_B or B_to_A.", call. = FALSE)
  blocks <- interaction(data$model, data$replicate, data$sequence, drop = TRUE)
  problems <- character()
  for (idx in split(seq_len(nrow(data)), blocks)) {
    z <- data[idx, , drop = FALSE]
    tag <- paste(z$model[1], z$replicate[1], z$sequence[1], sep = "/")
    if (!any(z$dose_A == 0 & z$dose_B == 0)) {
      problems <- c(problems, paste0(tag, ": missing (0,0) vehicle")); next
    }
    for (a in unique(z$dose_A[z$dose_A > 0])) {
      if (!any(z$dose_A == a & z$dose_B == 0)) problems <- c(problems, paste0(tag, ": missing A-only margin for dose_A=", a))
    }
    for (b in unique(z$dose_B[z$dose_B > 0])) {
      if (!any(z$dose_A == 0 & z$dose_B == b)) problems <- c(problems, paste0(tag, ": missing B-only margin for dose_B=", b))
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
collapse_technical_replicates <- function(data, value = "readout", id_cols = intersect(c("model", "replicate", "sequence", "dose_A", "dose_B", "time"), names(data)), fun = mean) {
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
  ans <- do.call(rbind, pieces); rownames(ans) <- NULL; ans
}

#' Normalize sequence-response data
#'
#' @param data A histPharm data frame.
#' @param value Response column.
#' @param blank_col Optional row-wise blank column.
#' @param baseline within_sequence, shared_replicate, or none.
#' @param output Output column name.
#' @param floor Positive numerical floor.
#' @return Data with a normalized response column.
#' @export
normalize_sequence_data <- function(data, value = "readout", blank_col = NULL, baseline = c("within_sequence", "shared_replicate", "none"), output = "response", floor = 1e-8) {
  baseline <- match.arg(baseline)
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", value))
  x <- as.numeric(data[[value]])
  if (!is.null(blank_col)) { .assert_columns(data, blank_col); x <- x - as.numeric(data[[blank_col]]) }
  x <- pmax(x, floor); data$.histpharm_work <- x
  if (baseline == "none") {
    data[[output]] <- x
  } else if (baseline == "within_sequence") {
    grp <- interaction(data$model, data$replicate, data$sequence, drop = TRUE)
    base <- ave(seq_len(nrow(data)), grp, FUN = function(i) {
      z <- data[i, , drop = FALSE]; v <- z$.histpharm_work[z$dose_A == 0 & z$dose_B == 0]
      if (!length(v)) return(rep(NA_real_, length(i))); rep(mean(v), length(i))
    })
    data[[output]] <- x / base
  } else {
    grp <- interaction(data$model, data$replicate, drop = TRUE)
    base <- ave(seq_len(nrow(data)), grp, FUN = function(i) {
      z <- data[i, , drop = FALSE]; v <- z$.histpharm_work[z$dose_A == 0 & z$dose_B == 0]
      if (!length(v)) return(rep(NA_real_, length(i))); rep(mean(v), length(i))
    })
    data[[output]] <- x / base
  }
  data$.histpharm_work <- NULL; data$baseline_method <- baseline; data
}

#' Estimate history-conditioned incremental drug effects
#'
#' @param data Matched sequence data.
#' @param value Response column.
#' @param floor Positive numerical floor.
#' @return Replicate-level history effects.
#' @export
estimate_history_effect <- function(data, value = "response", floor = 1e-8) {
  .assert_columns(data, c("model", "replicate", "sequence", "dose_A", "dose_B", value))
  .group_apply(data, c("model", "replicate", "sequence"), function(z) {
    pairs <- unique(z[z$dose_A > 0 & z$dose_B > 0, c("dose_A", "dose_B"), drop = FALSE])
    if (!nrow(pairs)) return(data.frame())
    out <- lapply(seq_len(nrow(pairs)), function(i) {
      a <- pairs$dose_A[i]; b <- pairs$dose_B[i]
      vals <- pmax(c(.lookup_value(z,0,0,value), .lookup_value(z,a,0,value), .lookup_value(z,0,b,value), .lookup_value(z,a,b,value)), floor)
      v00 <- vals[1]; va0 <- vals[2]; v0b <- vals[3]; vab <- vals[4]
      if (z$sequence[1] == "A_to_B") {
        with_history <- log(va0 / vab); no_history <- log(v00 / v0b); first_drug <- "A"; second_drug <- "B"
      } else {
        with_history <- log(v0b / vab); no_history <- log(v00 / va0); first_drug <- "B"; second_drug <- "A"
      }
      data.frame(model=z$model[1], replicate=z$replicate[1], sequence=z$sequence[1], dose_A=a, dose_B=b, first_drug=first_drug, second_drug=second_drug, second_effect_with_history=with_history, second_effect_without_history=no_history, history_effect=with_history-no_history, V00=v00, Va0=va0, V0b=v0b, Vab=vab, stringsAsFactors=FALSE)
    })
    do.call(rbind, out)
  })
}

#' Estimate reciprocal history asymmetry
#' @param data Reciprocal sequence data.
#' @param value Response column.
#' @param floor Positive numerical floor.
#' @return Replicate-level reciprocal contrasts.
#' @export
estimate_reciprocal_asymmetry <- function(data, value = "response", floor = 1e-8) {
  h <- estimate_history_effect(data, value=value, floor=floor)
  a <- h[h$sequence=="A_to_B", c("model","replicate","dose_A","dose_B","history_effect")]
  b <- h[h$sequence=="B_to_A", c("model","replicate","dose_A","dose_B","history_effect")]
  names(a)[5] <- "H_A_to_B"; names(b)[5] <- "H_B_to_A"
  out <- merge(a,b,by=c("model","replicate","dose_A","dose_B"),all=FALSE)
  out$delta_history <- out$H_A_to_B - out$H_B_to_A; out
}

#' Estimate absolute sequence effect
#' @param data Reciprocal sequence data.
#' @param value Cross-sequence-comparable response column.
#' @param floor Positive numerical floor.
#' @return Replicate-level absolute sequence contrasts.
#' @export
estimate_absolute_sequence_effect <- function(data, value="response", floor=1e-8) {
  .assert_columns(data,c("model","replicate","sequence","dose_A","dose_B",value))
  z <- data[data$dose_A>0 & data$dose_B>0, c("model","replicate","sequence","dose_A","dose_B",value)]
  a <- z[z$sequence=="A_to_B", c("model","replicate","dose_A","dose_B",value)]
  b <- z[z$sequence=="B_to_A", c("model","replicate","dose_A","dose_B",value)]
  names(a)[5] <- "V_A_to_B"; names(b)[5] <- "V_B_to_A"
  out <- merge(a,b,by=c("model","replicate","dose_A","dose_B"),all=FALSE)
  out$absolute_log_ratio <- .safe_log(out$V_B_to_A,floor)-.safe_log(out$V_A_to_B,floor)
  out$preferred_absolute_sequence <- ifelse(out$absolute_log_ratio>0,"A_to_B",ifelse(out$absolute_log_ratio<0,"B_to_A","tie")); out
}

#' Classify joint interaction-efficacy sequence regimes
#' @param reciprocal Reciprocal history table.
#' @param absolute Absolute sequence table.
#' @param tolerance Neutrality threshold.
#' @return Joint sequence-regime table.
#' @export
classify_sequence_regime <- function(reciprocal, absolute, tolerance=0) {
  key <- c("model","replicate","dose_A","dose_B")
  .assert_columns(reciprocal,c(key,"delta_history")); .assert_columns(absolute,c(key,"absolute_log_ratio"))
  z <- merge(reciprocal,absolute,by=key,all=FALSE)
  sh <- ifelse(z$delta_history>tolerance,1L,ifelse(z$delta_history< -tolerance,-1L,0L))
  sa <- ifelse(z$absolute_log_ratio>tolerance,1L,ifelse(z$absolute_log_ratio< -tolerance,-1L,0L))
  z$sequence_regime <- ifelse(sh==1L & sa==1L,"A_to_B_concordant",ifelse(sh==-1L & sa==-1L,"B_to_A_concordant",ifelse(sh==1L & sa==-1L,"A_to_B_interaction_B_to_A_absolute",ifelse(sh==-1L & sa==1L,"B_to_A_interaction_A_to_B_absolute",ifelse(sh==0L & sa==0L,"neutral","partially_neutral")))))
  z
}

#' Transform nominal doses to observed single-agent effect space
#' @param data Sequence-response data.
#' @param value Response column.
#' @param floor Positive numerical floor.
#' @return Effect-space coordinates.
#' @export
transform_effect_space <- function(data,value="response",floor=1e-8) {
  .assert_columns(data,c("model","replicate","sequence","dose_A","dose_B",value))
  .group_apply(data,c("model","replicate","sequence"),function(z) {
    comb <- z[z$dose_A>0 & z$dose_B>0,,drop=FALSE]; if(!nrow(comb)) return(data.frame())
    v00 <- .lookup_value(z,0,0,value)
    ans <- lapply(seq_len(nrow(comb)),function(i) {
      va0 <- .lookup_value(z,comb$dose_A[i],0,value); v0b <- .lookup_value(z,0,comb$dose_B[i],value)
      ar <- va0/pmax(v00,floor); br <- v0b/pmax(v00,floor)
      data.frame(model=comb$model[i],replicate=comb$replicate[i],sequence=comb$sequence[i],dose_A=comb$dose_A[i],dose_B=comb$dose_B[i],marginal_A_residual=ar,marginal_B_residual=br,marginal_A_suppression=1-ar,marginal_B_suppression=1-br,combination_response=comb[[value]][i],stringsAsFactors=FALSE)
    }); do.call(rbind,ans)
  })
}

#' Bootstrap uncertainty for sequence metrics
#' @param metrics Replicate-level metric table.
#' @param metric Numeric metric column.
#' @param strata Surface strata.
#' @param B Bootstrap draws.
#' @param conf Confidence level.
#' @param seed Optional random seed.
#' @return Bootstrap estimates, intervals, and sign probabilities.
#' @export
bootstrap_sequence_metrics <- function(metrics,metric,strata=intersect(c("model","dose_A","dose_B"),names(metrics)),B=2000,conf=0.95,seed=NULL) {
  .assert_columns(metrics,c(strata,"replicate",metric)); if(!is.null(seed)) set.seed(seed); alpha <- (1-conf)/2
  .group_apply(metrics,strata,function(z) {
    x <- z[[metric]][is.finite(z[[metric]])]; if(!length(x)) return(data.frame())
    boots <- replicate(B,mean(sample(x,length(x),replace=TRUE)))
    out <- z[1,strata,drop=FALSE]; out$metric <- metric; out$estimate <- mean(x)
    out$ci_low <- as.numeric(stats::quantile(boots,alpha,names=FALSE)); out$ci_high <- as.numeric(stats::quantile(boots,1-alpha,names=FALSE)); out$prob_gt_zero <- mean(boots>0); out$n_biological_replicates <- length(x); out
  })
}
