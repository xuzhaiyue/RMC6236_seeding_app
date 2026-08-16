#' Simulate reciprocal sequence-response data with known history effects
#'
#' The simulator generates complete reciprocal dose surfaces whose history
#' effects are known by construction. It is intended for estimator validation,
#' power exploration, and benchmark scenarios rather than biological inference.
#'
#' @param dose_A Numeric dose grid for drug A, including zero.
#' @param dose_B Numeric dose grid for drug B, including zero.
#' @param n_models Number of biological models.
#' @param n_replicates Number of biological replicates per model.
#' @param n_technical Number of technical repeats per condition.
#' @param h_A_to_B True A-to-B history effect.
#' @param h_B_to_A True B-to-A history effect.
#' @param k_A Base single-agent effect coefficient for A.
#' @param k_B Base single-agent effect coefficient for B.
#' @param a_early_multiplier Timing multiplier for A when given first.
#' @param a_late_multiplier Timing multiplier for A when given second.
#' @param b_early_multiplier Timing multiplier for B when given first.
#' @param b_late_multiplier Timing multiplier for B when given second.
#' @param noise_sd Log-normal measurement-noise standard deviation.
#' @param seed Optional random seed.
#' @return A long-form synthetic histPharm data frame.
#' @export
simulate_sequence_data <- function(
  dose_A = c(0, 0.25, 0.5, 1, 2),
  dose_B = c(0, 0.25, 0.5, 1, 2),
  n_models = 3,
  n_replicates = 3,
  n_technical = 2,
  h_A_to_B = 0.5,
  h_B_to_A = 0,
  k_A = 0.55,
  k_B = 0.45,
  a_early_multiplier = 1,
  a_late_multiplier = 1,
  b_early_multiplier = 1,
  b_late_multiplier = 1,
  noise_sd = 0.04,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (!0 %in% dose_A || !0 %in% dose_B) stop("Both dose grids must include zero.", call. = FALSE)
  if (n_models < 1 || n_replicates < 1 || n_technical < 1) stop("Replication counts must be positive.", call. = FALSE)

  grid <- expand.grid(
    model = paste0("M", seq_len(n_models)),
    replicate = seq_len(n_replicates),
    sequence = c("A_to_B", "B_to_A"),
    dose_A = dose_A,
    dose_B = dose_B,
    technical = seq_len(n_technical),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  model_A <- stats::setNames(stats::rlnorm(n_models, 0, 0.12), paste0("M", seq_len(n_models)))
  model_B <- stats::setNames(stats::rlnorm(n_models, 0, 0.12), paste0("M", seq_len(n_models)))
  plate_keys <- unique(grid[c("model", "replicate", "sequence")])
  plate_keys$plate_factor <- stats::rlnorm(nrow(plate_keys), 0, noise_sd / 3)
  grid <- merge(grid, plate_keys, by = c("model", "replicate", "sequence"), sort = FALSE)

  out <- numeric(nrow(grid))
  truth <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    z <- grid[i, ]
    a <- z$dose_A
    b <- z$dose_B
    if (z$sequence == "A_to_B") {
      eff_A <- exp(-k_A * model_A[z$model] * a * a_early_multiplier)
      eff_B <- exp(-k_B * model_B[z$model] * b * b_late_multiplier)
      h <- if (a > 0 && b > 0) h_A_to_B else 0
    } else {
      eff_A <- exp(-k_A * model_A[z$model] * a * a_late_multiplier)
      eff_B <- exp(-k_B * model_B[z$model] * b * b_early_multiplier)
      h <- if (a > 0 && b > 0) h_B_to_A else 0
    }
    mu <- z$plate_factor * eff_A * eff_B * exp(-h)
    out[i] <- mu * stats::rlnorm(1, 0, noise_sd)
    truth[i] <- h
  }
  grid$readout <- out
  grid$true_history_effect <- truth
  grid$plate_factor <- NULL
  grid
}

#' Run known-truth simulation benchmarks
#'
#' Benchmarks no-history, symmetric-history, and directional-history scenarios.
#' The summary reports bias, RMSE, and direction recovery for reciprocal
#' history asymmetry.
#'
#' @param n_sim Number of synthetic experiments per scenario.
#' @param dose_A Dose grid for A.
#' @param dose_B Dose grid for B.
#' @param n_models Number of biological models per synthetic experiment.
#' @param n_replicates Biological replicates per model.
#' @param noise_sd Log-normal noise standard deviation.
#' @param seed Random seed.
#' @return A list with simulation-level estimates and a scenario summary.
#' @export
run_simulation_benchmark <- function(
  n_sim = 100,
  dose_A = c(0, 0.25, 0.5, 1, 2),
  dose_B = c(0, 0.25, 0.5, 1, 2),
  n_models = 3,
  n_replicates = 3,
  noise_sd = 0.05,
  seed = 1
) {
  scenarios <- data.frame(
    scenario = c("no_history", "symmetric_history", "A_directional", "B_directional"),
    h_A_to_B = c(0, 0.4, 0.6, 0.1),
    h_B_to_A = c(0, 0.4, 0.1, 0.6),
    stringsAsFactors = FALSE
  )
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_sim * nrow(scenarios))
  cursor <- 1L
  res <- list()
  for (s in seq_len(nrow(scenarios))) {
    truth_delta <- scenarios$h_A_to_B[s] - scenarios$h_B_to_A[s]
    for (j in seq_len(n_sim)) {
      x <- simulate_sequence_data(
        dose_A = dose_A, dose_B = dose_B,
        n_models = n_models, n_replicates = n_replicates,
        n_technical = 2,
        h_A_to_B = scenarios$h_A_to_B[s],
        h_B_to_A = scenarios$h_B_to_A[s],
        noise_sd = noise_sd,
        seed = seeds[cursor]
      )
      cursor <- cursor + 1L
      x <- collapse_technical_replicates(x)
      x <- normalize_sequence_data(x, baseline = "shared_replicate")
      d <- estimate_reciprocal_asymmetry(x)
      res[[length(res) + 1L]] <- data.frame(
        scenario = scenarios$scenario[s],
        simulation = j,
        truth_delta = truth_delta,
        estimate_delta = mean(d$delta_history, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  sims <- do.call(rbind, res)
  summary <- .group_apply(sims, "scenario", function(z) {
    truth <- z$truth_delta[1]
    err <- z$estimate_delta - truth
    if (abs(truth) < 1e-12) {
      direction_recovery <- mean(abs(z$estimate_delta) < 0.15)
    } else {
      direction_recovery <- mean(sign(z$estimate_delta) == sign(truth))
    }
    data.frame(
      scenario = z$scenario[1],
      truth_delta = truth,
      mean_estimate = mean(z$estimate_delta),
      bias = mean(err),
      rmse = sqrt(mean(err^2)),
      direction_recovery = direction_recovery,
      stringsAsFactors = FALSE
    )
  })
  list(simulations = sims, summary = summary)
}
