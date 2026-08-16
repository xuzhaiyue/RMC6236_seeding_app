test_that("known history effects are recovered exactly without noise", {
  x <- simulate_sequence_data(
    dose_A = c(0, 1), dose_B = c(0, 1),
    model_names = c("HPAC", "SU.86.86"), n_replicates = 2, n_technical = 2,
    h_A_to_B = 0.6, h_B_to_A = 0.1,
    noise_sd = 0, plate_sd = 0, seed = 1
  )
  expect_silent(validate_sequence_data(x))
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "within_sequence")
  h <- estimate_history_effect(x)
  expect_equal(h$history_effect[h$sequence == "A_to_B"], rep(0.6, 4), tolerance = 1e-10)
  expect_equal(h$history_effect[h$sequence == "B_to_A"], rep(0.1, 4), tolerance = 1e-10)
  d <- estimate_reciprocal_asymmetry(x)
  expect_equal(d$delta_history, rep(0.5, 4), tolerance = 1e-10)
})

test_that("interaction asymmetry and absolute efficacy remain distinct", {
  x <- simulate_sequence_data(
    dose_A = c(0, 1), dose_B = c(0, 1),
    model_names = "HPAC", n_replicates = 2, n_technical = 1,
    h_A_to_B = 0.3, h_B_to_A = 0.3,
    a_early_multiplier = 1.5, a_late_multiplier = 0.5,
    noise_sd = 0, plate_sd = 0, seed = 2
  )
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "within_sequence")
  d <- estimate_reciprocal_asymmetry(x)
  a <- estimate_absolute_sequence_effect(x)
  expect_equal(d$delta_history, rep(0, nrow(d)), tolerance = 1e-10)
  expect_true(any(abs(a$absolute_log_ratio) > 1e-6))
})

test_that("validation catches missing matched vehicle", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0)
  x <- x[!(x$sequence == "A_to_B" & x$dose_A == 0 & x$dose_B == 0), ]
  expect_error(validate_sequence_data(x), "missing \\(0,0\\) vehicle")
})

test_that("technical replicates collapse without inflating biological n", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 3, noise_sd = 0, plate_sd = 0, seed = 3)
  expect_silent(validate_sequence_data(x))
  y <- collapse_technical_replicates(x)
  expect_true(all(y$technical_n == 3))
  expect_equal(nrow(y), 8)
})

test_that("duplicate collapsed biological conditions hard fail", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0)
  y <- collapse_technical_replicates(x)
  y <- rbind(y, y[1, ])
  expect_error(normalize_sequence_data(y), "Duplicate biological conditions")
  expect_error(estimate_history_effect(transform(y, response = readout)), "Duplicate biological conditions")
})

test_that("missing positive-positive coordinate hard fails complete-grid validation", {
  x <- simulate_sequence_data(dose_A = c(0, 1, 2), dose_B = c(0, 1, 2), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0)
  x <- x[!(x$sequence == "A_to_B" & x$dose_A == 2 & x$dose_B == 2), ]
  expect_error(validate_sequence_data(x), "incomplete dose grid")
})

test_that("reciprocal timing mismatch hard fails", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0)
  x$phase2_start[x$sequence == "B_to_A"] <- 60
  expect_error(validate_sequence_data(x), "not time matched|invalid phase timing order")
})

test_that("below-LLOQ values are censored rather than floored", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, h_A_to_B = 0.5, h_B_to_A = 0, noise_sd = 0, plate_sd = 0)
  y <- collapse_technical_replicates(x)
  hit <- y$sequence == "A_to_B" & y$dose_A == 1 & y$dose_B == 1
  y$readout[hit] <- 1e-8
  y <- normalize_sequence_data(y, baseline = "within_sequence", lloq = 0.01)
  h <- estimate_history_effect(y)
  z <- h[h$sequence == "A_to_B", ]
  expect_true(is.na(z$history_effect))
  expect_equal(z$history_effect_status, "censored_or_missing")
  expect_gt(z$n_censored_components, 0)
})

test_that("non-positive blank-corrected values require explicit LLOQ", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0)
  y <- collapse_technical_replicates(x)
  y$blank <- 0
  y$blank[1] <- y$readout[1] + 0.1
  expect_error(normalize_sequence_data(y, blank_col = "blank"), "Provide an explicit LLOQ")
})

test_that("absolute sequence effect requires calibration metadata", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0)
  y <- collapse_technical_replicates(x)
  y <- normalize_sequence_data(y, baseline = "within_sequence")
  z <- y
  z$calibration_group <- NULL
  expect_error(estimate_absolute_sequence_effect(z), "calibration_group")
  expect_silent(estimate_absolute_sequence_effect(y))
})

test_that("effect-space coordinates are finite for quantified margins", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0, seed = 4)
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "within_sequence")
  e <- transform_effect_space(x)
  expect_true(all(is.finite(e$marginal_A_suppression)))
  expect_true(all(is.finite(e$marginal_B_suppression)))
})

test_that("bootstrap is biological-replicate-level and labels n=3 exploratory", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 3, n_technical = 3, noise_sd = 0.02, plate_sd = 0, seed = 5)
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "within_sequence")
  d <- estimate_reciprocal_asymmetry(x)
  b <- bootstrap_sequence_metrics(d, "delta_history", B = 50, seed = 1)
  expect_equal(b$n_biological_replicates, 3)
  expect_true(all(b$inference_grade == "exploratory_small_n"))
})

test_that("four biological replicates receive descriptive bootstrap grade", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 4, n_technical = 2, noise_sd = 0.02, plate_sd = 0, seed = 51)
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "within_sequence")
  d <- estimate_reciprocal_asymmetry(x)
  b <- bootstrap_sequence_metrics(d, "delta_history", B = 50, seed = 1)
  expect_true(all(b$inference_grade == "bootstrap_descriptive"))
})

test_that("durable fate reports response change separately and validates time order", {
  x <- data.frame(
    model = rep("HPAC", 2), replicate = 1, sequence = "A_to_B",
    dose_A = 1, dose_B = 1, time = c(96, 144), response = c(0.2, 0.4)
  )
  z <- evaluate_durable_fate(x, 96, 144)
  expect_equal(z$fold_response_change, 2)
  expect_equal(z$log_response_change_per_hour, log(2) / 48)
  expect_error(evaluate_durable_fate(x, 144, 96), "strictly greater")
})

test_that("state versus clock rejects technical-well pseudoreplication", {
  set.seed(1)
  x <- expand.grid(model = paste0("M", 1:4), time = c(0, 12, 24, 36), replicate = 1:2)
  x$state1 <- x$time / 12 + as.numeric(factor(x$model)) * 0.1
  x$dose_A <- 1
  x$dose_B <- 1
  x$outcome <- 0.8 * x$state1 + rnorm(nrow(x), sd = 0.05)
  z <- compare_state_clock_models(x, outcome = "outcome", model_col = "model", state_vars = "state1", clock_var = "time")
  expect_true(all(c("clock", "dose", "state", "state_clock") %in% z$summary$predictor))
  expect_equal(length(unique(z$predictions$held_out_model)), 4)
  x$technical <- 1
  expect_error(compare_state_clock_models(x, outcome = "outcome", model_col = "model", state_vars = "state1"), "collapse technical")
})

test_that("simulation benchmark distinguishes directional truth", {
  z <- run_simulation_benchmark(n_sim = 3, dose_A = c(0, 1), dose_B = c(0, 1), n_models = 2, n_replicates = 2, noise_sd = 0.01, seed = 6)
  a <- z$summary[z$summary$scenario == "A_directional", ]
  b <- z$summary[z$summary$scenario == "B_directional", ]
  expect_gt(a$mean_estimate, 0)
  expect_lt(b$mean_estimate, 0)
})

test_that("d-chain adapter emits the published example column schema", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), model_names = "HPAC", n_replicates = 1, n_technical = 1, noise_sd = 0, plate_sd = 0, seed = 7)
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "within_sequence")
  z <- as_dchain_viability_data(x, sequence = "A_to_B", first_dose = 1, drug_A = "Gem", drug_B = "RMC")
  expect_equal(names(z), c("Experiment", "CellLine", "Run", "Plate", "Pretreatment", "Compound", "Concentration", "RelCount"))
  expect_true(all(c("AB", "A0", "A") %in% unique(z$Experiment)))
})

test_that("virtual HPAC and SU.86.86 reciprocal 8x8 data pass the full calibration pipeline", {
  gem <- c(0, 2.5, 5, 10, 15, 25, 50, 100)
  rmc <- c(0, 1, 2.5, 5, 10, 15, 30, 60)
  x <- simulate_sequence_data(
    dose_A = gem, dose_B = rmc,
    model_names = c("HPAC", "SU.86.86"),
    n_replicates = 3, n_technical = 2,
    h_A_to_B = 0.45, h_B_to_A = 0.10,
    noise_sd = 0.02, plate_sd = 0.01, seed = 8686
  )
  expect_equal(nrow(x), 2 * 3 * 2 * 8 * 8 * 2)
  expect_silent(validate_sequence_data(x))
  y <- collapse_technical_replicates(x)
  expect_equal(nrow(y), 2 * 3 * 2 * 8 * 8)
  y <- normalize_sequence_data(y, baseline = "within_sequence", lloq = 1e-5)
  h <- estimate_history_effect(y)
  d <- estimate_reciprocal_asymmetry(y)
  expect_true(all(c("HPAC", "SU.86.86") %in% unique(h$model)))
  expect_gt(mean(d$delta_history, na.rm = TRUE), 0)
  expect_true(all(d$delta_status == "exact"))
})
