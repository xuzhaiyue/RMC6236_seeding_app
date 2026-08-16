test_that("known history effects are recovered exactly without noise", {
  x <- simulate_sequence_data(
    dose_A = c(0, 1), dose_B = c(0, 1),
    n_models = 2, n_replicates = 2, n_technical = 2,
    h_A_to_B = 0.6, h_B_to_A = 0.1,
    noise_sd = 0, seed = 1
  )
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "shared_replicate")
  h <- estimate_history_effect(x)
  expect_equal(h$history_effect[h$sequence == "A_to_B"], rep(0.6, 4), tolerance = 1e-10)
  expect_equal(h$history_effect[h$sequence == "B_to_A"], rep(0.1, 4), tolerance = 1e-10)
  d <- estimate_reciprocal_asymmetry(x)
  expect_equal(d$delta_history, rep(0.5, 4), tolerance = 1e-10)
})

test_that("interaction asymmetry and absolute efficacy remain distinct", {
  x <- simulate_sequence_data(
    dose_A = c(0, 1), dose_B = c(0, 1),
    n_models = 1, n_replicates = 2, n_technical = 1,
    h_A_to_B = 0.3, h_B_to_A = 0.3,
    a_early_multiplier = 1.5, a_late_multiplier = 0.5,
    noise_sd = 0, seed = 2
  )
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "shared_replicate")
  d <- estimate_reciprocal_asymmetry(x)
  a <- estimate_absolute_sequence_effect(x)
  expect_equal(d$delta_history, rep(0, nrow(d)), tolerance = 1e-10)
  expect_true(any(abs(a$absolute_log_ratio) > 1e-6))
})

test_that("validation catches missing matched vehicle", {
  x <- data.frame(
    model = "M1", replicate = 1, sequence = "A_to_B",
    dose_A = 1, dose_B = 1, readout = 0.5
  )
  expect_error(validate_sequence_data(x), "missing \\(0,0\\) vehicle")
})

test_that("technical replicates collapse without inflating biological n", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), n_models = 1, n_replicates = 1, n_technical = 3, noise_sd = 0, seed = 3)
  y <- collapse_technical_replicates(x)
  expect_true(all(y$technical_n == 3))
  expect_equal(nrow(y), 8)
})

test_that("effect-space coordinates are finite", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), n_models = 1, n_replicates = 1, n_technical = 1, noise_sd = 0, seed = 4)
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "shared_replicate")
  e <- transform_effect_space(x)
  expect_true(all(is.finite(e$marginal_A_suppression)))
  expect_true(all(is.finite(e$marginal_B_suppression)))
})

test_that("bootstrap returns uncertainty without treating technical wells as n", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), n_models = 1, n_replicates = 4, n_technical = 3, noise_sd = 0.02, seed = 5)
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "shared_replicate")
  d <- estimate_reciprocal_asymmetry(x)
  b <- bootstrap_sequence_metrics(d, "delta_history", B = 50, seed = 1)
  expect_equal(b$n_biological_replicates, 4)
  expect_true(all(b$ci_low <= b$estimate & b$estimate <= b$ci_high))
})

test_that("durable fate reports regrowth separately", {
  x <- data.frame(
    model = rep("M1", 2), replicate = 1, sequence = "A_to_B",
    dose_A = 1, dose_B = 1, time = c(96, 144), response = c(0.2, 0.4)
  )
  z <- evaluate_durable_fate(x, 96, 144)
  expect_equal(z$fold_regrowth, 2)
  expect_equal(z$log_growth_per_hour, log(2) / 48)
})

test_that("state versus clock uses held-out biological models", {
  set.seed(1)
  x <- expand.grid(model = paste0("M", 1:4), time = c(0, 12, 24, 36), rep = 1:2)
  x$state1 <- x$time / 12 + as.numeric(factor(x$model)) * 0.1
  x$dose_A <- 1
  x$dose_B <- 1
  x$outcome <- 0.8 * x$state1 + rnorm(nrow(x), sd = 0.05)
  z <- compare_state_clock_models(x, outcome = "outcome", model_col = "model", state_vars = "state1", clock_var = "time")
  expect_true(all(c("clock", "dose", "state", "state_clock") %in% z$summary$predictor))
  expect_equal(length(unique(z$predictions$held_out_model)), 4)
})

test_that("simulation benchmark distinguishes directional truth", {
  z <- run_simulation_benchmark(n_sim = 3, dose_A = c(0, 1), dose_B = c(0, 1), n_models = 2, n_replicates = 2, noise_sd = 0.01, seed = 6)
  a <- z$summary[z$summary$scenario == "A_directional", ]
  b <- z$summary[z$summary$scenario == "B_directional", ]
  expect_gt(a$mean_estimate, 0)
  expect_lt(b$mean_estimate, 0)
})

test_that("d-chain adapter emits the published example column schema", {
  x <- simulate_sequence_data(dose_A = c(0, 1), dose_B = c(0, 1), n_models = 1, n_replicates = 1, n_technical = 1, noise_sd = 0, seed = 7)
  x <- collapse_technical_replicates(x)
  x <- normalize_sequence_data(x, baseline = "shared_replicate")
  z <- as_dchain_viability_data(x, sequence = "A_to_B", first_dose = 1, drug_A = "Gem", drug_B = "RMC")
  expect_equal(names(z), c("Experiment", "CellLine", "Run", "Plate", "Pretreatment", "Compound", "Concentration", "RelCount"))
  expect_true(all(c("AB", "A0", "A") %in% unique(z$Experiment)))
})
