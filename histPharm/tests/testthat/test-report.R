testthat::test_that("reciprocal report generates raw-first effect-space outputs", {
  x <- simulate_sequence_data(
    dose_A = c(0, 5, 15, 50),
    dose_B = c(0, 2.5, 10, 30),
    model_names = c("HPAC", "SU.86.86"),
    n_replicates = 3,
    n_technical = 2,
    h_A_to_B = 0.4,
    h_B_to_A = 0.1,
    noise_sd = 0.01,
    plate_sd = 0.005,
    phase1_start = 0,
    phase1_end = 48,
    phase2_start = 48,
    phase2_end = 96,
    endpoint_time = 96,
    seed = 1011
  )

  td <- tempfile("histpharm_report_")
  res <- run_reciprocal_report(x, td, lloq = 1e-6)

  testthat::expect_true(dir.exists(td))
  testthat::expect_true(file.exists(file.path(td, "report_manifest.csv")))
  testthat::expect_true(file.exists(file.path(td, "05_EffectSpace", "CROSS_MODEL_effect_space_overlay_DeltaH.pdf")))
  testthat::expect_true(all(c("HPAC", "SU.86.86") %in% unique(res$effect_space$model)))
  testthat::expect_true(mean(res$reciprocal$delta_history, na.rm = TRUE) > 0)
  testthat::expect_false(anyDuplicated(res$effect_space[c("model", "replicate", "dose_A", "dose_B")]))
})

testthat::test_that("effect-space report does not impute censored history values", {
  x <- simulate_sequence_data(
    dose_A = c(0, 5, 15), dose_B = c(0, 5, 15),
    model_names = "HPAC", n_replicates = 1, n_technical = 2,
    h_A_to_B = 0.3, h_B_to_A = 0.0,
    noise_sd = 0, plate_sd = 0,
    phase1_start = 0, phase1_end = 48,
    phase2_start = 48, phase2_end = 96, endpoint_time = 96,
    seed = 22
  )
  collapsed <- collapse_technical_replicates(x)
  hit <- collapsed$sequence == "A_to_B" & collapsed$dose_A == 15 & collapsed$dose_B == 15
  collapsed$readout[hit] <- 1e-8
  norm <- normalize_sequence_data(collapsed, lloq = 0.001)
  h <- estimate_history_effect(norm)
  target <- h$sequence == "A_to_B" & h$dose_A == 15 & h$dose_B == 15
  testthat::expect_true(is.na(h$history_effect[target]))
})
