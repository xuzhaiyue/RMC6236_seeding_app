# histPharm v0.1.1 virtual reciprocal validation
# Models: HPAC and SU.86.86
# Drug A: Gemcitabine-like dose axis
# Drug B: RMC-6236-like dose axis

library(histPharm)

gem <- c(0, 2.5, 5, 10, 15, 25, 50, 100)
rmc <- c(0, 1, 2.5, 5, 10, 15, 30, 60)

x <- simulate_sequence_data(
  dose_A = gem,
  dose_B = rmc,
  model_names = c("HPAC", "SU.86.86"),
  n_replicates = 3,
  n_technical = 2,
  h_A_to_B = 0.45,
  h_B_to_A = 0.10,
  noise_sd = 0.02,
  plate_sd = 0.01,
  phase1_start = 0,
  phase1_end = 48,
  phase2_start = 48,
  phase2_end = 96,
  endpoint_time = 96,
  seed = 8686
)

# Clean-data validation
validate_sequence_data(x)
stopifnot(nrow(x) == 2 * 3 * 2 * 8 * 8 * 2)

collapsed <- collapse_technical_replicates(x)
stopifnot(nrow(collapsed) == 2 * 3 * 2 * 8 * 8)

norm <- normalize_sequence_data(
  collapsed,
  baseline = "within_sequence",
  lloq = 1e-5
)

h <- estimate_history_effect(norm)
dh <- estimate_reciprocal_asymmetry(norm)
abs <- estimate_absolute_sequence_effect(norm)

stopifnot(all(c("HPAC", "SU.86.86") %in% unique(h$model)))
stopifnot(mean(dh$delta_history, na.rm = TRUE) > 0)
stopifnot(all(dh$delta_status == "exact"))
stopifnot(all(abs$absolute_status == "exact"))

# Adversarial 1: below-LLOQ combination must not produce a giant exact H.
lloq_case <- collapsed
hit <- lloq_case$model == "HPAC" & lloq_case$replicate == 1 &
  lloq_case$sequence == "A_to_B" & lloq_case$dose_A == 100 & lloq_case$dose_B == 60
lloq_case$readout[hit] <- 1e-8
lloq_case <- normalize_sequence_data(lloq_case, baseline = "within_sequence", lloq = 0.005)
h_lloq <- estimate_history_effect(lloq_case)
hit_h <- h_lloq$model == "HPAC" & h_lloq$replicate == 1 &
  h_lloq$sequence == "A_to_B" & h_lloq$dose_A == 100 & h_lloq$dose_B == 60
stopifnot(is.na(h_lloq$history_effect[hit_h]))
stopifnot(h_lloq$history_effect_status[hit_h] == "censored_or_missing")

# Adversarial 2: missing positive-positive coordinate must hard fail.
missing_case <- x[!(x$model == "HPAC" & x$replicate == 1 &
  x$sequence == "A_to_B" & x$dose_A == 100 & x$dose_B == 60), ]
stopifnot(inherits(try(validate_sequence_data(missing_case), silent = TRUE), "try-error"))

# Adversarial 3: reciprocal schedule timing mismatch must hard fail.
time_case <- x
time_case$phase2_start[time_case$sequence == "B_to_A"] <- 60
stopifnot(inherits(try(validate_sequence_data(time_case), silent = TRUE), "try-error"))

# Adversarial 4: duplicated collapsed biological condition must hard fail.
dup_case <- rbind(collapsed, collapsed[1, ])
stopifnot(inherits(try(normalize_sequence_data(dup_case), silent = TRUE), "try-error"))

message("histPharm HPAC/SU.86.86 virtual 8x8 validation: PASS")
