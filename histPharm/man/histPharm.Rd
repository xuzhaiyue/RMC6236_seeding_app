\name{histPharm}
\alias{histPharm}
\alias{validate_sequence_data}
\alias{collapse_technical_replicates}
\alias{normalize_sequence_data}
\alias{estimate_history_effect}
\alias{estimate_reciprocal_asymmetry}
\alias{estimate_absolute_sequence_effect}
\alias{classify_sequence_regime}
\alias{transform_effect_space}
\alias{bootstrap_sequence_metrics}
\alias{simulate_sequence_data}
\alias{run_simulation_benchmark}
\alias{compare_state_clock_models}
\alias{evaluate_durable_fate}
\alias{as_dchain_viability_data}
\alias{plot_history_surface}
\alias{plot_sequence_effect_map}
\title{histPharm: history-conditioned pharmacology for sequential drug experiments}
\description{Tools for reciprocal sequence-resolved pharmacology, explicit LLOQ/censoring, strict matched-grid validation, history-conditioned incremental-effect estimation, effect-space mapping, simulation benchmarks, held-out state prediction, and delayed-fate analysis.}
\usage{
validate_sequence_data(data, strict = TRUE, require_complete_grid = TRUE,
  require_timing = TRUE, require_plate = TRUE,
  technical_col = "technical")

collapse_technical_replicates(data, value = "readout",
  id_cols = intersect(c("model", "replicate", "sequence", "dose_A", "dose_B",
  "endpoint", "time", "plate", "calibration_group", "phase1_start",
  "phase1_end", "phase2_start", "phase2_end", "endpoint_time"), names(data)),
  fun = mean)

normalize_sequence_data(data, value = "readout", blank_col = NULL,
  baseline = c("within_sequence", "shared_replicate", "none"),
  output = "response", lloq = NULL, lloq_col = NULL,
  lloq_action = c("mask", "error"))

estimate_history_effect(data, value = "response")

estimate_reciprocal_asymmetry(data, value = "response")

estimate_absolute_sequence_effect(data, value = "response",
  calibration_col = "calibration_group", require_calibration = TRUE)

classify_sequence_regime(reciprocal, absolute, tolerance = 0)

transform_effect_space(data, value = "response")

bootstrap_sequence_metrics(metrics, metric,
  strata = intersect(c("model", "dose_A", "dose_B"), names(metrics)),
  B = 2000, conf = 0.95, seed = NULL, min_replicates = 4)

simulate_sequence_data(dose_A = c(0, 0.25, 0.5, 1, 2),
  dose_B = c(0, 0.25, 0.5, 1, 2), n_models = 3,
  model_names = NULL, n_replicates = 3, n_technical = 2,
  h_A_to_B = 0.5, h_B_to_A = 0, k_A = 0.55, k_B = 0.45,
  a_early_multiplier = 1, a_late_multiplier = 1,
  b_early_multiplier = 1, b_late_multiplier = 1,
  noise_sd = 0.04, plate_sd = noise_sd / 3,
  phase1_start = 0, phase1_end = 48, phase2_start = 48,
  phase2_end = 96, endpoint_time = 96, seed = NULL)

run_simulation_benchmark(n_sim = 100,
  dose_A = c(0, 0.25, 0.5, 1, 2),
  dose_B = c(0, 0.25, 0.5, 1, 2), n_models = 3,
  n_replicates = 3, noise_sd = 0.05, seed = 1)

compare_state_clock_models(data, outcome, model_col = "model",
  replicate_col = "replicate", state_vars, clock_var = "time",
  dose_vars = intersect(c("dose_A", "dose_B"), names(data)))

evaluate_durable_fate(data, acute_time, late_time, value = "response",
  time_col = "time")

as_dchain_viability_data(data, sequence = c("A_to_B", "B_to_A"),
  first_dose, drug_A = "DrugA", drug_B = "DrugB",
  value = "response", plate_col = NULL)

plot_history_surface(metrics, metric = "history_effect", main = NULL)

plot_sequence_effect_map(reciprocal, absolute,
  main = "Sequence effect map")
}
\arguments{
  \item{data}{A long-form data frame. Required canonical columns are described in the README and vignette.}
  \item{strict}{Whether validation problems are errors.}
  \item{require_complete_grid}{Require every dose_A by dose_B coordinate in each reciprocal block.}
  \item{require_timing}{Require phase timing metadata and matched reciprocal schedules.}
  \item{require_plate}{Require explicit plate identifiers.}
  \item{technical_col}{Column containing technical-replicate identifiers.}
  \item{value}{Name of the response column.}
  \item{id_cols}{Columns defining a biological condition for technical-replicate collapse.}
  \item{fun}{Aggregation function for technical repeats.}
  \item{blank_col}{Optional row-wise blank column.}
  \item{baseline}{Normalization reference: within schedule, shared biological replicate, or none.}
  \item{output}{Name of the normalized output column.}
  \item{lloq}{Optional scalar lower limit of quantification after blank correction.}
  \item{lloq_col}{Optional row-wise LLOQ column.}
  \item{lloq_action}{Whether below-LLOQ observations are masked or trigger an error.}
  \item{calibration_col}{Column linking reciprocal schedules to a shared calibration group.}
  \item{require_calibration}{Whether cross-sequence absolute comparison requires explicit calibration metadata.}
  \item{reciprocal}{Output of \code{estimate_reciprocal_asymmetry}.}
  \item{absolute}{Output of \code{estimate_absolute_sequence_effect}.}
  \item{tolerance}{Neutrality threshold around zero.}
  \item{metrics}{A biological-replicate-level metric table.}
  \item{metric}{Name of the metric column.}
  \item{strata}{Columns defining bootstrap strata.}
  \item{B}{Number of bootstrap draws.}
  \item{conf}{Confidence level.}
  \item{seed}{Optional random seed.}
  \item{min_replicates}{Minimum biological replicate count for non-small-n descriptive labelling.}
  \item{dose_A}{Dose grid for drug A, including zero.}
  \item{dose_B}{Dose grid for drug B, including zero.}
  \item{n_models}{Number of biological models when model_names is NULL.}
  \item{model_names}{Optional explicit biological model names.}
  \item{n_replicates}{Number of biological replicates per model.}
  \item{n_technical}{Number of technical repeats per condition.}
  \item{h_A_to_B}{Known-truth A-to-B history effect for simulation.}
  \item{h_B_to_A}{Known-truth B-to-A history effect for simulation.}
  \item{k_A}{Base single-agent effect coefficient for A.}
  \item{k_B}{Base single-agent effect coefficient for B.}
  \item{a_early_multiplier}{Timing multiplier for A when A is first.}
  \item{a_late_multiplier}{Timing multiplier for A when A is second.}
  \item{b_early_multiplier}{Timing multiplier for B when B is first.}
  \item{b_late_multiplier}{Timing multiplier for B when B is second.}
  \item{noise_sd}{Log-normal measurement-noise standard deviation.}
  \item{plate_sd}{Multiplicative plate-to-plate variation in simulation.}
  \item{phase1_start}{Start time of phase 1.}
  \item{phase1_end}{End time of phase 1.}
  \item{phase2_start}{Start time of phase 2.}
  \item{phase2_end}{End time of phase 2.}
  \item{endpoint_time}{Common assay endpoint time.}
  \item{n_sim}{Number of synthetic experiments per benchmark scenario.}
  \item{outcome}{Numeric outcome column for predictive comparison.}
  \item{model_col}{Biological-model column defining held-out folds.}
  \item{replicate_col}{Biological-replicate identifier used in predictions.}
  \item{state_vars}{Measured biological-state predictor columns.}
  \item{clock_var}{Elapsed-time predictor column.}
  \item{dose_vars}{Nominal-dose predictor columns.}
  \item{acute_time}{Acute endpoint time.}
  \item{late_time}{Delayed post-washout endpoint time.}
  \item{time_col}{Time column in hours.}
  \item{sequence}{Sequence to export, either \code{A_to_B} or \code{B_to_A}.}
  \item{first_dose}{Fixed first-drug dose for d-chain-style export.}
  \item{drug_A}{Drug A name.}
  \item{drug_B}{Drug B name.}
  \item{plate_col}{Optional plate identifier column.}
  \item{main}{Optional plot title.}
}
\value{Return values depend on the function. Estimators return biological-replicate-level data frames; simulation and prediction benchmark functions return lists with detailed results and summaries; plotting functions return plotting data invisibly.}
\details{A positive history effect indicates that prior first-drug exposure increased the relative incremental effect of the second drug on the log-multiplicative scale under matched timing controls. It is not, by itself, a claim of mechanism, durable killing, or clinical superiority. Values at or below the LLOQ are not replaced by tiny constants; exact history effects are withheld when a required component is censored. Absolute sequence efficacy is analyzed separately and requires explicit cross-sequence calibration metadata.}
\keyword{methods}
