# histPharm

**histPharm** is an R package for reciprocal, history-conditioned pharmacology in sequential drug-response experiments.

It is designed around a question that static combination scores do not directly estimate:

> Does exposure to a first treatment change the *incremental effect* of the second treatment, and can that history effect be predicted from a measurable biological state?

## Core estimands

For a matched `A_to_B` experiment at doses `a` and `b`, histPharm defines

`H_AtoB = log[(V_a0 * V_0b) / (V_ab * V_00)]`.

This is equivalently the second-drug effect after A history minus the time-matched second-drug effect without A history. The reciprocal experiment yields `H_BtoA`, and

`DeltaH = H_AtoB - H_BtoA`

quantifies directional history asymmetry.

histPharm intentionally keeps this separate from the cross-sequence burden contrast

`A_abs = log(V_BtoA / V_AtoB)`,

because a stronger history interaction does not necessarily imply a lower final residual burden.

## v0.1.1 expert-review remediation

v0.1.1 hardens the calibration pipeline after methodological review:

- **No arbitrary tiny-value flooring.** Values at/below an explicit LLOQ are masked/censored; exact `H` is not reported when a required component is censored.
- **Strict reciprocal-grid validation.** Missing positive-positive coordinates, missing margins, unmatched reciprocal dose grids, or incomplete reciprocal schedules fail validation.
- **Technical versus biological replication is enforced.** Raw technical replicates require explicit technical identifiers; collapsed biological conditions must be unique.
- **Timing metadata is enforced.** Publication-grade validation requires `phase1_start`, `phase1_end`, `phase2_start`, `phase2_end`, and `endpoint_time`, with reciprocal schedules time matched.
- **Plate metadata is enforced.** Raw reciprocal blocks require plate identifiers; cross-sequence burden comparisons require explicit `calibration_group` metadata.
- **Small-n bootstrap is labelled exploratory.** The package does not reinterpret bootstrap sign frequency as a p-value or posterior probability.
- **Durable-fate terminology is restricted.** Per-hour change is reported as `log_response_change_per_hour`; it is a growth rate only when the assay is validated as proportional to cell number.
- **Missing surface coordinates remain missing.** Plotting no longer silently renders absent grid cells as zero.

## Publication-grade input schema

The recommended long-form input contains:

```text
model
replicate
sequence
dose_A
dose_B
technical
plate
calibration_group
phase1_start
phase1_end
phase2_start
phase2_end
endpoint_time
endpoint
readout
```

where `sequence` is exactly `A_to_B` or `B_to_A`.

Each biological replicate and sequence must contain the complete matched dose grid, including `(0,0)`, A-only margins, B-only margins, and all positive-positive coordinates. Reciprocal sequences must use the same dose grid and phase timing.

## Safe calibration workflow

```r
library(histPharm)

x <- simulate_sequence_data(
  dose_A = c(0, 2.5, 5, 10, 15, 25, 50, 100),
  dose_B = c(0, 1, 2.5, 5, 10, 15, 30, 60),
  model_names = c("HPAC", "SU.86.86"),
  n_replicates = 3,
  n_technical = 2,
  h_A_to_B = 0.45,
  h_B_to_A = 0.10,
  seed = 8686
)

validate_sequence_data(x)

x <- collapse_technical_replicates(x)
x <- normalize_sequence_data(
  x,
  baseline = "within_sequence",
  lloq = 0.005
)

h <- estimate_history_effect(x)
dh <- estimate_reciprocal_asymmetry(x)
abs <- estimate_absolute_sequence_effect(x)

head(h)
head(dh)
head(abs)
```

## One-command reciprocal discovery report

`run_reciprocal_report()` turns a complete raw reciprocal matrix into a fixed, auditable report package. It first validates the experimental design, then collapses technical replicates, performs explicit LLOQ masking, calculates the core estimands, writes machine-readable tables, and exports PDF figures.

```r
report <- run_reciprocal_report(
  data = raw_data,
  outdir = "histPharm_report",
  lloq = YOUR_ASSAY_LLOQ
)
```

The output tree is fixed:

```text
histPharm_report/
  01_QC/
    collapsed_biological_conditions.csv
    normalized_response_with_LLOQ_flags.csv
  02_Response/
    <model>_single_agent_curves.pdf
    <model>_A_to_B_response_surface.pdf
    <model>_B_to_A_response_surface.pdf
  03_History/
    history_effects.csv
    <model>_A_to_B_H_surface.pdf
    <model>_B_to_A_H_surface.pdf
  04_Sequence/
    reciprocal_deltaH.csv
    absolute_sequence_contrast.csv   # only when calibration is valid
    <model>_DeltaH_surface.pdf
    <model>_absolute_sequence_surface.pdf
    <model>_interaction_vs_absolute.pdf
  05_EffectSpace/
    effect_space_projection.csv
    <model>_effect_space_DeltaH.pdf
    CROSS_MODEL_effect_space_overlay_DeltaH.pdf
  report_manifest.csv
```

### Effect-space output is a primary visualization guardrail

The cross-model effect-space overlay maps physical doses to observed single-agent residual fractions. It asks whether models that require different nominal concentrations nevertheless occupy similar biological-response coordinates when history asymmetry is strong.

The plot does **not** predefine a 50--70% "golden center" and does **not** infer an overlap island by smoothing. The 0.5 and 0.7 reference lines are visual landmarks only. All displayed locations are observed marginal-response coordinates, and point size encodes the magnitude of `DeltaH`.

This supports the intended claim only if the raw data actually show convergence:

> sequence-specific vulnerability may be organized by biological perturbation depth rather than absolute nominal dose.

If HPAC and SU.86.86 do not converge in effect space, the plot remains a valid negative or model-specific result.

## LLOQ discipline

histPharm does not convert an observation that is merely known to be below quantification into an invented value such as `1e-8`.

If a required component of the four-cell history contrast is censored, the exact history effect is returned as `NA` with `history_effect_status = "censored_or_missing"`.

This prevents a below-LLOQ combination response from generating an artificial very large history effect.

## Interpretation discipline

histPharm does **not** rename every positive interaction as biological synergy. `history_effect > 0` means that, on the chosen log-multiplicative scale and under matched timing controls, prior treatment increased the relative incremental effect of the second drug.

It does not by itself establish molecular memory, a specific mechanism, irreversible killing, or clinical superiority.

Raw observations determine conclusions. Smoothing, curve fitting, and interpolation are visualization or prediction aids and must not replace measured values.

## Benchmarking philosophy

The package is intended to be compared with, not rhetorically substituted for, established frameworks. d-chain is an important published Bayesian sequential dose-response benchmark. histPharm adds estimands enabled by a reciprocal 2D design: direct history-effect asymmetry, separation of interaction from cross-sequence burden, experimentally measured state prediction, and delayed regenerative fate.

## Installation from the current development location

```r
remotes::install_github(
  "xuzhaiyue/RMC6236_seeding_app",
  ref = "agent/histpharm-method-framework",
  subdir = "histPharm"
)
```

The package remains on an isolated development branch until expert-review remediation and real-data validation are complete.
