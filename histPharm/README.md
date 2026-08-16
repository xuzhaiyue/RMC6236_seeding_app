# histPharm

**histPharm** is an R package for reciprocal, history-conditioned pharmacology in sequential drug-response experiments.

It is designed around a question that static combination scores do not directly estimate:

> Does exposure to a first treatment change the *incremental effect* of the second treatment, and can that history effect be predicted from a measurable biological state?

## Core estimands

For a matched `A_to_B` experiment at doses `a` and `b`, histPharm defines the history-conditioned interaction

`H_AtoB = log[(V_a0 * V_0b) / (V_ab * V_00)]`.

This is equivalently the second-drug effect after A history minus the time-matched second-drug effect without A history. The reciprocal experiment yields `H_BtoA`, and

`DeltaH = H_AtoB - H_BtoA`

quantifies directional history asymmetry.

histPharm intentionally keeps this separate from the absolute sequence contrast

`A_abs = log(V_BtoA / V_AtoB)`,

because a stronger history interaction does not necessarily imply a lower final residual burden.

## What is implemented in v0.1.0

- validation of reciprocal 2D sequence data;
- technical-replicate collapse and explicit normalization choices;
- `H_AtoB`, `H_BtoA`, and reciprocal `DeltaH`;
- absolute residual-burden sequence contrast;
- nominal-dose to single-agent effect-space transformation;
- biological-replicate bootstrap confidence intervals and sign probabilities;
- synthetic benchmark generator with known interaction truth;
- automated no-history, symmetric, and directional simulation benchmark scenarios;
- joint interaction-versus-efficacy sequence-regime classification;
- leave-one-model-out comparison of state, clock, dose, and state+clock predictors;
- post-washout regrowth / durable-fate metrics;
- heatmap-style base-R visualizations;
- export of fixed-first sequential curves to the column schema used by the published d-chain example dataset.

## Input schema

A minimal long-form data frame contains:

```text
model  replicate  sequence  dose_A  dose_B  readout
```

where `sequence` is exactly `A_to_B` or `B_to_A`. Each biological replicate and sequence must contain the matched `(0,0)` vehicle, both single-drug margins, and the combination points to be estimated.

Technical replicates can be kept as individual rows and collapsed explicitly with `collapse_technical_replicates()`.

## Minimal example

```r
library(histPharm)

x <- simulate_sequence_data(
  dose_A = c(0, 0.25, 0.5, 1, 2),
  dose_B = c(0, 0.25, 0.5, 1, 2),
  n_models = 3,
  n_replicates = 3,
  h_A_to_B = 0.55,
  h_B_to_A = 0.10,
  seed = 1
)

x <- collapse_technical_replicates(x)
x <- normalize_sequence_data(x, baseline = "shared_replicate")

h <- estimate_history_effect(x)
dh <- estimate_reciprocal_asymmetry(x)
abs <- estimate_absolute_sequence_effect(x)

head(h)
head(dh)
head(abs)
```

## Interpretation discipline

histPharm does **not** rename every positive interaction as biological synergy. `history_effect > 0` means that, on the chosen log-multiplicative scale and under matched timing controls, prior treatment increased the relative incremental effect of the second drug. Biological mechanism, durable killing, and clinical superiority require separate evidence.

Raw observations determine conclusions. Smoothing, curve fitting, and interpolation are visualization or prediction aids and must not replace measured values.

## Benchmarking philosophy

The package is intended to be compared with, not rhetorically substituted for, established frameworks. d-chain is an important published Bayesian sequential dose-response benchmark. histPharm adds estimands enabled by a reciprocal 2D design: direct history-effect asymmetry, separation of interaction from absolute efficacy, experimentally measured state prediction, and delayed regenerative fate.

## Installation from the current development location

```r
remotes::install_github(
  "xuzhaiyue/RMC6236_seeding_app",
  ref = "agent/histpharm-method-framework",
  subdir = "histPharm"
)
```

The package is being developed on an isolated branch before it is split into its own `histPharm` repository.
