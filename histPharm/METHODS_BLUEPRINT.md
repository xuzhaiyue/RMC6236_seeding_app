# histPharm methods blueprint

## Methodological target

histPharm is not positioned as the first framework for sequential therapy, state-guided therapy, or schedule-dependent combination modeling. Those fields have established precedents.

The specific contribution targeted here is narrower and testable:

> An experimentally anchored reciprocal framework that estimates treatment-history-conditioned second-drug effects, separates interaction asymmetry from absolute efficacy, maps nominal dose to observed perturbation depth, and tests whether measured recovery state predicts sequence-specific vulnerability and delayed fate better than elapsed clock time.

## Primary estimands

For matched response `V` within the `A_to_B` schedule:

`H_AtoB(a,b) = log[(V_a0 * V_0b) / (V_ab * V_00)]`.

This is exactly:

`log(V_a0 / V_ab) - log(V_00 / V_0b)`,

the second-drug effect after A history minus the time-matched second-drug effect without A history.

The reciprocal surface defines `H_BtoA`; directional history asymmetry is

`DeltaH = H_AtoB - H_BtoA`.

Absolute efficacy is separate:

`A_abs = log(V_BtoA / V_AtoB)`.

Positive `A_abs` means the A-to-B sequence leaves lower measured burden at the same nominal pair.

## Non-negotiable design assumptions

1. The reciprocal schedules must be genuinely time matched.
2. Each schedule must retain the `(0,0)`, A-only, and B-only margins required for the counterfactual contrast.
3. Biological replicates, not technical wells, are the inferential unit.
4. Raw observations determine the primary conclusion. Smoothing is not allowed to replace data.
5. Absolute sequence comparisons require a normalization that remains interpretable across reciprocal plates/schedules.
6. A positive history effect is not automatically called synergy, mechanism, or durable killing.

## Validation ladder

### Layer 1 — Mathematical recovery
Use known-truth simulations to quantify bias, RMSE, sign recovery, false directional calls, sensitivity to measurement noise, plate drift, and timing-specific single-agent effects.

### Layer 2 — Published benchmarks
Use simultaneous-combination scores only as conventional comparators and d-chain as an established sequential dose-response comparator. The goal is not to claim universal superiority; the goal is to show which estimand each framework can and cannot recover.

### Layer 3 — Reciprocal experimental surfaces
Prospectively analyze complete A-to-B and B-to-A dose surfaces with frozen preprocessing and prespecified margins.

### Layer 4 — Effect-space generalization
Test whether history asymmetry converges across biological models when nominal doses are replaced by observed single-agent perturbation depth.

### Layer 5 — State versus clock
Measure recovery state independently of the response outcome. Compare clock-only, dose-only, state-only, and state+clock models with leave-one-biological-model-out validation.

### Layer 6 — Prospective prediction
Freeze the state model, predict an intervention window in a held-out model, and test that window experimentally without refitting first.

### Layer 7 — Delayed fate
Analyze acute burden and post-washout regenerative fitness as separate outcomes. Do not collapse them into an arbitrary weighted master score.

### Layer 8 — Cross-drug generalization
Develop on one drug pair and test on additional mechanistically distinct first-drug classes. Negative or reversed sequence effects are informative if correctly identified.

## Claim boundary

A complete histPharm methods claim requires more than a new heatmap. The minimum package-level method contribution is:

- formal estimands;
- reciprocal design requirements;
- uncertainty quantification;
- synthetic truth benchmarks;
- comparison with published methods;
- held-out biological-state prediction;
- at least one prospective prediction;
- delayed-fate analysis;
- cross-drug or cross-model generalization;
- open-source implementation and frozen analysis examples.
