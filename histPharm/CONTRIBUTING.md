# Contributing to histPharm

Contributions should preserve the distinction between measured observations, formal estimands, and model-assisted interpretation.

1. Add or modify tests for every change to an estimand.
2. Do not smooth or impute raw response surfaces inside primary estimators.
3. Keep biological replicates as the inferential unit; technical replicates are measurement repeats.
4. New benchmark adapters must cite and respect the source method's actual input assumptions.
5. New predictive models must report held-out performance and avoid leakage across the held-out biological model.
