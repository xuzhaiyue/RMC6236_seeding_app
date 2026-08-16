# histPharm 0.1.1

Expert-review remediation release.

- Removed arbitrary `1e-8` flooring from primary normalization/estimators.
- Added explicit LLOQ masking/error handling and censoring status fields.
- Added strict complete-grid, reciprocal-dose-grid, duplicate-condition, timing, and plate validation.
- Required explicit calibration metadata for cross-sequence absolute-burden comparisons.
- Enforced biological-replicate-level bootstrap inputs and small-n exploratory labelling.
- Hardened state-vs-clock and durable-fate functions against technical-well pseudoreplication and invalid time ordering.
- Prevented missing surface coordinates from being silently plotted as zero.
- Added HPAC/SU.86.86 reciprocal 8x8 adversarial validation tests.
- Fixed explicit Author/Maintainer metadata for local R CMD check compatibility.

# histPharm 0.1.0

Initial development release implementing reciprocal history-conditioned pharmacology estimands, simulation benchmarks, state-vs-clock comparison, durable-fate metrics, d-chain schema export, documentation, tests, and CI.
