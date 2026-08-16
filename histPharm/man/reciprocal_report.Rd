\name{run_reciprocal_report}
\alias{run_reciprocal_report}
\alias{plot_response_surface}
\alias{plot_metric_surface}
\alias{build_effect_space_projection}
\alias{plot_effect_space_projection}
\alias{plot_cross_model_effect_space}
\title{Raw-first reciprocal discovery report and effect-space visualization}
\description{
Generate publication-oriented reciprocal sequence-response reports without smoothing,
interpolation, or silent imputation. Effect-space plots use observed single-agent
marginal residual fractions and do not assume that a common high-effect region exists.
}
\usage{
run_reciprocal_report(data, outdir, value = "readout", blank_col = NULL,
  lloq = NULL, lloq_col = NULL, baseline = "within_sequence",
  require_absolute = TRUE, calibration_col = "calibration_group")

plot_response_surface(data, value = "response", model = NULL,
  sequence = NULL, main = NULL)

plot_metric_surface(metrics, metric, model = NULL,
  sequence = NULL, main = NULL)

build_effect_space_projection(data, reciprocal = NULL, value = "response")

plot_effect_space_projection(projection, metric = "delta_history",
  model = NULL, main = NULL)

plot_cross_model_effect_space(projection, metric = "delta_history",
  models = NULL, main = "Cross-model effect-space overlay")
}
\arguments{
  \item{data}{Long-form reciprocal sequence-response data.}
  \item{outdir}{Directory receiving machine-readable tables and PDF figures.}
  \item{value}{Readout or normalized-response column.}
  \item{blank_col}{Optional row-wise blank column used during normalization.}
  \item{lloq}{Explicit scalar lower limit of quantification after blank correction.}
  \item{lloq_col}{Optional row-wise lower-limit-of-quantification column.}
  \item{baseline}{Normalization baseline passed to \code{normalize_sequence_data}.}
  \item{require_absolute}{Whether calibrated cross-sequence absolute comparison is required. If unavailable, the report records a not-generated manifest entry rather than fabricating a value.}
  \item{calibration_col}{Column identifying a shared reciprocal calibration group.}
  \item{model}{Optional biological model selector.}
  \item{sequence}{Optional reciprocal sequence selector.}
  \item{main}{Optional plot title.}
  \item{metrics}{Replicate-level metric table.}
  \item{metric}{Numeric metric to visualize.}
  \item{reciprocal}{Optional output of \code{estimate_reciprocal_asymmetry}.}
  \item{projection}{Output of \code{build_effect_space_projection}.}
  \item{models}{Optional vector of models to include in a cross-model overlay.}
}
\value{
Plotting functions return plotted data invisibly. \code{build_effect_space_projection}
returns replicate-level observed effect-space coordinates merged with reciprocal history
asymmetry. \code{run_reciprocal_report} invisibly returns a list containing collapsed,
normalized, history, reciprocal, absolute, effect-space, manifest, and error-status tables.
}
\details{
The report is deliberately raw-first. Missing or below-LLOQ values remain missing;
measured grid points are drawn explicitly; no smoothing or interpolation is performed.
The cross-model effect-space coordinates are descriptive averages of observed reciprocal
single-agent marginal residual fractions and are intended for visualization rather than
replacement of replicate-level inference. The 0.5 and 0.7 reference lines are visual
landmarks only and are not predefined optimal or golden regions.
}
