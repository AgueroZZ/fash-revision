#!/usr/bin/env Rscript

# Summarize prespecified R3 Middle repair candidates using retained artifacts.
# This is an evidence table, not a formal simulation result.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), mustWork = TRUE)
diagnostic_dir <- file.path(
  project_root,
  "output/revision_simulations/diagnostics/r3_middle_calibration"
)

mixture <- read.csv(file.path(
  diagnostic_dir,
  "middle_mixture_sensitivity_summary.csv"
))
false_calls <- read.csv(file.path(
  diagnostic_dir,
  "false_call_decomposition.csv"
))
alpha005 <- read.csv(file.path(
  diagnostic_dir,
  "alpha005_all_target_mixture_sensitivity.csv"
))
odds_correction <- read.csv(file.path(
  diagnostic_dir,
  "alpha005_middle_inference_odds_correction.csv"
))
reference_mixture <- read.csv(file.path(
  diagnostic_dir,
  "reference_category_mixture.csv"
))

bf_mixture <- mixture[
  mixture$method == "FASH-IWP1-BF" & mixture$alpha >= 0.05,
  ,
  drop = FALSE
]
by_mechanism <- split(bf_mixture, bf_mixture$truth_mechanism)

metric_row <- function(mechanism) {
  rows <- by_mechanism[[mechanism]]
  data.frame(
    truth_mechanism = mechanism,
    current_max_excess = max(rows$observed_mean - rows$alpha),
    current_alpha_at_max =
      rows$alpha[which.max(rows$observed_mean - rows$alpha)],
    matched_mixture_max_excess =
      max(rows$matched_mixture_mean - rows$alpha),
    matched_mixture_alpha_at_max = rows$alpha[
      which.max(rows$matched_mixture_mean - rows$alpha)
    ],
    stringsAsFactors = FALSE
  )
}
mechanism_metrics <- do.call(rbind, lapply(
  names(by_mechanism),
  metric_row
))

far_false_calls <- sum(false_calls$false_calls[
  false_calls$truth_distance_bin %in%
    c("less_than_-0.50", "-0.50_to_-0.25")
])
total_false_calls <- unique(false_calls$total_false_calls)
if (length(total_false_calls) != 1L) {
  stop("False-call diagnostics contain inconsistent totals.")
}

middle_alpha005 <- alpha005[
  alpha005$truth_mechanism == "random_bspline" &
    alpha005$method == "FASH-IWP1-BF" &
    alpha005$target == "middle",
  ,
  drop = FALSE
]
middle_odds_alpha005 <- odds_correction[
  odds_correction$truth_mechanism == "random_bspline" &
    odds_correction$method == "FASH-IWP1-BF",
  ,
  drop = FALSE
]

candidate_table <- data.frame(
  candidate = c(
    "Current two-stage one-sided posterior",
    "Remove first-stage null calls only",
    "Uncorrected exhaustive three-category posterior",
    "Increase truth margin",
    "Equal-category prior-odds correction in inference",
    "Prior-geometry-matched truth-category mixture"
  ),
  estimand_changed = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
  formal_inference_changed = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
  evidence = c(
    paste0(
      "R3A BF maximum excess = ",
      sprintf(
        "%.3f",
        mechanism_metrics$current_max_excess[
          mechanism_metrics$truth_mechanism == "random_bspline"
        ]
      ),
      "."
    ),
    paste0(
      "Conditional-only R3A BF FSR remains approximately 0.263 at alpha 0.20; ",
      "first-stage null calls explain only a small part of the excess."
    ),
    paste0(
      "With exhaustive Early/Middle/Late support and negligible ties, ",
      "P(Middle functional > 0) equals the Middle category probability."
    ),
    paste0(
      far_false_calls,
      " of ",
      total_false_calls,
      " alpha-0.05 R3A BF Middle false calls are at least 0.25 below zero."
    ),
    paste0(
      "At alpha 0.05, re-thresholding to equal-category posterior odds ",
      "reduced mean R3A BF empirical FSR from ",
      sprintf("%.3f", mean(middle_alpha005$observed_empirical_fsr)),
      " to ",
      sprintf(
        "%.3f",
        mean(middle_odds_alpha005$odds_corrected_empirical_fsr)
      ),
      " and changes the real-data posterior interpretation."
    ),
    paste0(
      "Reference probabilities are Early ",
      sprintf(
        "%.4f",
        reference_mixture$canonical_reference_probability[
          reference_mixture$category == "early"
        ]
      ),
      ", Middle ",
      sprintf(
        "%.4f",
        reference_mixture$canonical_reference_probability[
          reference_mixture$category == "middle"
        ]
      ),
      ", Late ",
      sprintf(
        "%.4f",
        reference_mixture$canonical_reference_probability[
          reference_mixture$category == "late"
        ]
      ),
      ". ",
      "Importance-weighted maximum Middle excess is ",
      sprintf(
        "%.3f",
        mechanism_metrics$matched_mixture_max_excess[
          mechanism_metrics$truth_mechanism == "random_bspline"
        ]
      ),
      " for R3A and ",
      sprintf(
        "%.3f",
        mechanism_metrics$matched_mixture_max_excess[
          mechanism_metrics$truth_mechanism == "raised_cosine"
        ]
      ),
      " for the retained center-aligned R3B sensitivity."
    )
  ),
  status = c(
    "fails calibration gate",
    "insufficient",
    "equivalent; insufficient",
    "diagnostic only; insufficient",
    "eligible but changes interpretation and is conservative at alpha 0.05",
    "recommended for a new primary calibration scenario"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  candidate_table,
  file.path(diagnostic_dir, "candidate_comparison.csv"),
  row.names = FALSE
)
write.csv(
  mechanism_metrics,
  file.path(diagnostic_dir, "candidate_calibration_gate_metrics.csv"),
  row.names = FALSE
)

cat("Wrote candidate repair evidence and calibration-gate metrics to:\n")
cat(diagnostic_dir, "\n")
