# Regression checks for the fashr 0.1.43 R1/R2 reporting overlay.

source("code/revision_simulations/shared/simulation_functions.R")
source("code/revision_simulations/shared/real_genotype_one_per_gene.R")

metric_for <- function(table, method, column) {
  row <- table[table$method == method, , drop = FALSE]
  if (nrow(row) != 1L || !column %in% names(row)) {
    stop("Could not extract a unique metric for ", method, ".", call. = FALSE)
  }
  row[[column]]
}

require_close <- function(observed, expected, label, tolerance = 1e-14) {
  if (!isTRUE(all.equal(observed, expected, tolerance = tolerance))) {
    stop(label, " did not match its retained value.", call. = FALSE)
  }
}

require_data_equal <- function(observed, expected, label) {
  comparison <- all.equal(
    observed,
    expected,
    tolerance = 0,
    check.attributes = FALSE
  )
  if (!isTRUE(comparison)) {
    stop(label, " changed: ", paste(comparison, collapse = "; "), call. = FALSE)
  }
}

source("code/revision_simulations/r1_random_bspline/reporting.R")

stopifnot(
  identical(current_fash_manifest$result_id, "r1_r2_fashr0143"),
  identical(
    current_fash_manifest$package_provenance$version,
    R1_R2_FASHR0143_VERSION
  ),
  identical(
    current_fash_manifest$package_provenance$remote_sha,
    R1_R2_FASHR0143_REMOTE_SHA
  ),
  all(result_provenance_table$Status == c(
    "Recomputed on Midway3",
    "Reused; not rerun",
    "Strict pairing check passed"
  ))
)
require_data_equal(
  direct_alpha[grepl("^Direct-", direct_alpha$method), ],
  historical_direct_alpha[grepl("^Direct-", historical_direct_alpha$method), ],
  "R1 historical direct curve"
)
require_data_equal(
  direct_alpha_005[grepl("^Direct-", direct_alpha_005$method), ],
  historical_direct_alpha_005[
    grepl("^Direct-", historical_direct_alpha_005$method),
  ],
  "R1 historical direct alpha-0.05 summary"
)
require_close(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_discoveries"),
  1230,
  "R1 IWP1-BF mean discoveries"
)
require_close(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_power"),
  0.93066037735849105,
  "R1 IWP1-BF power"
)
require_close(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_fdr"),
  0.037563172411889202,
  "R1 IWP1-BF empirical FDR"
)
require_close(
  metric_for(fash_alpha_005, "FASH-linear-BF", "mean_power"),
  0.47201257861635199,
  "R1 linear-BF power"
)
stopifnot(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_power") >
    metric_for(
      historical_fash_alpha[historical_fash_alpha$alpha == 0.05, ],
      "FASH-IWP1-BF",
      "mean_power"
    )
)

source("code/revision_simulations/r2_spiky_transient/reporting.R")

stopifnot(
  identical(current_fash_manifest$result_id, "r1_r2_fashr0143")
)
require_data_equal(
  direct_alpha[grepl("^Direct-", direct_alpha$method), ],
  historical_direct_alpha[grepl("^Direct-", historical_direct_alpha$method), ],
  "R2 historical direct curve"
)
require_data_equal(
  direct_alpha_005[grepl("^Direct-", direct_alpha_005$method), ],
  historical_direct_alpha_005[
    grepl("^Direct-", historical_direct_alpha_005$method),
  ],
  "R2 historical direct alpha-0.05 summary"
)
require_data_equal(
  peak_alpha[grepl("^Direct-", peak_alpha$method), ],
  historical_peak_alpha[grepl("^Direct-", historical_peak_alpha$method), ],
  "R2 historical direct peak curve"
)
require_data_equal(
  peak_alpha_005[grepl("^Direct-", peak_alpha_005$method), ],
  historical_peak_alpha_005[
    grepl("^Direct-", historical_peak_alpha_005$method),
  ],
  "R2 historical direct peak alpha-0.05 summary"
)
require_close(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_discoveries"),
  659.8,
  "R2 IWP1-BF mean discoveries"
)
require_close(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_power"),
  0.51666666666666705,
  "R2 IWP1-BF power"
)
require_close(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_fdr"),
  0.0039456054476991404,
  "R2 IWP1-BF empirical FDR"
)
require_close(
  metric_for(fash_alpha_005, "FASH-linear-BF", "mean_power"),
  0.019025157232704398,
  "R2 linear-BF power"
)
stopifnot(
  metric_for(fash_alpha_005, "FASH-IWP1-BF", "mean_power") >
    metric_for(
      historical_fash_alpha[historical_fash_alpha$alpha == 0.05, ],
      "FASH-IWP1-BF",
      "mean_power"
    )
)

message("R1/R2 fashr 0.1.43 reporting overlay tests passed.")
