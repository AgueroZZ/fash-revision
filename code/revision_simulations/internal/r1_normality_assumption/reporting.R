# Load and validate the two caches used by the internal R1 diagnostic page.

find_workflowr_root <- function() {
  if (file.exists("analysis/index.Rmd")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/analysis/index.Rmd")) {
    return(normalizePath(
      "coderepo-local",
      winslash = "/",
      mustWork = TRUE
    ))
  }
  stop("Could not find the workflowr repository root.")
}

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required R1 diagnostic cache file is missing: ", path)
  }
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_required_rds <- function(path) {
  if (!file.exists(path)) {
    stop("Required R1 diagnostic cache file is missing: ", path)
  }
  readRDS(path)
}

require_reporting_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required for reporting.")
  }
}

format_integer <- function(x) {
  format(
    as.integer(round(x)),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

format_decimal <- function(x, digits = 3L) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

format_percent <- function(x, digits = 1L) {
  paste0(format_decimal(100 * x, digits), "%")
}

render_scrollable_table <- function(data,
                                    caption = NULL,
                                    align = NULL,
                                    minimum_width = "820px") {
  table_html <- knitr::kable(
    data,
    format = "html",
    escape = TRUE,
    caption = caption,
    align = align,
    row.names = FALSE
  )
  table_html <- kableExtra::kable_styling(
    table_html,
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  )
  knitr::asis_output(paste0(
    '<div class="normality-table-scroll">',
    '<div style="min-width:', minimum_width, ';">',
    as.character(table_html),
    "</div></div>"
  ))
}

normality_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey94"),
      legend.position = "bottom",
      plot.title.position = "plot"
    )
}

require_reporting_package("ggplot2")
require_reporting_package("knitr")
require_reporting_package("kableExtra")
require_reporting_package("patchwork")
require_reporting_package("scales")

workflowr_root <- find_workflowr_root()
simulation_cache_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_normality_assumption_seed12345_t5_v1"
)
residual_cache_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_real_data_standardized_residual_normality_days1_3_9_v1"
)

# The simulation cache supplies only the accepted Part 2 consequence analysis.
r1_configuration <- read_required_rds(file.path(
  simulation_cache_directory,
  "configuration.rds"
))
performance_alpha005 <- read_required_csv(file.path(
  simulation_cache_directory,
  "performance_alpha005.csv"
))
performance_alpha_curve <- read_required_csv(file.path(
  simulation_cache_directory,
  "performance_alpha_curve.csv"
))
simulation_validation <- read_required_csv(file.path(
  simulation_cache_directory,
  "validation.csv"
))

if (!is.list(r1_configuration) ||
    !identical(as.integer(r1_configuration$seed), 12345L) ||
    !identical(as.integer(r1_configuration$n_units), 6362L) ||
    !identical(as.integer(r1_configuration$n_donors), 19L) ||
    !identical(as.integer(r1_configuration$n_time), 16L) ||
    !identical(as.integer(r1_configuration$n_covariates), 5L) ||
    !identical(as.integer(r1_configuration$residual_df), 12L) ||
    nrow(simulation_validation) != 11L ||
    !all(simulation_validation$pass)) {
  stop("The retained R1 simulation cache is invalid.")
}

expected_error_distributions <- c("Gaussian", "Standardized t5")
expected_se_scales <- c("raw regression SE", "t-adjusted SE")
expected_fit_stages <- c("Raw", "BF")
if (nrow(performance_alpha005) != 8L ||
    !setequal(
      performance_alpha005$error_distribution,
      expected_error_distributions
    ) ||
    !setequal(performance_alpha005$se_scale, expected_se_scales) ||
    !setequal(performance_alpha005$fit_stage, expected_fit_stages) ||
    any(abs(performance_alpha005$alpha - 0.05) > 1e-12) ||
    any(performance_alpha005$n_units != 6362L) ||
    any(performance_alpha005$realized_fdp < 0 |
          performance_alpha005$realized_fdp > 1)) {
  stop("The retained alpha-0.05 performance table is invalid.")
}
performance_group <- interaction(
  performance_alpha_curve$error_distribution,
  performance_alpha_curve$se_scale,
  performance_alpha_curve$fit_stage,
  drop = TRUE
)
if (length(unique(table(performance_group))) != 1L ||
    !setequal(
      performance_alpha_curve$error_distribution,
      expected_error_distributions
    ) ||
    !setequal(performance_alpha_curve$se_scale, expected_se_scales) ||
    !setequal(performance_alpha_curve$fit_stage, expected_fit_stages) ||
    any(!is.finite(performance_alpha_curve$alpha)) ||
    any(!is.finite(performance_alpha_curve$realized_fdp))) {
  stop("The retained performance curves are invalid.")
}

gaussian_performance_alpha005 <- performance_alpha005[
  performance_alpha005$error_distribution == "Gaussian",
  ,
  drop = FALSE
]
gaussian_performance_alpha_curve <- performance_alpha_curve[
  performance_alpha_curve$error_distribution == "Gaussian",
  ,
  drop = FALSE
]

# The real-data cache supplies the traditional Day 1/3/9 residual diagnostic.
residual_configuration <- read_required_rds(file.path(
  residual_cache_directory,
  "configuration.rds"
))
selected_units <- read_required_csv(file.path(
  residual_cache_directory,
  "selected_units.csv"
))
standardized_residuals <- read_required_rds(file.path(
  residual_cache_directory,
  "standardized_residuals.rds"
))
fit_status <- read_required_csv(file.path(
  residual_cache_directory,
  "fit_status.csv"
))
residual_qq <- read_required_csv(file.path(
  residual_cache_directory,
  "residual_qq.csv"
))
residual_histogram <- read_required_csv(file.path(
  residual_cache_directory,
  "residual_histogram.csv"
))
residual_day_summary <- read_required_csv(file.path(
  residual_cache_directory,
  "residual_day_summary.csv"
))
residual_tail_summary <- read_required_csv(file.path(
  residual_cache_directory,
  "residual_tail_summary.csv"
))
residual_validation <- read_required_csv(file.path(
  residual_cache_directory,
  "validation.csv"
))

expected_days <- c(1L, 3L, 9L)
expected_day_labels <- paste("Day", expected_days)
if (!is.list(residual_configuration) ||
    !identical(as.integer(residual_configuration$analysis_days), expected_days) ||
    !identical(as.integer(residual_configuration$n_available_units), 6362L) ||
    !identical(as.integer(residual_configuration$n_selected_units), 6362L) ||
    !identical(as.integer(residual_configuration$n_donors_per_day), 19L) ||
    !identical(as.integer(residual_configuration$n_parameters), 7L) ||
    !identical(as.integer(residual_configuration$residual_df), 12L) ||
    !identical(
      residual_configuration$residual_definition,
      "e_i / {s * sqrt(1 - h_ii)}"
    ) ||
    nrow(residual_validation) != 9L ||
    !all(residual_validation$pass)) {
  stop("The real-data residual configuration is invalid.")
}

if (nrow(selected_units) != 6362L ||
    anyDuplicated(selected_units$unit_key) ||
    anyDuplicated(selected_units$gene_id) ||
    !setequal(unique(fit_status$day), expected_days) ||
    nrow(fit_status) != 3L * 6362L ||
    !all(fit_status$valid_fit) ||
    nrow(standardized_residuals) != 3L * 6362L * 19L ||
    !setequal(unique(standardized_residuals$day), expected_days) ||
    !setequal(unique(standardized_residuals$day_label), expected_day_labels) ||
    any(standardized_residuals$residual_df != 12L) ||
    any(!is.finite(standardized_residuals$standardized_residual))) {
  stop("The cached real-data standardized residuals are invalid.")
}
residual_fit_group <- interaction(
  standardized_residuals$day,
  standardized_residuals$unit_key,
  drop = TRUE
)
if (length(residual_fit_group) != nrow(standardized_residuals) ||
    any(table(residual_fit_group) != 19L)) {
  stop("The residual denominators do not match 19 donors per regression.")
}

qq_group <- split(residual_qq$probability, residual_qq$day)
if (!setequal(unique(residual_qq$day), expected_days) ||
    length(unique(vapply(qq_group, length, integer(1L)))) != 1L ||
    max(abs(
      residual_qq$normal_reference_quantile -
        stats::qnorm(residual_qq$probability)
    )) > 1e-12 ||
    any(!is.finite(residual_qq$empirical_quantile))) {
  stop("The cached residual QQ quantiles are invalid.")
}
if (!setequal(unique(residual_histogram$day), expected_days) ||
    any(!is.finite(residual_histogram$density)) ||
    any(residual_histogram$count < 0) ||
    !all(vapply(expected_days, function(day) {
      sum(residual_histogram$count[residual_histogram$day == day]) ==
        sum(standardized_residuals$day == day)
    }, logical(1L)))) {
  stop("The cached residual histograms are invalid.")
}
if (nrow(residual_day_summary) != 3L ||
    !setequal(residual_day_summary$day, expected_days) ||
    any(residual_day_summary$n_valid_fits != 6362L) ||
    any(residual_day_summary$n_residuals != 6362L * 19L) ||
    nrow(residual_tail_summary) != 9L ||
    !setequal(residual_tail_summary$day, expected_days) ||
    !setequal(
      residual_tail_summary$nominal_two_sided_rate,
      c(0.05, 0.01, 0.001)
    ) ||
    any(residual_tail_summary$observed_two_sided_rate < 0 |
          residual_tail_summary$observed_two_sided_rate > 1)) {
  stop("The cached residual summaries are invalid.")
}
