# Rendered-page checks for the selected final R1-R3 simulation reports.

pages <- c(
  "revision_genotype_level_simulation.html",
  "revision_combined_spiky_genotype_simulation.html",
  "revision_functional_testing_simulation.html"
)

for (page in pages) {
  path <- file.path("docs", page)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  stopifnot(
    grepl("rmd-show-all-code", html, fixed = TRUE),
    !grepl("Unable to load file", html, fixed = TRUE)
  )
  image_matches <- regmatches(
    html,
    gregexpr('src="figure/[^"]+', html, perl = TRUE)
  )[[1L]]
  image_paths <- sub('^src="', "", image_matches)
  stopifnot(
    length(image_paths) > 0L,
    all(file.exists(file.path("docs", image_paths)))
  )
  message(
    page,
    ": ",
    length(image_paths),
    " local figure references verified."
  )
}

r3_html <- paste(
  readLines(
    file.path("docs", "revision_functional_testing_simulation.html"),
    warn = FALSE
  ),
  collapse = "\n"
)
r3_source <- paste(
  readLines(
    file.path("analysis", "revision_functional_testing_simulation.rmd"),
    warn = FALSE
  ),
  collapse = "\n"
)
stopifnot(
  !grepl("stephenslab/fashr-paper", r3_html, fixed = TRUE),
  !grepl(
    "center_aligned_relative_clearance_main_effect_fashr0143_pilot5",
    r3_html,
    fixed = TRUE
  ),
  grepl("through nominal alpha 0.20", r3_html, fixed = TRUE),
  grepl(
    "directly over all[[:space:]]*<code>6,362</code>",
    r3_html,
    perl = TRUE
  ),
  grepl("22_r3_full_universe_functional_fashr0143", r3_html, fixed = TRUE),
  !grepl("Apply global dynamic-eQTL FDR control", r3_html, fixed = TRUE),
  !grepl("For each dynamically selected variant", r3_html, fixed = TRUE),
  grepl(
    paste0(
      "BF ribbons are pointwise minima and maxima",
      "[[:space:]]+across five[[:space:]]+replications"
    ),
    r3_html,
    perl = TRUE
  ),
  grepl(
    "interval_methods = &quot;FASH-IWP1-BF&quot;",
    r3_html,
    fixed = TRUE
  ),
  grepl("power_replication_min", r3_html, fixed = TRUE),
  grepl("empirical_fsr_replication_min", r3_html, fixed = TRUE),
  !grepl("Shading shows 95% intervals across five seeds", r3_html, fixed = TRUE),
  !grepl("five-seed Monte Carlo means and 95% t-based", r3_html, fixed = TRUE),
  !grepl("reviewer-facing curves display", r3_html, fixed = TRUE),
  !grepl("finite-category", r3_html, fixed = TRUE),
  !grepl("prespecified", r3_html, fixed = TRUE),
  !grepl("scientific-validation", r3_html, fixed = TRUE),
  !grepl('<div id="interpretation"', r3_html, fixed = TRUE)
)

stopifnot(
  lengths(regmatches(
    r3_source,
    gregexpr('candidate_scope = "full_universe"', r3_source, fixed = TRUE)
  )) == 2L,
  grepl(
    'r3a_alpha$method == "FASH-IWP1-BF"',
    r3_source,
    fixed = TRUE
  ),
  grepl(
    'r3b_alpha$method == "FASH-IWP1-BF"',
    r3_source,
    fixed = TRUE
  )
)

for (figure_number in 1:6) {
  pattern <- paste0('alt="Figure ', figure_number, "[.]")
  match_positions <- gregexpr(pattern, r3_html, perl = TRUE)[[1L]]
  stopifnot(length(match_positions) == 1L, match_positions[[1L]] > 0L)
}
for (table_number in 1:5) {
  stopifnot(grepl(
    paste0("<caption>Table ", table_number, "[.]"),
    r3_html,
    perl = TRUE
  ))
}

message("Selected R1-R3 rendered-page checks passed.")
