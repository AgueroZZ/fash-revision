#!/usr/bin/env Rscript

# Tests for the reporting layer of the temporal-class Roadmap enrichment page.
#
# Sourcing reporting.R runs the cache-validation gate, so reaching the first
# check is itself an assertion. What follows pins the accessors the page's prose
# calls, the class-label mapping that silently produced NAs once already, and
# that every figure and table builds and renders.
#
#   Rscript --vanilla test_reporting.R

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "temporal_class_roadmap_enrichment", "reporting.R"
))

failures <- 0L
check <- function(description, condition) {
  passed <- isTRUE(condition)
  if (!passed) failures <<- failures + 1L
  message(if (passed) "PASS  " else "FAIL  ", description)
  invisible(passed)
}

check("the cache validation gate passed", TRUE)

# ---- Scope ------------------------------------------------------------------

check(
  "the design is 676 estimates over 13 epigenomes and 2 definitions",
  nrow(enrichment) == 676L &&
    length(unique(enrichment$epigenome_id)) == 13L &&
    length(unique(enrichment$enhancer_definition)) == 2L
)

check(
  "all four classes plus the non-switch contrast are present",
  setequal(as.character(unique(enrichment$class)),
           c("early", "middle", "late", "switch", "no_switch"))
)

check(
  "both comparators and both lead rules are present",
  setequal(unique(enrichment$comparator),
           c("Tested background", "Other FASH discoveries")) &&
    setequal(unique(enrichment$lead_rule), c("lfsr", "lfdr"))
)

# This is the bug that shipped once: unname() stripped the names the display
# lookup indexes by, so every class_label came out NA and every figure went grey.
check(
  "every row has a non-NA class display label",
  !anyNA(enrichment$class_label) &&
    setequal(as.character(unique(enrichment$class_label)),
             c("Early", "Middle", "Late", "Switch", "Non-switch"))
)

check(
  "class labels map one-to-one onto classes",
  {
    mapping <- unique(enrichment[, c("class", "class_label")])
    nrow(mapping) == 5L && !anyDuplicated(mapping$class_label)
  }
)

# ---- Consistency with the sister pages -------------------------------------

check(
  "pooling the three timing classes reproduces the 13-epigenome current_all",
  {
    gap <- abs(pooled_crosscheck$fold_enrichment_pooled -
                 pooled_crosscheck$fold_enrichment_published)
    nrow(pooled_crosscheck) == 26L && max(gap, na.rm = TRUE) < 1e-10
  }
)

check(
  "the timing class counts match the three-epigenome sister page",
  {
    counts <- stats::setNames(class_summary$pairs, class_summary$class)
    counts[["early"]] == 1157L && counts[["middle"]] == 1943L &&
      counts[["late"]] == 6114L
  }
)

check(
  "the discovery set is R6's",
  discovered_pairs == 9214L && discovered_variants == 9148L
)

# ---- Switch is orthogonal, not a fourth bucket -----------------------------

check(
  "the switch class is flagged as not mutually exclusive",
  !class_summary$mutually_exclusive[class_summary$class == "switch"] &&
    all(class_summary$mutually_exclusive[class_summary$class != "switch"])
)

check(
  "switch overlaps every timing window at a similar rate",
  {
    inner <- switch_crosstab[switch_crosstab$timing_class != "all", ,
                             drop = FALSE]
    nrow(inner) == 3L && all(inner$switch_pairs > 100L) &&
      max(inner$switch_share) - min(inner$switch_share) < 0.05
  }
)

check(
  "the timing pair counts sum to the discovery set",
  sum(switch_crosstab$pairs[switch_crosstab$timing_class != "all"]) ==
    discovered_pairs
)

# ---- Accessors agree with the cache -----------------------------------------

check(
  "describe_fold reproduces the cached estimate",
  {
    row <- enrichment[
      as.character(enrichment$class) == "early" &
        enrichment$strategy == "all" & enrichment$epigenome_id == "E020" &
        enrichment$comparator == "Tested background" &
        as.character(enrichment$enhancer_definition) == "R6" &
        enrichment$lead_rule == "lfsr",
    ]
    nrow(row) == 1L &&
      identical(
        describe_fold("early", "E020"),
        paste0(
          formatC(row$fold_enrichment, format = "f", digits = 2),
          "-fold (95% CI ",
          formatC(row$ci_lower_fold, format = "f", digits = 2), "-",
          formatC(row$ci_upper_fold, format = "f", digits = 2), ")"
        )
      )
  }
)

check(
  "count_below_one counts intervals strictly below one",
  {
    rows <- select_rows("early", strategy = "all")
    rows <- rows[rows$estimable, ]
    identical(
      count_below_one("early", strategy = "all"),
      paste0(sum(rows$ci_upper_fold < 1), " of ", nrow(rows))
    )
  }
)

check(
  "count_above_one counts intervals strictly above one",
  {
    rows <- select_rows("late", strategy = "all", group = "Heart")
    rows <- rows[rows$estimable, ]
    identical(
      count_above_one("late", strategy = "all", group = "Heart"),
      paste0(sum(rows$ci_lower_fold > 1), " of ", nrow(rows))
    )
  }
)

check(
  "count_estimable reports the power limit rather than hiding it",
  {
    rows <- select_rows("early", strategy = "lead")
    identical(
      count_estimable("early", strategy = "lead"),
      paste0(sum(rows$estimable), " of ", nrow(rows))
    ) && sum(rows$estimable) < nrow(rows)
  }
)

check(
  "a non-estimable cell reports itself as such, not as a number",
  identical(describe_fold("early", "E065", strategy = "all"), "not estimable")
)

check(
  "fold_range brackets the estimable members of a group",
  {
    rows <- select_rows("late", strategy = "all", group = "Heart")
    rows <- rows[rows$estimable, ]
    identical(
      fold_range("late", strategy = "all", group = "Heart"),
      paste0(
        formatC(min(rows$fold_enrichment), format = "f", digits = 2), "-",
        formatC(max(rows$fold_enrichment), format = "f", digits = 2), "-fold"
      )
    )
  }
)

check(
  "fold_list names every member of a group",
  {
    text <- fold_list("late", strategy = "all", group = "Heart")
    all(vapply(c("E083", "E095", "E104", "E105"), grepl, logical(1),
               x = text, fixed = TRUE))
  }
)

check(
  "an empty query is refused rather than returning nothing",
  inherits(try(select_rows("early", epigenome = "E999"), silent = TRUE),
           "try-error")
)

check(
  "set sizes agree with the pre-specified roster",
  identical(class_size("early", "all"), "1,157") &&
    identical(class_size("switch", "all"), "981") &&
    identical(class_size("late", "lead"), "932") &&
    identical(class_size("switch", "lead"), "250")
)

# ---- The two sensitivity axes ----------------------------------------------

definition <- definition_agreement()
lead_rule <- lead_rule_sensitivity()

check(
  "the state definitions agree closely",
  definition$n_cells > 250L && definition$correlation > 0.95
)

check(
  "the lead rules agree much less well than the state definitions do",
  lead_rule$correlation < definition$correlation &&
    lead_rule$max_absolute_shift > 0.2 &&
    lead_rule$min_jaccard < 0.6
)

check(
  "the lead-rule comparison only uses cells estimable under both rules",
  lead_rule$n_cells <= lead_rule$n_total && lead_rule$n_cells > 0L
)

# ---- Figures and tables ----------------------------------------------------

figure_builders <- list(
  function() plot_classes_vs_background("R6", "lfsr"),
  function() plot_classes_vs_background("Strober", "lfsr"),
  function() plot_classes_vs_discoveries("R6"),
  function() plot_switch_contrast("R6"),
  function() plot_lead_rule_sensitivity("R6"),
  function() plot_definition_sensitivity("all"),
  function() plot_fold_heatmap("Tested background", "R6", "lfsr"),
  function() plot_fold_heatmap("Other FASH discoveries", "R6", "lfsr")
)

check(
  "every figure builder returns a ggplot",
  all(vapply(figure_builders, function(builder) {
    inherits(builder(), "ggplot")
  }, logical(1)))
)

check(
  "every figure renders without error",
  {
    device <- tempfile(fileext = ".png")
    grDevices::png(device, width = 1200, height = 1400, res = 120)
    on.exit(
      {
        if (!is.null(grDevices::dev.list())) grDevices::dev.off()
        unlink(device)
      },
      add = TRUE
    )
    for (builder in figure_builders) print(builder())
    TRUE
  }
)

check(
  "no figure silently drops a class it was asked to draw",
  {
    built <- plot_classes_vs_background("R6", "lfsr")
    setequal(
      as.character(unique(built$data$series)),
      c("Early", "Middle", "Late", "Switch")
    )
  }
)

check(
  "every table builds",
  all(vapply(
    list(
      class_table, switch_crosstab_table, estimability_table, lead_rule_table,
      not_estimable_table,
      function() full_summary_table("Tested background"),
      function() full_summary_table("Other FASH discoveries"),
      function() full_summary_table("Tested background", "Strober")
    ),
    function(builder) !is.null(builder()),
    logical(1)
  ))
)

check(
  "the full table has one row per class x strategy x epigenome",
  length(capture.output(full_summary_table("Tested background"))) >=
    4L * 2L * 13L
)

message("")
if (failures > 0L) {
  stop(failures, " test(s) failed.")
}
message("All tests passed.")
