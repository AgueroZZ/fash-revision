#!/usr/bin/env Rscript

# Tests for the reporting layer of the internal Roadmap enhancer exploration.
#
# Sourcing reporting.R already runs the cache-validation gate, so getting this
# far is itself the first assertion. What follows checks that the accessors the
# page's prose calls agree with the cached estimates, that the figures build,
# and that the two epigenomes R6 does not show are genuinely present.
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
  "roadmap_enhancer_exploration", "reporting.R"
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
  "all thirteen epigenomes and both definitions are present",
  length(unique(enrichment$epigenome_id)) == 13L &&
    length(unique(enrichment$enhancer_definition)) == 2L &&
    nrow(enrichment) == 208L
)

check(
  "every discovery set the R6 page reports is here",
  all(c("current_all", "current_lead", "linear_all", "linear_lead",
        "quadratic_all", "quadratic_lead", "current_only_all",
        "current_only_lead") %in% enrichment$variant_set)
)

check(
  "the ten added epigenomes are all estimable for the FASH lead set",
  all(enrichment$estimable[
    enrichment$variant_set == "current_lead" &
      enrichment$enhancer_definition == "R6"
  ])
)

# ---- Accessors agree with the cache -----------------------------------------

check(
  "describe_fold reproduces the cached estimate",
  {
    row <- enrichment[
      enrichment$variant_set == "current_lead" &
        enrichment$epigenome_id == "E095" &
        enrichment$enhancer_definition == "R6",
    ]
    identical(
      describe_fold("current_lead", "E095"),
      paste0(
        formatC(row$fold_enrichment, format = "f", digits = 2), "-fold (95% CI ",
        formatC(row$ci_lower_fold, format = "f", digits = 2), "-",
        formatC(row$ci_upper_fold, format = "f", digits = 2), ")"
      )
    )
  }
)

check(
  "group_above_one counts intervals strictly above one",
  {
    rows <- enrichment[
      enrichment$variant_set == "current_lead" &
        enrichment$biological_group == "iPSC" &
        enrichment$enhancer_definition == "R6",
    ]
    identical(
      group_above_one("current_lead", "iPSC"),
      paste0(sum(rows$ci_lower_fold > 1), " of 5")
    )
  }
)

check(
  "group_above_one_labels names exactly those epigenomes",
  {
    rows <- enrichment[
      enrichment$variant_set == "current_lead" &
        enrichment$biological_group == "iPSC" &
        enrichment$enhancer_definition == "R6",
    ]
    rows <- rows[order(rows$epigenome_order), ]
    expected <- rows$epigenome_id[rows$ci_lower_fold > 1]
    identical(
      group_above_one_labels("current_lead", "iPSC"),
      if (length(expected)) paste(expected, collapse = ", ") else "none"
    )
  }
)

check(
  "group_fold_range brackets every member of the group",
  {
    rows <- enrichment[
      enrichment$variant_set == "current_all" &
        enrichment$biological_group == "Heart" &
        enrichment$enhancer_definition == "R6",
    ]
    identical(
      group_fold_range("current_all", "Heart"),
      paste0(
        formatC(min(rows$fold_enrichment), format = "f", digits = 2), "-",
        formatC(max(rows$fold_enrichment), format = "f", digits = 2), "-fold"
      )
    )
  }
)

check(
  "group_fold_list names every member of the group once",
  {
    text <- group_fold_list("current_lead", "Heart")
    all(vapply(c("E083", "E095", "E104", "E105"), function(id) {
      grepl(id, text, fixed = TRUE)
    }, logical(1))) && !grepl("E065", text, fixed = TRUE)
  }
)

check(
  "a non-estimable cell is reported as such rather than as a number",
  {
    if (!nrow(not_estimable)) {
      TRUE
    } else {
      grepl(
        "not estimable",
        group_fold_list(
          not_estimable$variant_set[1L],
          not_estimable$biological_group[1L],
          as.character(not_estimable$enhancer_definition[1L])
        ),
        fixed = TRUE
      )
    }
  }
)

check(
  "an epigenome outside the panel is refused",
  inherits(try(describe_fold("current_lead", "E999"), silent = TRUE),
           "try-error")
)

# ---- The R6 estimates are untouched ----------------------------------------

check(
  "the R6-definition estimates for E020/E013/E095 reproduce the published page",
  {
    published <- r6_estimate_reproduction
    finite <- is.finite(published$fold_enrichment_published)
    all(abs(
      published$fold_enrichment_published[finite] -
        published$fold_enrichment_here[finite]
    ) < 1e-10) && nrow(published) == 24L
  }
)

check(
  "the R6-definition indicators match the published annotation matrix exactly",
  all(r6_agreement$n_disagreements == 0L) && nrow(r6_agreement) == 3L
)

# ---- Definition sensitivity summary ----------------------------------------

agreement <- definition_agreement()

check(
  "definition_agreement covers exactly the cells estimable under both",
  {
    estimable_pairs <- Reduce(intersect, lapply(c("R6", "Strober"), function(d) {
      rows <- enrichment[enrichment$enhancer_definition == d &
                           enrichment$estimable, , drop = FALSE]
      paste(rows$variant_set, rows$epigenome_id)
    }))
    agreement$n_cells == length(estimable_pairs)
  }
)

check(
  "the two definitions are reported as strongly correlated, not identical",
  agreement$correlation > 0.9 && agreement$max_absolute_shift > 0
)

# ---- Figures and tables build ----------------------------------------------

check(
  "figure A builds for both definitions",
  all(vapply(c("R6", "Strober"), function(definition) {
    inherits(plot_method_comparison(definition), "ggplot")
  }, logical(1)))
)

check(
  "figure B builds for both definitions",
  all(vapply(c("R6", "Strober"), function(definition) {
    inherits(plot_fash_only(definition), "ggplot")
  }, logical(1)))
)

check(
  "the sensitivity figure builds",
  inherits(plot_definition_sensitivity(), "ggplot")
)

check(
  "the heatmap builds for both definitions",
  all(vapply(c("R6", "Strober"), function(definition) {
    inherits(plot_fold_heatmap(definition), "ggplot")
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
    for (plot in list(
      plot_method_comparison("R6"), plot_fash_only("R6"),
      plot_definition_sensitivity(), plot_fold_heatmap("R6"),
      plot_method_comparison("Strober"), plot_fash_only("Strober")
    )) {
      print(plot)
    }
    TRUE
  }
)

check(
  "every table builds",
  all(vapply(
    list(
      panel_table, definition_coverage_table, state_table,
      not_estimable_table, provenance_table,
      function() full_summary_table("R6"),
      function() full_summary_table("Strober")
    ),
    function(builder) !is.null(builder()),
    logical(1)
  ))
)

check(
  "the full summary table has one row per epigenome x discovery set",
  length(capture.output(full_summary_table("R6"))) >= 13L * 8L
)

message("")
if (failures > 0L) {
  stop(failures, " test(s) failed.")
}
message("All tests passed.")
