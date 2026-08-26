#!/usr/bin/env Rscript

# Unit tests for the temporal-class Roadmap enrichment helpers.
#
# The things that could silently change a number: the switch call, the
# orthogonality of switch to the timing partition, the two lead-selection rules,
# and the key-based join to the R6 discovery set. Synthetic fixtures first, then
# a handful of assertions against the real cached inputs.
#
#   Rscript --vanilla test_temporal_class_roadmap_helpers.R

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
  "temporal_class_roadmap_enrichment", "temporal_class_roadmap_helpers.R"
))

failures <- 0L
check <- function(description, condition) {
  passed <- isTRUE(condition)
  if (!passed) failures <<- failures + 1L
  message(if (passed) "PASS  " else "FAIL  ", description)
  invisible(passed)
}

# ---- Constants --------------------------------------------------------------

check(
  "the four classes are early, middle, late, switch",
  identical(TCR_CLASSES, c("early", "middle", "late", "switch"))
)

check(
  "the three timing classes are inherited from the sister page unchanged",
  identical(TEMPORAL_CLASSES, c("early", "middle", "late")) &&
    identical(TCR_CLASSES[1:3], TEMPORAL_CLASSES)
)

check(
  "the switch threshold matches code/02_dyn_lfsr.R",
  identical(SWITCH_THRESHOLD, 0.25)
)

check(
  "both lead rules are declared",
  identical(sort(names(LEAD_RULES)), c("lfdr", "lfsr"))
)

# ---- Synthetic fixture ------------------------------------------------------

# Six pairs over three genes. Window probabilities are built so that the argmax
# is known by construction: pairs 1-2 early, 3-4 middle, 5-6 late.
fixture_keys <- c("G1_rsA", "G1_rsB", "G2_rsC", "G2_rsD", "G3_rsE", "G3_rsF")
make_result <- function(lfsr, cfsr = lfsr) {
  data.frame(indices = seq_along(lfsr), lfsr = lfsr, cfsr = cfsr,
             row.names = fixture_keys)
}
fixture <- list(
  #                     1     2     3     4     5     6
  early  = make_result(c(0.10, 0.20, 0.90, 0.90, 0.90, 0.90)),
  middle = make_result(c(0.80, 0.80, 0.10, 0.30, 0.95, 0.95)),
  late   = make_result(c(0.85, 0.85, 0.85, 0.85, 0.05, 0.15)),
  # switch is called on cfsr; pairs 1, 4 and 5 are switches.
  switch = make_result(
    lfsr = c(0.01, 0.60, 0.70, 0.02, 0.03, 0.80),
    cfsr = c(0.01, 0.60, 0.70, 0.04, 0.02, 0.80)
  )
)
fixture_pairs <- data.frame(
  pair_key = fixture_keys,
  gene_id = c("G1", "G1", "G2", "G2", "G3", "G3"),
  variant_id = c("rsA", "rsB", "rsC", "rsD", "rsE", "rsF"),
  # Overall lfdr deliberately ranks the SECOND variant of each gene best, so the
  # two lead rules must disagree.
  score = c(0.04, 0.01, 0.04, 0.01, 0.04, 0.01),
  stringsAsFactors = FALSE
)

membership <- build_class_membership(fixture_pairs, fixture)

check(
  "the timing argmax is as constructed",
  identical(
    as.character(membership$timing_class[order(membership$variant_id)]),
    c("early", "early", "middle", "middle", "late", "late")
  )
)

check(
  "the timing classes partition the pairs",
  all(rowSums(as.matrix(
    membership[, paste0("in_", TEMPORAL_CLASSES), drop = FALSE]
  )) == 1L)
)

check(
  "switch is called on cfsr, not lfsr",
  {
    ordered <- membership[order(membership$variant_id), , drop = FALSE]
    # rsD has lfsr 0.02 and cfsr 0.04 -> in; rsC has cfsr 0.70 -> out.
    identical(ordered$in_switch, c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))
  }
)

check(
  "switch overlaps all three timing classes rather than partitioning them",
  {
    crosstab <- crosstab_switch_by_timing(membership)
    inner <- crosstab[crosstab$timing_class != "all", , drop = FALSE]
    all(inner$switch_pairs > 0) && sum(inner$switch_pairs) == 3L &&
      crosstab$switch_pairs[crosstab$timing_class == "all"] == 3L
  }
)

check(
  "a class can hold a pair that another class also holds",
  any(membership$in_switch & membership$in_late)
)

# ---- Lead selection ---------------------------------------------------------

lfsr_sets <- build_class_variant_sets(membership, lead_rule = "lfsr")
lfdr_sets <- build_class_variant_sets(membership, lead_rule = "lfdr")

check(
  "every class yields an all set and a lead set",
  identical(
    sort(names(lfsr_sets)),
    sort(paste0(rep(TCR_CLASSES, each = 2L), c("_all", "_lead")))
  )
)

check(
  "the lfsr rule picks the variant most confident for its own class",
  # Gene G1 is early; rsA has early lfsr 0.10 against rsB's 0.20.
  identical(lfsr_sets$early_lead, "rsA") &&
    # Gene G2 is middle; rsC has middle lfsr 0.10 against rsD's 0.30.
    identical(lfsr_sets$middle_lead, "rsC") &&
    # Gene G3 is late; rsE has late lfsr 0.05 against rsF's 0.15.
    identical(lfsr_sets$late_lead, "rsE")
)

check(
  "the lfdr rule picks the lowest overall score instead, and disagrees here",
  identical(lfdr_sets$early_lead, "rsB") &&
    identical(lfdr_sets$middle_lead, "rsD") &&
    identical(lfdr_sets$late_lead, "rsF") &&
    !identical(lfsr_sets$early_lead, lfdr_sets$early_lead)
)

check(
  "the all sets are identical under both rules",
  identical(
    lfsr_sets[paste0(TCR_CLASSES, "_all")],
    lfdr_sets[paste0(TCR_CLASSES, "_all")]
  )
)

check(
  "a lead set never exceeds its own all set",
  all(vapply(TCR_CLASSES, function(class_name) {
    lead <- lfsr_sets[[paste0(class_name, "_lead")]]
    all(lead %in% lfsr_sets[[paste0(class_name, "_all")]])
  }, logical(1)))
)

check(
  "one lead per gene per class",
  {
    inside <- membership[membership$in_switch, , drop = FALSE]
    length(lfsr_sets$switch_lead) <= length(unique(inside$gene_id))
  }
)

check(
  "an unknown lead rule is refused",
  inherits(
    try(build_class_variant_sets(membership, lead_rule = "pvalue"),
        silent = TRUE),
    "try-error"
  )
)

no_switch <- build_no_switch_sets(membership, lead_rule = "lfsr")

check(
  "the no-switch set is the exact complement of the switch set",
  setequal(
    union(no_switch$no_switch_all, lfsr_sets$switch_all),
    membership$variant_id
  ) && !length(intersect(no_switch$no_switch_all, lfsr_sets$switch_all))
)

check(
  "set names round-trip to class and strategy",
  {
    parsed <- parse_set_name(c("switch_lead", "no_switch_all", "middle_all"))
    identical(parsed$class, c("switch", "no_switch", "middle")) &&
      identical(parsed$strategy, c("lead", "all", "all"))
  }
)

# ---- Guards -----------------------------------------------------------------

check(
  "a discovery set that does not match the classification is refused",
  inherits(
    try(
      build_class_membership(fixture_pairs[1:3, , drop = FALSE], fixture),
      silent = TRUE
    ),
    "try-error"
  )
)

check(
  "misaligned classification row names are refused",
  {
    broken <- fixture
    rownames(broken$switch) <- rev(fixture_keys)
    inherits(
      try(build_class_membership(fixture_pairs, broken), silent = TRUE),
      "try-error"
    )
  }
)

check(
  "an out-of-range switch alpha is refused",
  inherits(
    try(build_class_membership(fixture_pairs, fixture, switch_alpha = 1.5),
        silent = TRUE),
    "try-error"
  )
)

# ---- Summary tables ---------------------------------------------------------

summary_table <- summarise_classes(membership, lead_rule = "lfsr")

check(
  "the summary table flags which classes are mutually exclusive",
  identical(
    summary_table$mutually_exclusive,
    c(TRUE, TRUE, TRUE, FALSE)
  )
)

check(
  "the summary table counts agree with the sets",
  all(vapply(seq_len(nrow(summary_table)), function(index) {
    class_name <- summary_table$class[index]
    summary_table$variants[index] ==
      length(lfsr_sets[[paste0(class_name, "_all")]]) &&
      summary_table$leads[index] ==
        length(lfsr_sets[[paste0(class_name, "_lead")]])
  }, logical(1)))
)

# ---- Against the real cached inputs ----------------------------------------

real_classifications <- load_all_classifications()
real_pairs <- readRDS(file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison_fashr0143", "discovery_pair_tables.rds"
))$current_all
real_membership <- build_class_membership(real_pairs, real_classifications)

check(
  "the real classification covers exactly the R6 discovered pairs",
  nrow(real_membership) == 9214L
)

check(
  "the real timing counts match the sister page",
  {
    counts <- table(as.character(real_membership$timing_class))
    counts[["early"]] == 1157L && counts[["middle"]] == 1943L &&
      counts[["late"]] == 6114L
  }
)

check(
  "the real switch call selects 981 pairs at cfsr <= 0.05",
  sum(real_membership$in_switch) == 981L
)

check(
  "switch really is spread across all three timing windows",
  {
    crosstab <- crosstab_switch_by_timing(real_membership)
    inner <- crosstab[crosstab$timing_class != "all", , drop = FALSE]
    all(inner$switch_pairs > 100L) &&
      max(inner$switch_share) - min(inner$switch_share) < 0.05
  }
)

check(
  "the real lead rules disagree enough to be worth reporting",
  {
    by_lfsr <- build_class_variant_sets(real_membership, "lfsr")
    by_lfdr <- build_class_variant_sets(real_membership, "lfdr")
    shared <- length(intersect(by_lfsr$late_lead, by_lfdr$late_lead))
    shared < length(by_lfsr$late_lead) && shared > 0.4 * length(by_lfsr$late_lead)
  }
)

# ---- Which middle window? ---------------------------------------------------

# Historical caches may use either the current open (3, 12) Middle window or
# the superseded closed [4, 11] window. They cannot be identified from the
# saved lfsr filename alone, so the guard tests the partition identity that
# only the open window satisfies.
null_probability <- load_posterior_null_probability()

check(
  "the cached pi0 lines up with the classifications",
  identical(
    unname(null_probability$indices),
    unname(real_classifications$early$indices)
  ) && length(null_probability$posterior_null_probability) == 9214L
)

middle_definition <- verify_middle_definition(
  real_classifications, null_probability
)

check(
  "the cache satisfies P(early)+P(middle)+P(late) = 1 - pi0",
  abs(middle_definition$mean_residual_open) < 1e-3 &&
    abs(middle_definition$median_residual_open) < 1e-3
)

check(
  "so the cache is the open (3, 12) middle window",
  identical(middle_definition$window, "open (3, 12)")
)

check(
  "comparing against 1 instead of 1 - pi0 is the trap, and it is real",
  # This is the mistake that made a correct cache look like a stale one: the
  # residual against 1 is centred at -pi0, not at zero.
  abs(middle_definition$mean_residual_versus_one) >
    50 * abs(middle_definition$mean_residual_open) &&
    abs(middle_definition$mean_residual_versus_one +
          middle_definition$median_pi0) < 0.02
)

check(
  "a closed [4, 11] cache would be rejected",
  {
    # Emulate the stale window by removing the mass the open window assigns to
    # (3, 4) and (11, 12): shrink P(middle) by a plausible gap share.
    stale <- real_classifications
    stale$middle$lfsr <- pmin(1, stale$middle$lfsr + 0.04)
    inherits(
      try(verify_middle_definition(stale, null_probability), silent = TRUE),
      "try-error"
    )
  }
)

message("")
if (failures > 0L) {
  stop(failures, " test(s) failed.")
}
message("All tests passed.")
