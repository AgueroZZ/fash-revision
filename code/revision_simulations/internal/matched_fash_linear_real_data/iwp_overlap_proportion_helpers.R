# Helpers for directional overlap summaries used by the internal page.

directional_overlap_summary_matched <- function(reference_sets,
                                                comparison_sets) {
  reference_sets <- normalize_named_sets_matched(reference_sets)
  comparison_sets <- normalize_named_sets_matched(comparison_sets)
  expected_units <- c("Gene-variant pairs", "Genes", "Variants")
  if (!identical(names(reference_sets), expected_units) ||
      !identical(names(comparison_sets), expected_units)) {
    stop("Directional overlap sets must use the three reporting levels.")
  }

  rows <- lapply(expected_units, function(unit) {
    reference <- reference_sets[[unit]]
    comparison <- comparison_sets[[unit]]
    intersection_count <- length(intersect(reference, comparison))
    data.frame(
      unit = unit,
      reference_count = length(reference),
      comparison_count = length(comparison),
      intersection_count = intersection_count,
      comparison_covered_by_reference = intersection_count / length(comparison),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
