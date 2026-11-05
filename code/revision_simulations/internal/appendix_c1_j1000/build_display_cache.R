# Build the author-approved deterministic posterior examples at J = 1000.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("code/revision_simulations/internal/appendix_c1_j1000/reporting.R")
suppressPackageStartupMessages(library(fashr))
invisible(require_appendix_b_fashr())
report <- c1j_load_report()
focused <- report$focused
ids <- c1j_select_examples(focused)
predictions <- observations <- list()
for (k in seq_along(ids)) {
  model_order <- if (k <= 2L) 1L else 2L
  fit <- focused$fit_bundle[[paste0("order", model_order)]]$bf
  index <- match(ids[[k]], focused$data_bundle$truth$unit_id)
  stopifnot(index %in% cumulative_lfdr_calls(fit$lfdr, 0.05)$indices)
  set.seed(20260918L + k)
  prediction <- predict(fit, index = index, smooth_var = seq(1, 16, by = 0.1), M = 3000L)
  predictions[[k]] <- cbind(panel = LETTERS[k], unit_id = ids[k], prediction)
  observations[[k]] <- cbind(panel = LETTERS[k], unit_id = ids[k], focused$data_bundle$data_list[[index]])
}
display <- list(predictions = do.call(rbind, predictions), observations = do.call(rbind, observations),
  ids = ids, seeds = 20260918L + 1:4, draws = 3000L,
  selection = "First distinct numeric ID in each prespecified truth class among BF-adjusted calls at alpha = 0.05",
  focused_sha256 = c1_sha(file.path(report$directory, "focused_seed12345.rds")))
c1_write(display, file.path(report$directory, "posterior_display.rds"))
writeLines(c1_sha(file.path(report$directory, "posterior_display.rds")),
           file.path(report$directory, "posterior_display.sha256"))
cat("Selected posterior examples:", paste(ids, collapse = ", "), "\n")
