# Export S3/S4 directly from the same validated chunks used by the formal page.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
destination <- arguments[[1L]]
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
lines <- readLines("analysis/appendixB.rmd", warn = FALSE)
chunks <- list(); current <- NULL; buffer <- character()
for (line in lines) {
  if (grepl("^```\\{r[ ,]", line)) {
    current <- sub("^```\\{r[ ]*([^,} ]*).*$", "\\1", line); buffer <- character()
  } else if (!is.null(current) && identical(line, "```")) {
    chunks[[current]] <- paste(buffer, collapse = "\n"); current <- NULL
  } else if (!is.null(current)) buffer <- c(buffer, line)
}
environment <- new.env(parent = globalenv())
eval(parse(text = chunks$setup), environment)
plots <- list()
for (name in c("power", "empirical-fdr", "posterior-examples")) {
  plots[[name]] <- eval(parse(text = chunks[[name]]), environment)
  stopifnot(inherits(plots[[name]], "ggplot"))
  built <- ggplot2::ggplot_build(plots[[name]])
  if (name %in% c("power", "empirical-fdr")) {
    # Every one of the 160 means must be represented by both a line and a point.
    stopifnot(nrow(built$data[[1L]]) == 160L, nrow(built$data[[2L]]) == 160L,
      all(is.finite(built$data[[1L]]$y)), all(is.finite(built$data[[2L]]$y)))
    for (panel in seq_along(built$layout$panel_params)) {
      values <- built$data[[2L]]$y[as.integer(built$data[[2L]]$PANEL) == panel]
      limits <- built$layout$panel_params[[panel]]$y.range
      stopifnot(all(values > min(limits)), all(values < max(limits)))
    }
  }
}
plots$power$data$panel <- factor(plots$power$data$panel,
  levels = levels(plots$power$data$panel), labels = c("C  Dynamic effects", "D  Nonlinear effects"))
withCallingHandlers({
  grDevices::cairo_pdf(file.path(destination, "FigureS3.pdf"), width = 10.4, height = 8.5, family = "Helvetica")
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
  print(plots[["empirical-fdr"]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plots$power, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grDevices::dev.off()
  ggplot2::ggsave(file.path(destination, "FigureS4.pdf"), plots[["posterior-examples"]],
    width = 10.4, height = 7.8, device = grDevices::cairo_pdf, family = "Helvetica")
}, warning = function(warning) stop(conditionMessage(warning), call. = FALSE))
saveRDS(plots, file.path(destination, "plots.rds"))
write.csv(environment$at005, file.path(destination, "alpha005.csv"), row.names = FALSE)
cat("Exported J=1000 Figures S3/S4; all 160 means per S3 metric are visible.\n")
