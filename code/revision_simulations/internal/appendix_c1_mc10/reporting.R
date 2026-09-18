# Load only completed, validated ten-replicate C.1 results.
source("code/revision_simulations/internal/appendix_c1_mc10/mc10_helpers.R")

c1_load_report <- function(root = ".") {
  contract <- c1_contract()
  directory <- file.path(root, "output/revision_simulations/internal", contract$result_id)
  stopifnot(file.exists(file.path(directory, "complete.flag")))
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  stopifnot(identical(manifest$schema, contract$schema),
            identical(manifest$configuration$contract, contract))
  paths <- file.path(directory, names(manifest$artifact_sha256))
  stopifnot(identical(unname(vapply(paths, c1_sha, character(1))),
                      unname(manifest$artifact_sha256)))
  source_paths <- c(helper = "code/revision_simulations/appendix_b/appendix_b_helpers.R",
                   runner = "code/revision_simulations/internal/appendix_c1_mc10/run_appendix_c1_mc10.R",
                   extension = "code/revision_simulations/internal/appendix_c1_mc10/mc10_helpers.R",
                   integration = "code/revision_simulations/internal/appendix_c1_mc10/gaussian_likelihood.R")
  stopifnot(identical(vapply(setNames(file.path(root, source_paths), names(source_paths)), c1_sha, character(1)),
                      manifest$configuration$source_sha256))
  tables <- readRDS(file.path(directory, "summary_tables.rds"))
  c1_validate(tables$grid_by_seed, tables$alpha_by_seed, contract)
  reconstructed <- c1_summary(tables$alpha_by_seed, c("order", "stage", "alpha"),
                               c("discoveries", "realized_fdp", "power"))
  names(reconstructed) <- sub("realized_fdp", "empirical_fdr", names(reconstructed), fixed = TRUE)
  stopifnot(isTRUE(all.equal(reconstructed, tables$alpha_summary)),
            isTRUE(all.equal(c1_summary(c1_grid_long(tables$grid_by_seed),
              c("setting", "rho_dynamic", "order", "stage", "true_pi0"), "pi0"), tables$grid_summary)))
  list(directory = directory, manifest = manifest, contract = contract, tables = tables,
       focused = readRDS(file.path(directory, "focused_seed12345.rds")))
}

c1_plot_theme <- function() ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom", strip.background = ggplot2::element_rect(fill = "grey95"))

c1_load_display <- function(report) {
  path <- file.path(report$directory, "posterior_display.rds")
  stopifnot(identical(c1_sha(path), readLines(file.path(report$directory, "posterior_display.sha256"))))
  display <- readRDS(path)
  stopifnot(identical(display$focused_sha256, report$manifest$artifact_sha256[["focused_seed12345.rds"]]),
            identical(display$ids, c("C109", "B75", "C24", "C115")),
            nrow(display$predictions) == 604L, nrow(display$observations) == 64L)
  display
}

c1_plot_grid <- function(data) {
  data$stage <- factor(data$stage, levels = c("Raw EB", "BF updated"))
  panel_order <- c("original:IWP1", "original:IWP2", "denser:IWP1", "denser:IWP2")
  data$panel <- factor(paste(data$setting, data$order, sep = ":"), levels = panel_order,
    labels = c("A  IWP1: original grid", "B  IWP2: original grid",
               "C  IWP1: denser grid", "D  IWP2: denser grid"))
  stopifnot(!anyNA(data$stage), !anyNA(data$panel))
  ggplot2::ggplot(data, ggplot2::aes(true_pi0, mean_pi0, color = stage, fill = stage)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pi0_min, ymax = pi0_max),
      alpha = 0.10, color = NA, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.45) +
    ggplot2::facet_wrap(~panel, ncol = 2, scales = "free_x") +
    ggplot2::scale_color_manual(values = c("Raw EB" = "black", "BF updated" = "#0072B2")) +
    ggplot2::scale_fill_manual(values = c("Raw EB" = "black", "BF updated" = "#0072B2")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = expression("True null proportion " * pi[0]),
                  y = expression("Mean estimated null proportion " * hat(pi)[0]), color = NULL, fill = NULL) +
    c1_plot_theme()
}

c1_plot_alpha <- function(data, metric = c("power", "fdr")) {
  metric <- match.arg(metric)
  data$stage <- factor(data$stage, levels = c("Raw EB", "BF updated"))
  data$panel <- factor(data$order, levels = c("IWP1", "IWP2"),
    labels = if (metric == "fdr") c("A  Dynamic effects", "B  Nonlinear effects") else
      c("C  Dynamic effects", "D  Nonlinear effects"))
  columns <- if (metric == "power") c("mean_power", "power_min", "power_max") else
    c("mean_empirical_fdr", "empirical_fdr_min", "empirical_fdr_max")
  data$mean <- data[[columns[1]]]; data$lower <- data[[columns[2]]]; data$upper <- data[[columns[3]]]
  stopifnot(all(is.finite(as.matrix(data[c("mean", "lower", "upper")]))), !anyNA(data$stage))
  upper <- if (metric == "power") 1 else max(0.2, max(data$upper) * 1.05)
  plot <- ggplot2::ggplot(data, ggplot2::aes(alpha, mean, color = stage, fill = stage)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.10, color = NA, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0.05, color = "grey50", linetype = "dashed", linewidth = 0.4) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_color_manual(values = c("Raw EB" = "#0072B2", "BF updated" = "#D55E00")) +
    ggplot2::scale_fill_manual(values = c("Raw EB" = "#0072B2", "BF updated" = "#D55E00")) +
    ggplot2::scale_x_continuous(breaks = c(0.005, 0.05, 0.10, 0.15, 0.20)) +
    ggplot2::coord_cartesian(ylim = c(0, upper)) +
    ggplot2::labs(x = "Nominal FDR level", y = if (metric == "power") "Mean power" else "Empirical FDR",
                  color = NULL, fill = NULL) + c1_plot_theme()
  if (metric == "fdr") plot <- plot + ggplot2::geom_abline(slope = 1, intercept = 0,
    color = "grey35", linetype = "dotted", linewidth = 0.4)
  plot
}
