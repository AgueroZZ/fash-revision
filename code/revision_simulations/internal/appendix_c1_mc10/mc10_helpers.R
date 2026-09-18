# Ten-seed extension of the retained Appendix B / manuscript C.1 experiment.
source("code/revision_simulations/appendix_b/appendix_b_helpers.R")
source("code/revision_simulations/internal/appendix_c1_mc10/gaussian_likelihood.R")

c1_contract <- function() list(
  schema = "appendix-c1-mc10-v1",
  result_id = "appendix_c1_mc10_fashr0143_20260918",
  seeds = seq(12345L, by = 10000L, length.out = 10L),
  rho = seq(0.05, 0.50, by = 0.01),
  spacing = c(original = 0.2, denser = 0.1),
  grid_J = 1000L, focused_J = 1200L,
  grid_penalty = 1, focused_penalty = 10, num_basis = 20L,
  sigma = c(0.1, 0.3, 0.5), alpha = seq(0.005, 0.20, by = 0.005),
  interval = "pointwise replicate minimum and maximum",
  likelihood_engine = "exact Gaussian residual-contrast integration validated against TMB",
  version = APPENDIX_B_FASHR_VERSION, sha = APPENDIX_B_FASHR_REMOTE_SHA
)

c1_sha <- function(path) digest::digest(file = path, algo = "sha256")

c1_write <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, version = 3)
  stopifnot(file.rename(temporary, path))
  invisible(path)
}

c1_summary <- function(data, keys, metrics) {
  groups <- split(seq_len(nrow(data)), interaction(data[keys], drop = TRUE))
  result <- lapply(groups, function(indices) {
    x <- data[indices, , drop = FALSE]
    stopifnot(nrow(x) == 10L, length(unique(x$seed)) == 10L)
    out <- x[1L, keys, drop = FALSE]
    out$n_replicates <- nrow(x)
    for (metric in metrics) {
      values <- x[[metric]]
      stopifnot(all(is.finite(values)))
      out[[paste0("mean_", metric)]] <- mean(values)
      out[[paste0(metric, "_min")]] <- min(values)
      out[[paste0(metric, "_max")]] <- max(values)
    }
    out
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result[do.call(base::order, result[keys]), , drop = FALSE]
}

c1_grid_long <- function(wide) {
  result <- list()
  for (order in c("IWP1", "IWP2")) for (stage in c("Raw EB", "BF updated")) {
    x <- wide[c("seed", "setting", "rho_dynamic", "input_sha256")]
    x$order <- order
    x$stage <- stage
    x$true_pi0 <- wide[[paste0("true_pi0_", tolower(order))]]
    x$pi0 <- wide[[paste0(if (stage == "Raw EB") "raw" else "bf", "_pi0_", tolower(order))]]
    result[[length(result) + 1L]] <- x
  }
  do.call(rbind, result)
}

c1_validate <- function(grid, alpha, contract = c1_contract()) {
  stopifnot(nrow(grid) == 920L, nrow(alpha) == 1600L,
            setequal(unique(grid$seed), contract$seeds),
            setequal(unique(alpha$seed), contract$seeds),
            !anyDuplicated(grid[c("seed", "setting", "rho_dynamic")]),
            !anyDuplicated(alpha[c("seed", "order", "stage", "alpha")]),
            all(grid$J == contract$grid_J), all(alpha$J == contract$focused_J),
            all(grid$n_nondynamic + grid$n_linear + grid$n_nonlinear == grid$J),
            all(alpha$discoveries == alpha$true_discoveries + alpha$false_discoveries),
            all(abs(alpha$realized_fdp - alpha$false_discoveries / pmax(alpha$discoveries, 1)) < 1e-14))
  positives <- ifelse(alpha$order == "IWP1", 240L, 120L)
  stopifnot(all(abs(alpha$power - alpha$true_discoveries / positives) < 1e-14))
  a <- grid[grid$setting == "original", ]; b <- grid[grid$setting == "denser", ]
  key <- function(x) paste(x$seed, x$rho_dynamic)
  stopifnot(identical(a$input_sha256, b$input_sha256[match(key(a), key(b))]))
  probabilities <- unlist(grid[c("raw_pi0_iwp1", "bf_pi0_iwp1", "raw_pi0_iwp2", "bf_pi0_iwp2")])
  stopifnot(all(is.finite(probabilities)), all(probabilities >= 0 & probabilities <= 1))
  invisible(TRUE)
}

c1_run_seed <- function(seed, directory, configuration) {
  contract <- configuration$contract
  token <- digest::digest(configuration, algo = "sha256")
  prefix <- file.path(directory, "replicates", paste0("seed_", seed))
  logfile <- paste0(prefix, ".log")
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  log <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ",
                            paste0(..., collapse = ""), "\n", file = logfile, append = TRUE)
  output <- paste0(prefix, ".rds")
  if (file.exists(output)) {
    result <- readRDS(output)
    stopifnot(identical(result$configuration_token, token))
    return(output)
  }
  checkpoint_path <- paste0(prefix, "_checkpoint.rds")
  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    stopifnot(identical(checkpoint$configuration_token, token))
  } else {
    checkpoint <- list(configuration_token = token, seed = seed, grid_rows = list())
  }
  if (is.null(checkpoint$alpha)) {
    log("Starting focused dataset")
    set.seed(seed)
    data <- build_appendix_b_datasets(contract$focused_J, 0.2, 0.1, contract$sigma)
    fit <- c1_fit_models(data, sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2)))),
                                contract$focused_penalty, contract$num_basis, 1L)
    checkpoint$alpha <- focused_alpha_curve(fit, contract$alpha)
    checkpoint$alpha$seed <- seed
    checkpoint$alpha$J <- contract$focused_J
    checkpoint$focused_input_sha256 <- digest::digest(data, algo = "sha256")
    checkpoint$focused_warnings <- unique(c(fit$order1$warnings, fit$order2$warnings))
    # Retain the first prespecified dataset and fits for the fixed S4 examples.
    if (seed == contract$seeds[[1L]]) c1_write(list(data_bundle = data, fit_bundle = fit),
                                             file.path(directory, "focused_seed12345.rds"))
    c1_write(checkpoint, checkpoint_path)
    rm(data, fit); gc()
    log("Completed focused dataset")
  }
  if (is.null(checkpoint$rng_state)) set.seed(seed) else
    assign(".Random.seed", checkpoint$rng_state, envir = .GlobalEnv)
  next_index <- length(checkpoint$grid_rows) + 1L
  if (next_index <= length(contract$rho)) for (index in seq.int(next_index, length(contract$rho))) {
    rho <- contract$rho[[index]]
    log("Starting grid rho=", sprintf("%.2f", rho), " (", index, "/46)")
    data <- build_appendix_b_datasets(contract$grid_J, rho, rho / 2, contract$sigma)
    input_hash <- digest::digest(data, algo = "sha256")
    generator_rng <- .Random.seed
    rows <- lapply(names(contract$spacing), function(setting) {
      started <- proc.time()[["elapsed"]]
      spacing <- contract$spacing[[setting]]
      fit <- c1_fit_models(data, sort(c(0, exp(-0.5 * seq(0, 10, by = spacing)))),
                                  contract$grid_penalty, contract$num_basis, 1L)
      row <- summarize_appendix_b_grid_fit(fit, setting, spacing, rho, rho / 2,
                                           contract$grid_J, seed, proc.time()[["elapsed"]] - started)
      row$seed <- seed
      row$input_sha256 <- input_hash
      row$warning_text <- paste(unique(c(fit$order1$warnings, fit$order2$warnings)), collapse = " | ")
      row
    })
    checkpoint$grid_rows[[index]] <- do.call(rbind, rows)
    checkpoint$rng_state <- generator_rng
    c1_write(checkpoint, checkpoint_path)
    assign(".Random.seed", generator_rng, envir = .GlobalEnv)
    log("Completed grid rho=", sprintf("%.2f", rho))
    rm(data, rows); gc()
  }
  checkpoint$grid <- do.call(rbind, checkpoint$grid_rows)
  checkpoint$grid_rows <- NULL
  checkpoint$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  stopifnot(nrow(checkpoint$grid) == 92L, nrow(checkpoint$alpha) == 160L)
  c1_write(checkpoint, output)
  log("Completed seed ", seed)
  output
}
