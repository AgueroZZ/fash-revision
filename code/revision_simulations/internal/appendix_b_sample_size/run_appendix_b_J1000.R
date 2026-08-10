# REV-001: rerun the Appendix B / Figure S2 grid settings at J = 1000.
#
# The manuscript states J = 1000 while analysis/appendixB.rmd calls
# get_result_once(J = 300, ...). At J = 300 the group sizes J*rho/2 are
# fractional for 23 of the 46 rho values, so lapply(1:size, ...) silently
# truncates and the realized null proportion drifts from the plotted value.
# At J = 1000 every group size is an exact integer, so J = 1000 is the
# self-consistent design and matches the manuscript.
#
# Existing J = 300 results are retained untouched; this writes new files.

suppressMessages(library(fashr))

stopifnot(file.exists("analysis/appendixB.rmd"))

J_TARGET   <- 1000L
SEED       <- 12345L
SIGMA_VEC  <- c(0.1, 0.3, 0.5)
PENALTY    <- 1
NUM_BASIS  <- 20
NUM_CORES  <- 5
RHO_VEC    <- seq(0.05, 0.5, by = 0.01)

OUT_DIR <- "data/appendixB"
SETTINGS <- list(
  list(name    = "denser_grid",
       spacing = 0.1,
       outfile = file.path(OUT_DIR, "simulation_result_denser_grid_J1000.RData")),
  list(name    = "dense_grid",
       spacing = 0.2,
       outfile = file.path(OUT_DIR, "simulation_result_dense_grid_J1000.RData"))
)

for (s in SETTINGS) {
  if (file.exists(s$outfile)) {
    stop("Refusing to overwrite existing output: ", s$outfile)
  }
}

# Same construction as analysis/appendixB.rmd, with round() on the group sizes.
# At J = 1000 every size is an exact integer in real arithmetic, so round()
# only removes IEEE representation error such as 1000 * (1 - 0.29) =
# 709.9999999999999, which 1:size would otherwise truncate to 709.
get_one_set_of_datasets <- function(J, pho0, pho1, sigma_vec) {
  if (pho0 <= pho1) stop("pho0 must be greater than pho1")

  sizeA <- round(J * (1 - pho0))
  sizeB <- round(J * (pho0 - pho1))
  sizeC <- round(J * pho1)
  stopifnot(sizeA + sizeB + sizeC == J)

  mk <- function(n, f) if (n > 0) lapply(seq_len(n), function(i) f()) else list()

  data_A <- mk(sizeA, function()
    simulate_process(sd_poly = 1, type = "nondynamic", sd = sigma_vec, normalize = FALSE))
  data_B <- mk(sizeB, function()
    simulate_process(sd_poly = 1, type = "linear", sd = sigma_vec, normalize = FALSE))
  data_C <- mk(sizeC, function()
    simulate_process(sd_poly = 0, type = "nonlinear", sd = sigma_vec,
                     sd_fun = 5, p = 2, normalize = FALSE))

  datasets <- c(data_A, data_B, data_C)
  names(datasets) <- c(
    if (sizeA > 0) paste0("A", seq_len(sizeA)),
    if (sizeB > 0) paste0("B", seq_len(sizeB)),
    if (sizeC > 0) paste0("C", seq_len(sizeC))
  )
  datasets
}

get_result_once <- function(J, pho0, pho1, sigma_vec, grid, penalty, num_basis, num_cores) {
  datasets <- get_one_set_of_datasets(J, pho0, pho1, sigma_vec)

  fit1 <- fash(Y = "y", smooth_var = "x", S = "sd", data_list = datasets, order = 1,
               verbose = FALSE, num_cores = num_cores,
               grid = grid, num_basis = num_basis, penalty = penalty)
  fit2 <- fash(Y = "y", smooth_var = "x", S = "sd", data_list = datasets, order = 2,
               verbose = FALSE, num_cores = num_cores,
               grid = grid, num_basis = num_basis, penalty = penalty)

  data.frame(
    pi_00       = 1 - pho0,
    pi_01       = 1 - pho1,
    hat_pi_00   = fit1$prior_weights$prior_weight[1],
    hat_pi_01   = fit2$prior_weights$prior_weight[1],
    tilde_pi_00 = BF_update(fit1, plot = FALSE)$prior_weights$prior_weight[1],
    tilde_pi_01 = BF_update(fit2, plot = FALSE)$prior_weights$prior_weight[1]
  )
}

for (s in SETTINGS) {
  grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = s$spacing))))
  cat(sprintf("[%s] setting=%s  spacing=%.1f  grid=%d components  J=%d\n",
              format(Sys.time(), "%H:%M:%S"), s$name, s$spacing, length(grid), J_TARGET))

  set.seed(SEED)
  t_start <- proc.time()[3]

  result_all <- lapply(seq_along(RHO_VEC), function(i) {
    pho0 <- RHO_VEC[i]
    r <- get_result_once(J = J_TARGET, pho0 = pho0, pho1 = pho0 / 2,
                         sigma_vec = SIGMA_VEC, grid = grid,
                         penalty = PENALTY, num_basis = NUM_BASIS,
                         num_cores = NUM_CORES)
    cat(sprintf("  [%s] %s rho=%.2f (%d/%d) elapsed=%.0fs\n",
                format(Sys.time(), "%H:%M:%S"), s$name, pho0, i, length(RHO_VEC),
                proc.time()[3] - t_start))
    flush.console()
    r
  })

  result_df <- do.call(rbind, result_all)
  attr(result_df, "provenance") <- list(
    J = J_TARGET, seed = SEED, spacing = s$spacing, penalty = PENALTY,
    num_basis = NUM_BASIS, sigma_vec = SIGMA_VEC, rho = RHO_VEC,
    generated = format(Sys.time()), script = "run_appendix_b_J1000.R"
  )
  save(result_df, file = s$outfile)
  cat(sprintf("[%s] wrote %s  (%.1f min)\n\n", format(Sys.time(), "%H:%M:%S"),
              s$outfile, (proc.time()[3] - t_start) / 60))
}

cat("REV-001 J=1000 rerun complete.\n")
