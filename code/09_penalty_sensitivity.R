library(fashr)
num_cores <- 16
load("./results/datasets_corrected.RData")


### Original grid:
log_prec <- seq(0,10, by = 0.2)
grid <- sort(c(0, exp(-0.5*log_prec)))


fash_fit1_unpenalized <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                              num_basis = 20, order = 1, betaprec = 0,
                              pred_step = 1, penalty = 1, grid = grid,
                              num_cores = num_cores, verbose = TRUE)
save(fash_fit1_unpenalized, file = "./results/fash_fit1_unpenalized_all.RData")


fash_fit2_unpenalized <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                  num_basis = 20, order = 2, betaprec = 0,
                  pred_step = 1, penalty = 1, grid = grid,
                  num_cores = num_cores, verbose = TRUE)
save(fash_fit2_unpenalized, file = "./results/fash_fit2_unpenalized_all.RData")
