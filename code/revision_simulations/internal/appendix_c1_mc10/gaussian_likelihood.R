# Exact Gaussian integration for the fixed C.1 design, with a flat baseline.
# The basis, precision, EB optimizer, and BF update come from pinned fashr.
# Residual contrasts eliminate the unpenalized polynomial baseline. One
# eigendecomposition per unit then evaluates every shrinkage component.
c1_gaussian_setup <- function(time, order, num_basis = 20L, pred_step = 1) {
  template <- data.frame(x = time, y = 0, offset = 0)
  matrices <- fashr::fash_set_tmbdat(template, Si = rep(1, length(time)),
                                   num_basis = num_basis, betaprec = 0, order = order)
  X <- as.matrix(matrices$X)
  B <- as.matrix(matrices$B)
  P <- as.matrix(matrices$P)
  Q <- qr.Q(qr(X), complete = TRUE)[, seq.int(ncol(X) + 1L, nrow(X)), drop = FALSE]
  kernel <- B %*% solve(P, t(B)) *
    ((2 * order - 1) * factorial(order - 1)^2 / pred_step^(2 * order - 1))
  list(time = time, order = order, num_basis = num_basis, pred_step = pred_step,
       Q = Q, kernel = crossprod(Q, kernel %*% Q),
       baseline_logdet = as.numeric(determinant(crossprod(X), logarithm = TRUE)$modulus))
}

c1_gaussian_likelihood <- function(fash_data, grid, setup) {
  Q <- setup$Q
  m <- ncol(Q)
  matrix_rows <- lapply(seq_along(fash_data$data_list), function(index) {
    data <- fash_data$data_list[[index]]
    se <- rep(fash_data$S[[index]], length.out = nrow(data))
    stopifnot(identical(as.numeric(data$x), as.numeric(setup$time)), all(se > 0))
    noise <- crossprod(Q, Q * se^2)
    root <- chol(noise)
    whitening <- solve(t(root))
    covariance <- whitening %*% setup$kernel %*% t(whitening)
    decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
    stopifnot(min(decomposition$values) > -1e-8)
    lambda <- pmax(decomposition$values, 0)
    residual <- as.numeric(crossprod(Q, data$y - data$offset))
    projected <- as.numeric(crossprod(decomposition$vectors, whitening %*% residual))
    variance <- 1 + outer(lambda, grid^2)
    -0.5 * (m * log(2 * pi) + 2 * sum(log(diag(root))) + setup$baseline_logdet +
              colSums(log(variance)) + colSums(projected^2 / variance))
  })
  result <- do.call(rbind, matrix_rows)
  rownames(result) <- names(fash_data$data_list)
  result
}

c1_fit_from_likelihood <- function(L, data, grid, penalty, order, num_basis) {
  eb <- fashr::fash_eb_est(L, penalty = penalty, grid = grid)
  rownames(eb$posterior_weight) <- names(data$data_list)
  null <- which(eb$prior_weight$psd == 0)
  lfdr <- if (length(null)) eb$posterior_weight[, null] else rep(0, nrow(L))
  fit <- structure(list(prior_weights = eb$prior_weight,
                        posterior_weights = eb$posterior_weight,
                        psd_grid = grid, lfdr = lfdr,
                        settings = list(num_basis = num_basis, betaprec = 0,
                          order = order, pred_step = 1, likelihood = "gaussian", penalty = penalty),
                        fash_data = data, L_matrix = L, eb_result = eb), class = "fash")
  updated <- capture_warnings(fashr::BF_update(fit, plot = FALSE))
  warnings <- updated$warnings
  if (!length(null)) warnings <- unique(c(warnings,
    "The estimated prior weight of the null component (PSD = 0) is zero; lfdr is set to 0 for all datasets."))
  list(raw = fit, bf = updated$value, warnings = warnings)
}

c1_fit_models <- function(data_bundle, grid, penalty, num_basis = 20L, num_cores = 1L) {
  stopifnot(num_cores == 1L)
  data <- fashr::fash_set_data(data_list = data_bundle$data_list, Y = "y", smooth_var = "x", S = "sd")
  fits <- lapply(1:2, function(order) {
    setup <- c1_gaussian_setup(data$data_list[[1L]]$x, order, num_basis)
    L <- c1_gaussian_likelihood(data, grid, setup)
    c1_fit_from_likelihood(L, data, grid, penalty, order, num_basis)
  })
  list(order1 = fits[[1L]], order2 = fits[[2L]], truth = data_bundle$truth,
       class_counts = data_bundle$class_counts,
       settings = list(grid = grid, penalty = penalty, num_basis = num_basis, num_cores = 1L))
}
