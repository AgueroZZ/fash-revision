make_centered_linear_deviation <- function(time_grid,
                                           amplitude = 2,
                                           direction = 1) {
  time_grid <- as.numeric(time_grid)
  if (length(time_grid) < 2L ||
      any(!is.finite(time_grid)) ||
      length(unique(time_grid)) < 2L ||
      length(amplitude) != 1L ||
      !is.finite(amplitude) ||
      amplitude <= 0 ||
      length(direction) != 1L ||
      !is.finite(direction) ||
      !direction %in% c(-1, 1)) {
    stop("Invalid centered-linear deviation settings.")
  }

  centered_time <- time_grid - mean(time_grid)
  deviation <- direction * amplitude * centered_time / max(abs(centered_time))
  if (abs(mean(deviation)) > 1e-12 ||
      abs(max(abs(deviation)) - amplitude) > 1e-12) {
    stop("The centered-linear deviation failed its scale invariants.")
  }
  deviation
}

validate_base_r1_effect_simulation <- function(base_effect_sim, time_grid) {
  if (!is.list(base_effect_sim) ||
      is.null(base_effect_sim$beta_matrix) ||
      is.null(base_effect_sim$unit_info)) {
    stop("base_effect_sim must contain beta_matrix and unit_info.")
  }
  required_columns <- c(
    "unit_index", "unit_id", "variant_id", "effect_class",
    "genetic_main_effect", "scenario"
  )
  if (!is.matrix(base_effect_sim$beta_matrix) ||
      ncol(base_effect_sim$beta_matrix) != length(time_grid) ||
      nrow(base_effect_sim$beta_matrix) != nrow(base_effect_sim$unit_info) ||
      !all(required_columns %in% names(base_effect_sim$unit_info)) ||
      any(!is.finite(base_effect_sim$beta_matrix))) {
    stop("The base R1 effect simulation is incomplete or misaligned.")
  }
  dynamic <- base_effect_sim$unit_info$effect_class == "dynamic_bspline"
  if (!any(dynamic) ||
      any(!is.finite(base_effect_sim$unit_info$genetic_main_effect[dynamic]))) {
    stop("The base R1 effect simulation has invalid dynamic main effects.")
  }
  invisible(TRUE)
}

make_r1_linear_truth_scenarios <- function(base_effect_sim,
                                           time_grid,
                                           linear_amplitude = 2,
                                           linear_sign_seed,
                                           mixture_seed) {
  validate_base_r1_effect_simulation(base_effect_sim, time_grid)
  if (length(linear_sign_seed) != 1L ||
      !is.finite(linear_sign_seed) ||
      linear_sign_seed != as.integer(linear_sign_seed) ||
      length(mixture_seed) != 1L ||
      !is.finite(mixture_seed) ||
      mixture_seed != as.integer(mixture_seed)) {
    stop("linear_sign_seed and mixture_seed must be finite integers.")
  }

  dynamic_indices <- which(
    base_effect_sim$unit_info$effect_class == "dynamic_bspline"
  )
  dynamic_count <- length(dynamic_indices)
  mixture_counts <- exact_proportional_counts(
    dynamic_count,
    c(dynamic_linear = 0.90, dynamic_bspline = 0.10)
  )

  set.seed(as.integer(linear_sign_seed))
  directions <- sample(c(-1, 1), dynamic_count, replace = TRUE)
  names(directions) <- as.character(dynamic_indices)

  linear_beta <- base_effect_sim$beta_matrix
  for (index in dynamic_indices) {
    linear_beta[index, ] <-
      base_effect_sim$unit_info$genetic_main_effect[index] +
      make_centered_linear_deviation(
        time_grid = time_grid,
        amplitude = linear_amplitude,
        direction = directions[[as.character(index)]]
      )
  }

  set.seed(as.integer(mixture_seed))
  retained_bspline_indices <- sort(sample(
    dynamic_indices,
    size = unname(mixture_counts[["dynamic_bspline"]]),
    replace = FALSE
  ))
  mixed_linear_indices <- setdiff(dynamic_indices, retained_bspline_indices)

  make_scenario <- function(scenario_name,
                            beta_matrix,
                            effect_classes) {
    effect_sim <- base_effect_sim
    effect_sim$beta_matrix <- beta_matrix
    effect_sim$unit_info$effect_class <- effect_classes
    effect_sim$unit_info$scenario <- scenario_name
    effect_sim$unit_info$unit_id <- sprintf("unit_%04d", seq_along(effect_classes))
    effect_sim$settings$truth_composition <- list(
      scenario = scenario_name,
      linear_amplitude = linear_amplitude,
      linear_sign_seed = as.integer(linear_sign_seed),
      mixture_seed = as.integer(mixture_seed),
      dynamic_count = dynamic_count,
      mixture_counts = mixture_counts
    )
    effect_sim
  }

  all_linear_classes <- base_effect_sim$unit_info$effect_class
  all_linear_classes[dynamic_indices] <- "dynamic_linear"
  all_linear <- make_scenario(
    scenario_name = "r1_all_linear_dynamic_truth",
    beta_matrix = linear_beta,
    effect_classes = all_linear_classes
  )

  mixed_beta <- linear_beta
  mixed_beta[retained_bspline_indices, ] <-
    base_effect_sim$beta_matrix[retained_bspline_indices, , drop = FALSE]
  mixed_classes <- all_linear_classes
  mixed_classes[retained_bspline_indices] <- "dynamic_bspline"
  mixed <- make_scenario(
    scenario_name = "r1_linear90_bspline10_dynamic_truth",
    beta_matrix = mixed_beta,
    effect_classes = mixed_classes
  )

  membership <- do.call(rbind, lapply(
    list(all_linear = all_linear, linear90_bspline10 = mixed),
    function(effect_sim) {
      data.frame(
        scenario = effect_sim$unit_info$scenario,
        unit_index = effect_sim$unit_info$unit_index,
        variant_id = effect_sim$unit_info$variant_id,
        base_effect_class = base_effect_sim$unit_info$effect_class,
        effect_class = effect_sim$unit_info$effect_class,
        linear_direction = ifelse(
          effect_sim$unit_info$unit_index %in% dynamic_indices,
          directions[as.character(effect_sim$unit_info$unit_index)],
          NA_real_
        ),
        selected_for_mixed_bspline = effect_sim$unit_info$unit_index %in%
          retained_bspline_indices,
        retained_bspline = effect_sim$unit_info$effect_class ==
          "dynamic_bspline",
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(membership) <- NULL

  list(
    scenarios = list(
      all_linear = all_linear,
      linear90_bspline10 = mixed
    ),
    membership = membership,
    dynamic_indices = dynamic_indices,
    retained_bspline_indices = retained_bspline_indices,
    mixed_linear_indices = mixed_linear_indices,
    directions = directions,
    mixture_counts = mixture_counts,
    settings = list(
      linear_amplitude = linear_amplitude,
      linear_sign_seed = as.integer(linear_sign_seed),
      mixture_seed = as.integer(mixture_seed)
    )
  )
}
