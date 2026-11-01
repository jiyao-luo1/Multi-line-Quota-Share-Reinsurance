### European data analysis: baseline Case 1

script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L && nzchar(sub("^--file=", "", file_arg))) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg),
                                  winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

SCRIPT_DIR <- script_directory()
local_library <- file.path(SCRIPT_DIR, ".r_libs")
if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}
if (!requireNamespace("nleqslv", quietly = TRUE)) {
  stop("Package 'nleqslv' is required.")
}

load_loss_matrix <- function(data_file, seed = 9L) {
  data_environment <- new.env(parent = emptyenv())
  load(data_file, envir = data_environment)
  if (!exists("final_df", envir = data_environment, inherits = FALSE)) {
    stop("The data file must contain an object named 'final_df'.")
  }
  final_df <- get("final_df", envir = data_environment, inherits = FALSE)
  required_columns <- c("1", "2", "3")
  if (!all(required_columns %in% names(final_df))) {
    stop("final_df must contain loss columns named '1', '2', and '3'.")
  }
  X <- as.matrix(final_df[, required_columns])
  storage.mode(X) <- "double"
  if (any(!is.finite(X)) || any(X < 0)) {
    stop("Loss observations must be finite and nonnegative.")
  }
  set.seed(seed)
  X[sample.int(nrow(X)), , drop = FALSE]
}

weighted_col_means <- function(X, weights = NULL) {
  if (is.null(weights)) {
    return(colMeans(X))
  }
  drop(crossprod(weights, X)) / sum(weights)
}

weighted_average <- function(values, weights = NULL) {
  if (is.null(weights)) {
    mean(values)
  } else {
    sum(weights * values) / sum(weights)
  }
}

smooth_tail_state <- function(values, alpha, bandwidth, weights = NULL,
                              initial = NULL, tolerance = 1e-9,
                              max_iterations = 12L) {
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("The smoothing bandwidth must be positive.")
  }
  target <- 1 - alpha
  total_weight <- if (is.null(weights)) length(values) else sum(weights)

  evaluate_at <- function(x) {
    standardized <- (x - values) / bandwidth
    tail_probability <- pnorm(standardized, lower.tail = FALSE)
    kernel_density <- dnorm(standardized) / bandwidth
    if (is.null(weights)) {
      equation <- mean(tail_probability) - target
      derivative <- -mean(kernel_density)
    } else {
      equation <- sum(weights * tail_probability) / total_weight - target
      derivative <- -sum(weights * kernel_density) / total_weight
    }
    list(x = x, equation = equation, derivative = derivative,
         tail_probability = tail_probability,
         kernel_density = kernel_density)
  }

  if (is.null(initial) || !is.finite(initial)) {
    initial <- unname(quantile(values, probs = alpha, type = 8,
                               names = FALSE))
  }
  state <- evaluate_at(initial)
  converged <- abs(state$equation) <= tolerance

  for (iteration in seq_len(max_iterations)) {
    if (converged) break
    if (!is.finite(state$derivative) ||
        abs(state$derivative) < .Machine$double.eps^0.5) break
    proposal <- state$x - state$equation / state$derivative
    if (!is.finite(proposal)) break
    next_state <- evaluate_at(proposal)
    if (abs(next_state$equation) > abs(state$equation)) {
      proposal <- (proposal + state$x) / 2
      next_state <- evaluate_at(proposal)
    }
    state <- next_state
    converged <- abs(state$equation) <= tolerance
  }

  if (!converged) {
    lower <- min(values) - 10 * bandwidth
    upper <- max(values) + 10 * bandwidth
    root <- uniroot(
      function(x) {
        standardized <- (x - values) / bandwidth
        tail_probability <- pnorm(standardized, lower.tail = FALSE)
        if (is.null(weights)) {
          mean(tail_probability) - target
        } else {
          sum(weights * tail_probability) / total_weight - target
        }
      },
      interval = c(lower, upper),
      tol = tolerance
    )$root
    state <- evaluate_at(root)
  }

  state
}

joint_ceded_proportion <- function(beta_vec, tilde_beta, mu0, mu_hat) {
  (sum(beta_vec * mu_hat) + mu0) / (tilde_beta * sum(mu_hat))
}

solve_individual_proportions <- function(
    X, beta_vec, tilde_beta_vec, mu0, alpha, bandwidth,
    c_start, lambda_start = 1, weights = NULL, mu_hat = NULL) {
  d <- ncol(X)
  if (is.null(mu_hat)) mu_hat <- weighted_col_means(X, weights)
  total_weight <- if (is.null(weights)) nrow(X) else sum(weights)
  c_start <- pmin(pmax(c_start, 1e-4), 1 - 1e-4)

  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$state <- NULL
  cache$gamma <- NULL

  state_at <- function(par) {
    if (!is.null(cache$par) && identical(par, cache$par)) {
      return(cache$state)
    }
    c_vec <- par[seq_len(d)]
    lambda <- par[d + 1L]
    retained_loss <- drop(X %*% (1 - c_vec))
    quantile_state <- smooth_tail_state(
      retained_loss, alpha, bandwidth, weights, initial = cache$gamma
    )

    tail_weight <- quantile_state$tail_probability
    if (!is.null(weights)) tail_weight <- weights * tail_weight
    tail_x_mean <- drop(crossprod(tail_weight, X)) / total_weight
    equation_1 <- -tail_x_mean / (1 - alpha) + mu_hat +
      lambda * tilde_beta_vec * mu_hat
    equation_2 <- sum((c_vec * tilde_beta_vec - beta_vec) * mu_hat) - mu0

    density_weight <- quantile_state$kernel_density
    if (!is.null(weights)) density_weight <- weights * density_weight
    density_sum <- sum(density_weight)
    if (!is.finite(density_sum) || density_sum <= 0) {
      stop("Degenerate smoothed quantile derivative.")
    }
    density_x_sum <- drop(crossprod(density_weight, X))
    weighted_second_moment <- crossprod(X, X * density_weight)
    jacobian_cc <- (
      weighted_second_moment -
        tcrossprod(density_x_sum) / density_sum
    ) / ((1 - alpha) * total_weight)

    jacobian_1lambda <- tilde_beta_vec * mu_hat
    jacobian_2c <- tilde_beta_vec * mu_hat
    jacobian <- rbind(
      cbind(jacobian_cc, jacobian_1lambda),
      c(jacobian_2c, 0)
    )

    result <- list(
      c_vec = c_vec,
      lambda = lambda,
      retained_loss = retained_loss,
      tail_probability = quantile_state$tail_probability,
      gamma = quantile_state$x,
      equation = c(equation_1, equation_2),
      jacobian = jacobian
    )
    # Store an explicit copy because some nonlinear-solver internals may
    # reuse and mutate their parameter vector between function evaluations.
    cache$par <- par + 0
    cache$state <- result
    cache$gamma <- quantile_state$x
    result
  }

  objective <- function(par) state_at(par)$equation
  jacobian <- function(par) state_at(par)$jacobian
  start <- c(c_start, lambda_start)

  solve_once <- function(method) {
    nleqslv::nleqslv(
      start, objective, jacobian,
      method = method, global = "dbldog", xscalm = "auto",
      control = list(ftol = 1e-8, xtol = 1e-8, maxit = 100)
    )
  }

  solution <- tryCatch(solve_once("Newton"), error = identity)
  if (inherits(solution, "error") ||
      max(abs(objective(solution$x))) > 1e-3) {
    cache$par <- NULL
    cache$state <- NULL
    cache$gamma <- NULL
    solution <- tryCatch(solve_once("Broyden"), error = identity)
  }
  if (inherits(solution, "error")) stop(conditionMessage(solution))

  final_state <- state_at(solution$x)
  max_residual <- max(abs(final_state$equation))
  if (!is.finite(max_residual) || max_residual > 1e-3) {
    stop(sprintf("The estimating equations did not converge; residual=%g.",
                 max_residual))
  }
  if (any(final_state$c_vec <= 0) || any(final_state$c_vec >= 1)) {
    stop("An individual ceded proportion is outside (0,1).")
  }

  list(
    c_hat = final_state$c_vec,
    lambda_hat = final_state$lambda,
    convergence = solution$termcd,
    max_residual = max_residual,
    state = final_state
  )
}

estimate_diversification <- function(
    X, X_upper, S, S_lower, c_joint, c_individual, mu0, alpha,
    bandwidth, individual_state, mu_hat, weights = NULL,
    weights_lower = NULL, weights_upper = NULL, joint_quantile_initial = NULL) {
  half_bandwidth <- 2^(1 / 3) * bandwidth
  total_weight <- if (is.null(weights)) nrow(X) else sum(weights)
  lower_weight <- if (is.null(weights_lower)) length(S_lower) else
    sum(weights_lower)
  upper_weight <- if (is.null(weights_upper)) nrow(X_upper) else
    sum(weights_upper)

  retained_upper <- drop(X_upper %*% (1 - c_individual))
  individual_upper_state <- smooth_tail_state(
    retained_upper, alpha, half_bandwidth, weights_upper,
    initial = individual_state$gamma
  )
  joint_full_state <- smooth_tail_state(
    S, alpha, bandwidth, weights, initial = joint_quantile_initial
  )
  joint_lower_state <- smooth_tail_state(
    S_lower, alpha, half_bandwidth, weights_lower,
    initial = joint_full_state$x
  )

  full_individual_weight <- individual_state$tail_probability
  if (!is.null(weights)) {
    full_individual_weight <- weights * full_individual_weight
  }
  upper_individual_weight <- individual_upper_state$tail_probability
  if (!is.null(weights_upper)) {
    upper_individual_weight <- weights_upper * upper_individual_weight
  }
  full_joint_weight <- joint_full_state$tail_probability
  if (!is.null(weights)) full_joint_weight <- weights * full_joint_weight
  lower_joint_weight <- joint_lower_state$tail_probability
  if (!is.null(weights_lower)) {
    lower_joint_weight <- weights_lower * lower_joint_weight
  }

  risk_individual_1 <- sum(
    full_individual_weight * individual_state$retained_loss
  ) / ((1 - alpha) * total_weight) -
    sum((1 - c_individual) * mu_hat) + mu0
  risk_individual_2 <- sum(
    upper_individual_weight * retained_upper
  ) / ((1 - alpha) * upper_weight) -
    sum((1 - c_individual) * mu_hat) + mu0
  risk_joint_1 <- (1 - c_joint) * sum(
    full_joint_weight * S
  ) / ((1 - alpha) * total_weight) -
    (1 - c_joint) * sum(mu_hat) + mu0
  risk_joint_2 <- (1 - c_joint) * sum(
    lower_joint_weight * S_lower
  ) / ((1 - alpha) * lower_weight) -
    (1 - c_joint) * sum(mu_hat) + mu0

  risks <- c(risk_individual_1, risk_individual_2,
             risk_joint_1, risk_joint_2)
  if (any(!is.finite(risks)) || any(risks <= 0)) {
    stop("The risk values used in the log-ratio must be positive.")
  }
  theta_1 <- log(risk_joint_1) - log(risk_individual_1)
  theta_2 <- log(risk_joint_2) - log(risk_individual_2)

  c(theta_es = theta_1,
    theta_es2 = theta_2,
    theta_es3 = (theta_1 + theta_2) / 2)
}

fit_configuration <- function(
    configuration, X, X_upper, S, S_lower, sd_S, joint_quantile_initial,
    weights = NULL, mu_hat = NULL, c_start = NULL, lambda_start = 1,
    fixed_bandwidth = NULL) {
  if (is.null(mu_hat)) mu_hat <- weighted_col_means(X, weights)
  c_joint <- joint_ceded_proportion(
    configuration$beta_vec, configuration$tilde_beta,
    configuration$mu0, mu_hat
  )
  if (!is.finite(c_joint) || c_joint <= 0 || c_joint >= 1) {
    stop(sprintf("Joint ceded proportion %.6f is outside (0,1).", c_joint))
  }
  bandwidth <- fixed_bandwidth
  if (is.null(bandwidth)) {
    bandwidth <- sd_S * (1 - c_joint) * nrow(X)^(-1 / 3)
  }
  if (is.null(c_start)) c_start <- rep(c_joint, ncol(X))

  individual_solution <- solve_individual_proportions(
    X = X,
    beta_vec = configuration$beta_vec,
    tilde_beta_vec = configuration$tilde_beta_vec,
    mu0 = configuration$mu0,
    alpha = configuration$alpha,
    bandwidth = bandwidth,
    c_start = c_start,
    lambda_start = lambda_start,
    weights = weights,
    mu_hat = mu_hat
  )

  m <- length(S_lower)
  weights_lower <- if (is.null(weights)) NULL else weights[seq_len(m)]
  weights_upper <- if (is.null(weights)) NULL else
    weights[m + seq_len(m)]
  theta <- estimate_diversification(
    X, X_upper, S, S_lower, c_joint, individual_solution$c_hat,
    configuration$mu0, configuration$alpha, bandwidth,
    individual_solution$state, mu_hat, weights,
    weights_lower, weights_upper, joint_quantile_initial
  )

  list(
    c_es_I = individual_solution$c_hat,
    c_es_J = c_joint,
    lambda = individual_solution$lambda_hat,
    theta_es = unname(theta["theta_es"]),
    theta_es2 = unname(theta["theta_es2"]),
    theta_es3 = unname(theta["theta_es3"]),
    bandwidth = bandwidth,
    convergence = individual_solution$convergence,
    max_residual = individual_solution$max_residual
  )
}

build_configurations <- function(settings, N_sim = 51L) {
  if (N_sim %% 2L != 1L) stop("N_sim must be odd.")
  midpoint <- (N_sim + 1L) %/% 2L
  changes <- seq(-0.01, 0.01, length.out = N_sim)
  mu0_changes <- seq(-1, 1, length.out = N_sim)
  family_names <- c(
    "tilde_beta_1", "tilde_beta_2", "tilde_beta_3",
    "tilde_beta", "mu0"
  )
  baseline <- list(
    alpha = settings$alpha,
    beta_vec = settings$beta_vec,
    tilde_beta = settings$tilde_beta,
    tilde_beta_vec = settings$tilde_beta_vec,
    mu0 = settings$mu0
  )
  configurations <- list(baseline)
  configuration_map <- matrix(
    1L, nrow = length(family_names), ncol = N_sim,
    dimnames = list(family_names, NULL)
  )
  parameter_values <- vector("list", length(family_names))
  names(parameter_values) <- family_names

  for (family in seq_along(family_names)) {
    if (family <= 3L) {
      parameter_values[[family]] <-
        settings$tilde_beta_vec[family] + changes
    } else if (family == 4L) {
      parameter_values[[family]] <- settings$tilde_beta + changes
    } else {
      parameter_values[[family]] <- settings$mu0 + mu0_changes
    }

    for (j in seq_len(N_sim)) {
      if (j == midpoint) next
      configuration <- baseline
      if (family <= 3L) {
        configuration$tilde_beta_vec[family] <-
          parameter_values[[family]][j]
      } else if (family == 4L) {
        configuration$tilde_beta <- parameter_values[[family]][j]
      } else {
        configuration$mu0 <- parameter_values[[family]][j]
      }
      configurations[[length(configurations) + 1L]] <- configuration
      configuration_map[family, j] <- length(configurations)
    }
  }

  list(
    configurations = configurations,
    map = configuration_map,
    parameter_values = parameter_values,
    family_names = family_names,
    midpoint = midpoint
  )
}

solve_all_configurations <- function(
    configurations, configuration_map, midpoint,
    X, X_upper, S, S_lower, sd_S, joint_quantile_initial,
    weights = NULL, mu_hat = NULL, point_fits = NULL) {
  number_of_configurations <- length(configurations)
  fits <- vector("list", number_of_configurations)
  if (is.null(mu_hat)) mu_hat <- weighted_col_means(X, weights)

  point_bandwidth <- function(id) {
    if (is.null(point_fits)) NULL else point_fits[[id]]$bandwidth
  }
  fallback_fit <- function(id, c_start, lambda_start) {
    tryCatch(
      fit_configuration(
        configurations[[id]], X, X_upper, S, S_lower, sd_S,
        joint_quantile_initial, weights, mu_hat, c_start, lambda_start,
        point_bandwidth(id)
      ),
      error = function(first_error) {
        if (is.null(point_fits)) stop(first_error)
        fit_configuration(
          configurations[[id]], X, X_upper, S, S_lower, sd_S,
          joint_quantile_initial, weights, mu_hat,
          point_fits[[id]]$c_es_I, point_fits[[id]]$lambda,
          point_bandwidth(id)
        )
      }
    )
  }

  baseline_start <- if (is.null(point_fits)) NULL else
    point_fits[[1L]]$c_es_I
  baseline_lambda <- if (is.null(point_fits)) 1 else
    point_fits[[1L]]$lambda
  fits[[1L]] <- fallback_fit(1L, baseline_start, baseline_lambda)

  for (family in seq_len(nrow(configuration_map))) {
    previous <- fits[[1L]]
    if (midpoint > 1L) {
      for (j in seq.int(midpoint - 1L, 1L)) {
        id <- configuration_map[family, j]
        fits[[id]] <- fallback_fit(id, previous$c_es_I, previous$lambda)
        previous <- fits[[id]]
      }
    }
    previous <- fits[[1L]]
    if (midpoint < ncol(configuration_map)) {
      for (j in seq.int(midpoint + 1L, ncol(configuration_map))) {
        id <- configuration_map[family, j]
        fits[[id]] <- fallback_fit(id, previous$c_es_I, previous$lambda)
        previous <- fits[[id]]
      }
    }
  }
  fits
}

pack_fits <- function(fits, d) {
  list(
    c_es_I = do.call(rbind, lapply(fits, `[[`, "c_es_I")),
    c_es_J = vapply(fits, `[[`, numeric(1), "c_es_J"),
    lambda = vapply(fits, `[[`, numeric(1), "lambda"),
    theta_es = vapply(fits, `[[`, numeric(1), "theta_es"),
    theta_es2 = vapply(fits, `[[`, numeric(1), "theta_es2"),
    theta_es3 = vapply(fits, `[[`, numeric(1), "theta_es3"),
    max_residual = vapply(fits, `[[`, numeric(1), "max_residual")
  )
}

bootstrap_worker <- function(bootstrap_index) {
  set.seed(worker_seed + bootstrap_index)
  weights <- rexp(worker_n)
  mu_hat <- weighted_col_means(worker_X, weights)
  fits <- solve_all_configurations(
    worker_configurations, worker_map, worker_midpoint,
    worker_X, worker_X_upper, worker_S, worker_S_lower,
    worker_sd_S, worker_joint_quantile,
    weights = weights, mu_hat = mu_hat, point_fits = worker_point_fits
  )
  pack_fits(fits, ncol(worker_X))
}

assemble_output <- function(
    point_fits, bootstrap_results, setup, settings, B, ncore) {
  point <- pack_fits(point_fits, length(settings$beta_vec))
  number_of_configurations <- nrow(point$c_es_I)
  d <- ncol(point$c_es_I)

  if (B > 0L) {
    boot_c_I <- array(
      unlist(lapply(bootstrap_results, `[[`, "c_es_I"), use.names = FALSE),
      dim = c(number_of_configurations, d, B)
    )
    boot_c_J <- do.call(cbind, lapply(bootstrap_results, `[[`, "c_es_J"))
    boot_theta_1 <- do.call(
      cbind, lapply(bootstrap_results, `[[`, "theta_es")
    )
    boot_theta_2 <- do.call(
      cbind, lapply(bootstrap_results, `[[`, "theta_es2")
    )
    boot_theta_3 <- do.call(
      cbind, lapply(bootstrap_results, `[[`, "theta_es3")
    )
  }

  output <- vector("list", length(setup$family_names))
  names(output) <- setup$family_names
  for (family in seq_along(output)) {
    output[[family]] <- vector("list", ncol(setup$map))
    for (j in seq_len(ncol(setup$map))) {
      id <- setup$map[family, j]
      c_boot <- if (B > 0L) {
        t(matrix(boot_c_I[id, , ], nrow = d, ncol = B))
      } else {
        matrix(numeric(0), nrow = 0L, ncol = d)
      }
      if (B == 1L) c_boot <- matrix(c_boot, nrow = 1L)
      output[[family]][[j]] <- list(
        c_es_I = point$c_es_I[id, ],
        c_es_J = point$c_es_J[id],
        theta_es = point$theta_es[id],
        theta_es2 = point$theta_es2[id],
        theta_es3 = point$theta_es3[id],
        c_es_I.boot = c_boot,
        c_es_J.boot = if (B > 0L) boot_c_J[id, ] else numeric(0),
        theta_es.boot = if (B > 0L) boot_theta_1[id, ] else numeric(0),
        theta_es2.boot = if (B > 0L) boot_theta_2[id, ] else numeric(0),
        theta_es3.boot = if (B > 0L) boot_theta_3[id, ] else numeric(0)
      )
    }
    attr(output[[family]], "parameter_values") <-
      setup$parameter_values[[family]]
  }
  attr(output, "settings") <- c(
    settings,
    list(B = B, ncore = ncore, N_sim = ncol(setup$map),
         sensitivity_families = setup$family_names)
  )
  output
}

compute_diagnostics <- function(output, confidence_level = 0.90) {
  point_c_I <- point_c_J <- numeric(0)
  boot_c_I <- boot_c_J <- numeric(0)
  basic_interval_values <- numeric(0)
  theta_baseline <- NULL
  midpoint <- (length(output[[1L]]) + 1L) %/% 2L

  for (family in seq_along(output)) {
    for (j in seq_along(output[[family]])) {
      result <- output[[family]][[j]]
      point_c_I <- c(point_c_I, result$c_es_I)
      point_c_J <- c(point_c_J, result$c_es_J)
      boot_c_I <- c(boot_c_I, result$c_es_I.boot)
      boot_c_J <- c(boot_c_J, result$c_es_J.boot)
      if (nrow(result$c_es_I.boot) >= 2L) {
        probabilities <- c((1 - confidence_level) / 2,
                           1 - (1 - confidence_level) / 2)
        individual_quantiles <- apply(
          result$c_es_I.boot, 2L, quantile,
          probs = probabilities, names = FALSE
        )
        joint_quantiles <- quantile(
          result$c_es_J.boot, probs = probabilities, names = FALSE
        )
        basic_interval_values <- c(
          basic_interval_values,
          2 * result$c_es_I - individual_quantiles[2L, ],
          2 * result$c_es_I - individual_quantiles[1L, ],
          2 * result$c_es_J - joint_quantiles[2L],
          2 * result$c_es_J - joint_quantiles[1L]
        )
      }
    }
    if (family == 1L) {
      baseline <- output[[family]][[midpoint]]
      theta_baseline <- c(
        theta_es = baseline$theta_es,
        theta_es2 = baseline$theta_es2,
        theta_es3 = baseline$theta_es3
      )
    }
  }

  list(
    baseline_theta = theta_baseline,
    point_c_range = range(c(point_c_I, point_c_J), finite = TRUE),
    bootstrap_c_range = if (length(boot_c_I) > 0L)
      range(c(boot_c_I, boot_c_J), finite = TRUE) else c(NA_real_, NA_real_),
    basic_90_interval_range = if (length(basic_interval_values) > 0L)
      range(basic_interval_values, finite = TRUE) else c(NA_real_, NA_real_)
  )
}

run_case_analysis <- function(X, settings, B = 14L, ncore = 14L,
                              N_sim = 51L, seed = 91000L) {
  n <- nrow(X)
  m <- floor(n / 2)
  X_upper <- X[m + seq_len(m), , drop = FALSE]
  S <- rowSums(X)
  S_lower <- S[seq_len(m)]
  sd_S <- sd(S)
  joint_quantile_initial <- unname(
    quantile(S, probs = settings$alpha, type = 8, names = FALSE)
  )
  setup <- build_configurations(settings, N_sim)

  point_fits <- solve_all_configurations(
    setup$configurations, setup$map, setup$midpoint,
    X, X_upper, S, S_lower, sd_S, joint_quantile_initial
  )

  bootstrap_results <- list()
  if (B > 0L) {
    worker_count <- max(1L, min(as.integer(ncore), as.integer(B)))
    cluster <- parallel::makePSOCKcluster(worker_count)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(
      cluster,
      function(paths) {
        .libPaths(paths)
        NULL
      },
      .libPaths()
    )

    export_names <- c(
      "weighted_col_means", "smooth_tail_state",
      "joint_ceded_proportion", "solve_individual_proportions",
      "estimate_diversification", "fit_configuration",
      "solve_all_configurations", "pack_fits", "bootstrap_worker",
      "worker_seed", "worker_n", "worker_X", "worker_X_upper",
      "worker_S", "worker_S_lower", "worker_sd_S",
      "worker_joint_quantile", "worker_configurations", "worker_map",
      "worker_midpoint", "worker_point_fits"
    )
    worker_seed <- seed
    worker_n <- n
    worker_X <- X
    worker_X_upper <- X_upper
    worker_S <- S
    worker_S_lower <- S_lower
    worker_sd_S <- sd_S
    worker_joint_quantile <- joint_quantile_initial
    worker_configurations <- setup$configurations
    worker_map <- setup$map
    worker_midpoint <- setup$midpoint
    worker_point_fits <- point_fits
    parallel::clusterExport(
      cluster, export_names, envir = environment()
    )
    bootstrap_results <- parallel::parLapplyLB(
      cluster, seq_len(B), bootstrap_worker
    )
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  }

  output <- assemble_output(
    point_fits, bootstrap_results, setup, settings, B, ncore
  )
  attr(output, "diagnostics") <- compute_diagnostics(output)
  output
}

main <- function() {
  # Production bootstrap default requested for the revised data analysis.
  B <- as.integer(Sys.getenv("REINSURANCE_B", unset = "1000"))
  ncore <- as.integer(Sys.getenv("REINSURANCE_NCORE", unset = "14"))
  N_sim <- as.integer(Sys.getenv("REINSURANCE_NSIM", unset = "51"))
  detected_cores <- parallel::detectCores(logical = TRUE)
  if (is.finite(detected_cores)) ncore <- min(ncore, detected_cores)

  settings <- list(
    alpha = 0.95,
    beta_vec = c(0.04, 0.03, 0.02),
    tilde_beta = 0.05341,
    tilde_beta_vec = c(0.07, 0.06, 0.04),
    mu0 = -1
  )

  X <- load_loss_matrix(file.path(SCRIPT_DIR, "data_eur.Rda"), seed = 9L)
  start_time <- proc.time()[["elapsed"]]
  sim.list.1 <- run_case_analysis(
    X, settings, B = B, ncore = ncore, N_sim = N_sim, seed = 91000L
  )
  elapsed <- proc.time()[["elapsed"]] - start_time
  attr(sim.list.1, "elapsed_seconds") <- elapsed
  save(sim.list.1, file = file.path(SCRIPT_DIR, "real_eur1_R1.Rda"))

  diagnostics <- attr(sim.list.1, "diagnostics")
  cat("Case 1 completed in", round(elapsed, 1), "seconds.\n")
  cat("Baseline theta estimates:",
      paste(round(diagnostics$baseline_theta, 6), collapse = ", "), "\n")
  cat("Point-estimate c range:",
      paste(round(diagnostics$point_c_range, 6), collapse = " to "), "\n")
  cat("Bootstrap c range:",
      paste(round(diagnostics$bootstrap_c_range, 6), collapse = " to "), "\n")
  cat("Basic 90% interval range:",
      paste(round(diagnostics$basic_90_interval_range, 6),
            collapse = " to "), "\n")
}

if (Sys.getenv("REINSURANCE_SKIP_MAIN", unset = "0") != "1") {
  main()
}
