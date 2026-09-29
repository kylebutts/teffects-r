## Helpers ----
.teffects_osample <- function(n, sample, excluded = NULL) {
  if (is.logical(sample)) {
    if (length(sample) != n) {
      stop("A logical sample mask must have length `n`.", call. = FALSE)
    }
    result <- !sample
  } else {
    result <- rep.int(TRUE, n)
    result[sample] <- FALSE
  }
  if (!is.null(excluded)) {
    sample_rows <- if (is.logical(sample)) which(sample) else sample
    result[sample_rows[excluded]] <- TRUE
  }
  result
}

# Keep one full-rank parameterization of a design with full factor indicators.
.teffects_independent_columns <- function(x) {
  if (!ncol(x)) {
    return(x)
  }
  decomposition <- suppressWarnings(Matrix::qr(x))
  diagonal <- abs(Matrix::diag(decomposition@R))
  tolerance <- max(dim(x)) * .Machine$double.eps * max(diagonal, 1)
  keep <- decomposition@q[diagonal > tolerance] + 1L
  x[, sort(keep), drop = FALSE]
}

.teffects_solve <- function(a, b) {
  as.matrix(Matrix::solve(a, b))
}

# Newton-Raphson driver shared by the GLM and balancing fits. `step` is called
# as `step(parameters, data)` and returns `list(step =, state =)`: the update to
# add and any quantity from the final iteration the caller needs (for example
# the working `bread` for a lazy covariance). Model data is passed as a list
# rather than through `...` so the step functions stay at package scope and
# never close over a caller's frame (which would pin the whole frame, and copy
# it if the closure were ever serialized), and so no data name can be partially
# matched against a formal. Convergence is on max(abs(step)).
.teffects_newton <- function(
  parameters,
  step,
  data,
  what,
  maxit = 100L,
  tolerance = 1e-10
) {
  for (iteration in seq_len(maxit)) {
    update <- step(parameters, data)
    parameters <- parameters + update$step
    if (max(abs(update$step)) < tolerance) {
      return(list(parameters = parameters, state = update$state))
    }
  }
  stop(sprintf("The %s did not converge.", what), call. = FALSE)
}

# Implements the influence-function convention of the vignette: with
# phi_i = -J^{-1} s_i stored as rows, `positions` selects the columns of
# -moment %*% t(J^{-1}) without forming the full n-by-parameter matrix.
.teffects_influence_weights <- function(jacobian, positions) {
  selector <- matrix(0, nrow(jacobian), length(positions))
  selector[cbind(positions, seq_along(positions))] <- 1
  .teffects_solve(t(jacobian), selector)
}

# phi_i = -moment_i J^{-1}, on the package's per-observation scale.
.teffects_apply_influence_weights <- function(moment, weights) {
  -as.matrix(moment %*% weights)
}

# Var(theta_hat) = n^{-2} sum_i phi_i phi_i', the convention used throughout
# the vignette.
.teffects_influence_vcov <- function(moment, weights) {
  meat_times_weights <- Matrix::crossprod(moment) %*% weights
  as.matrix(Matrix::crossprod(weights, meat_times_weights)) / nrow(moment)^2
}

.teffects_row_multiply <- function(x, values) {
  x * as.numeric(values)
}

.teffects_weighted_colmeans <- function(x, weights) {
  as.numeric(Matrix::crossprod(x, weights)) / nrow(x)
}

.teffects_observation_mean <- function(x, weights = NULL) {
  if (is.null(weights)) {
    return(mean(x))
  }
  stats::weighted.mean(x, weights)
}


## Linear model ----
.teffects_lm_solve <- function(model, rhs) {
  as.matrix(Matrix::solve(model$factor, rhs))
}

# Fit mu_k(X_i) = R_i' beta_k in every arm by ordinary least squares. This is
# The residuals and coefficient factorizations are retained so that any target
# loading can form loading' phi_{beta_k,i}.
.teffects_fit_outcome_models <- function(
  R,
  y,
  treatment,
  groups,
  weights = NULL
) {
  n <- length(y)
  K <- length(groups)
  models <- residuals <- vector("list", K)
  mu_k <- matrix(0, n, K)

  for (k in seq_len(K)) {
    rows <- treatment == groups[[k]]
    model <- if (is.null(weights)) {
      .teffects_lm(R[rows, , drop = FALSE], y[rows])
    } else {
      .teffects_lm(R[rows, , drop = FALSE], y[rows], weights[rows])
    }
    fitted <- as.numeric(R %*% model$coefficients)
    models[[k]] <- model
    mu_k[, k] <- fitted
    residuals[[k]] <- y[rows] - fitted[rows]
  }
  names(models) <- as.character(groups)

  list(
    models = models,
    predictions = mu_k,
    residuals = residuals,
    treatment = treatment,
    groups = groups
  )
}

# Solve a weighted least-squares problem, keeping the Cholesky factorization for
# later solves.
.teffects_lm <- function(x, y, weights = rep(1, length(y))) {
  bread <- Matrix::crossprod(.teffects_row_multiply(x, sqrt(weights)))
  factor <- Matrix::Cholesky(
    Matrix::forceSymmetric(bread),
    LDL = FALSE
  )
  rhs <- Matrix::crossprod(x, weights * y)
  beta <- as.numeric(
    Matrix::solve(factor, rhs)
  )
  names(beta) <- colnames(x)

  structure(
    list(
      coefficients = beta,
      factor = factor
    ),
    class = "teffects_nuisance"
  )
}

# Evaluate loading' phi_{beta_k,i}. Only arm-k rows are nonzero.
.teffects_outcome_influence <- function(outcome_fit, R, k, loading) {
  rows <- outcome_fit$treatment == outcome_fit$groups[[k]]
  direction <- .teffects_lm_solve(outcome_fit$models[[k]], loading)
  influence <- numeric(nrow(R))
  influence[rows] <- nrow(R) *
    outcome_fit$residuals[[k]] *
    as.numeric(R[rows, , drop = FALSE] %*% direction)
  influence
}

## Generalized Linear model ----

.teffects_treatment_fit <- function(
  coefficients,
  probability,
  score,
  jacobian,
  dlog_probability,
  method
) {
  bread <- .teffects_solve(jacobian, diag(nrow(jacobian)))
  vcov <- bread %*% crossprod(score) %*% t(bread) / nrow(score)^2
  vcov <- (vcov + t(vcov)) / 2
  dimnames(vcov) <- list(names(coefficients), names(coefficients))
  structure(
    list(
      coefficients = coefficients,
      fitted.values = probability,
      probability = probability,
      score = score,
      jacobian = jacobian,
      dlog_probability = dlog_probability,
      vcov = vcov,
      method = method
    ),
    class = "teffects_nuisance"
  )
}


# One IRLS update for a binary GLM: the score, the working `bread`, and the
# Newton step at `beta`.
.teffects_glm_update <- function(beta, data) {
  x <- data$x
  y <- data$y
  weights <- data$weights
  family <- data$family
  eta <- as.numeric(x %*% beta)
  probability <- pmin(pmax(family$linkinv(eta), 1e-10), 1 - 1e-10)
  derivative <- family$mu.eta(eta)
  variance <- probability * (1 - probability)
  score <- Matrix::crossprod(
    x,
    weights * (y - probability) * derivative / variance
  )
  bread <- Matrix::crossprod(
    .teffects_row_multiply(x, sqrt(weights * derivative^2 / variance))
  )
  list(
    step = as.numeric(.teffects_solve(bread, score)),
    state = list(bread = bread)
  )
}

# IRLS for a binary GLM. Used for the binary propensity score whose score is
# score, for probit and logit links in `telasso()`, and with weights for the
# experimental LATE estimators.
.teffects_glm <- function(
  x,
  y,
  link,
  weights = rep(1, length(y)),
  maxit = 100L,
  tolerance = 1e-10,
  compute_vcov = TRUE,
  coefficients_only = FALSE
) {
  family <- stats::binomial(link = link)

  fit <- .teffects_newton(
    numeric(ncol(x)),
    .teffects_glm_update,
    list(x = x, y = y, weights = weights, family = family),
    what = "treatment model",
    maxit = maxit,
    tolerance = tolerance
  )
  beta <- fit$parameters
  names(beta) <- colnames(x)
  if (coefficients_only) {
    return(beta)
  }
  eta <- as.numeric(x %*% beta)
  probability <- family$linkinv(eta)

  structure(
    list(
      coefficients = beta,
      linear.predictors = eta,
      fitted.values = probability,
      family = family,
      vcov = if (compute_vcov) {
        .teffects_solve(fit$state$bread, diag(ncol(x)))
      } else {
        NULL
      },
      bread = if (!compute_vcov) fit$state$bread
    ),
    class = "teffects_nuisance"
  )
}

# One Newton update for the multinomial-logit score and expected information
# at `beta`.
.teffects_multinomial_update <- function(beta, data) {
  x <- data$x
  treatment <- data$treatment
  groups <- data$groups
  probability <- .teffects_multinomial_components(
    x,
    beta,
    groups,
    derivatives = FALSE
  )$probability
  fit_components <- .teffects_multinomial_diagnostics(
    x,
    treatment,
    groups,
    probability
  )
  list(
    step = as.numeric(solve(
      fit_components$information,
      fit_components$score
    ))
  )
}

# Multinomial-logit fit for a multivalued treatment. The score is the union of
# the arm-specific scores, and `dlog_probability` holds
# l_k(X_i) = d log pi_k(X_i; gamma) / d gamma for every arm.
.teffects_multinomial_logit <- function(
  x,
  treatment,
  groups,
  maxit = 100L,
  tolerance = 1e-10,
  diagnostics = TRUE,
  coefficients_only = FALSE
) {
  n <- nrow(x)
  p <- ncol(x)
  K <- length(groups)

  beta <- .teffects_newton(
    numeric((K - 1L) * p),
    .teffects_multinomial_update,
    list(x = x, treatment = treatment, groups = groups),
    what = "multinomial treatment model",
    maxit = maxit,
    tolerance = tolerance
  )$parameters

  coefficient_names <- unlist(lapply(2:K, function(j) {
    sprintf("%s::%s", groups[[j]], colnames(x))
  }))
  names(beta) <- coefficient_names
  if (coefficients_only) {
    return(beta)
  }

  components <- .teffects_multinomial_components(
    x,
    beta,
    groups,
    derivatives = diagnostics
  )
  probability <- components$probability
  if (!diagnostics) {
    return(structure(
      list(
        coefficients = beta,
        probability = probability,
        method = "multinomial logit"
      ),
      class = "teffects_nuisance"
    ))
  }

  fit_components <- .teffects_multinomial_diagnostics(
    x,
    treatment,
    groups,
    probability,
    individual = TRUE
  )
  score <- fit_components$score
  information <- fit_components$information
  colnames(score) <- coefficient_names

  .teffects_treatment_fit(
    coefficients = beta,
    probability = probability,
    score = score,
    jacobian = -information / n,
    dlog_probability = components$dlog_probability,
    method = "multinomial logit"
  )
}

# Multinomial-logit score and expected information, the ingredients of
# for a multivalued treatment model.
.teffects_multinomial_diagnostics <- function(
  x,
  treatment,
  groups,
  probability,
  individual = FALSE
) {
  n <- nrow(x)
  p <- ncol(x)
  K <- length(groups)
  number_of_coefficients <- (K - 1L) * p
  score <- if (individual) {
    matrix(0, n, number_of_coefficients)
  } else {
    numeric(number_of_coefficients)
  }
  information <- matrix(
    0,
    number_of_coefficients,
    number_of_coefficients
  )

  for (j in 2:K) {
    j_position <- (j - 2L) * p + seq_len(p)
    residual <- as.numeric(treatment == groups[[j]]) - probability[, j]
    if (individual) {
      score[, j_position] <- as.matrix(.teffects_row_multiply(x, residual))
    } else {
      score[j_position] <- as.numeric(Matrix::crossprod(x, residual))
    }
    for (h in j:K) {
      h_position <- (h - 2L) * p + seq_len(p)
      weight <- probability[, j] *
        (as.numeric(j == h) - probability[, h])
      block <- as.matrix(Matrix::crossprod(
        x,
        .teffects_row_multiply(x, weight)
      ))
      information[j_position, h_position] <- block
      if (h != j) {
        information[h_position, j_position] <- t(block)
      }
    }
  }

  list(score = score, information = information)
}

# Multinomial-logit probabilities pi_k(X_i) and their log derivatives l_k(X_i).
.teffects_multinomial_components <- function(
  x,
  coefficients,
  groups,
  derivatives = TRUE
) {
  n <- nrow(x)
  p <- ncol(x)
  K <- length(groups)
  beta <- matrix(coefficients, p, K - 1L)
  linear_predictor <- cbind(0, as.matrix(x %*% beta))
  linear_predictor <- linear_predictor - apply(linear_predictor, 1L, max)
  probability <- exp(linear_predictor)
  probability <- probability / rowSums(probability)
  colnames(probability) <- as.character(groups)
  if (!derivatives) {
    return(list(probability = probability))
  }

  dlog_probability <- vector("list", K)
  for (j in seq_len(K)) {
    derivative <- matrix(0, n, (K - 1L) * p)
    for (h in 2:K) {
      positions <- (h - 2L) * p + seq_len(p)
      derivative[, positions] <- as.matrix(.teffects_row_multiply(
        x,
        as.numeric(j == h) - probability[, h]
      ))
    }
    dlog_probability[[j]] <- derivative
  }
  list(
    probability = probability,
    dlog_probability = dlog_probability
  )
}

## Covariate-balancing propensity scores ----

# Binary covariate-balancing propensity scores. `z` is
# D_i, `pi` is pi(X_i), and `score` and `derivative` are the balancing moment
# and its derivative with respect to gamma.
#
#   ATE or POM: Q_i (D_i - pi) / {pi (1 - pi)}
#   ATT:        Q_i (D_i - pi) / (1 - pi)
#   ATC:        Q_i (D_i - pi) / pi
.teffects_cbps_binary_weights <- function(pi, z, stat) {
  complement <- 1 - pi
  if (identical(stat, "att")) {
    return(list(
      score = (z - pi) / complement,
      derivative = -(1 - z) * pi / complement
    ))
  }
  if (identical(stat, "atc")) {
    return(list(
      score = (z - pi) / pi,
      derivative = -z * complement / pi
    ))
  }
  list(
    score = (z - pi) / (pi * complement),
    derivative = -z * complement / pi - (1 - z) * pi / complement
  )
}

# Multivalued covariate-balancing propensity scores. Balancing is imposed on
# adjacent arms:
#
#   P_n[ Q_i { Z_{ki} / pi_k(X_i) - Z_{k-1,i} / pi_{k-1}(X_i) } ] = 0,
#   k = 1, ..., K - 1.
.teffects_cbps_multinomial_system <- function(
  Q,
  indicator,
  pi_k,
  keep_score = TRUE
) {
  n <- nrow(Q)
  p <- ncol(Q)
  K <- ncol(indicator)
  parameter_count <- (K - 1L) * p
  score <- if (keep_score) matrix(0, n, parameter_count) else NULL
  moment <- if (keep_score) NULL else numeric(parameter_count)
  jacobian <- matrix(0, parameter_count, parameter_count)

  for (j in 2:K) {
    equation <- (j - 2L) * p + seq_len(p)
    upper <- indicator[, j] / pi_k[, j]
    lower <- indicator[, j - 1L] / pi_k[, j - 1L]
    score_weight <- upper - lower
    if (keep_score) {
      score[, equation] <- as.matrix(.teffects_row_multiply(Q, score_weight))
    } else {
      moment[equation] <- .teffects_weighted_colmeans(Q, score_weight)
    }
    for (h in 2:K) {
      parameter <- (h - 2L) * p + seq_len(p)
      derivative <- -upper *
        (as.numeric(h == j) - pi_k[, h]) +
        lower * (as.numeric(h == j - 1L) - pi_k[, h])
      jacobian[equation, parameter] <- as.matrix(Matrix::crossprod(
        Q,
        .teffects_row_multiply(Q, derivative)
      )) /
        n
    }
  }
  list(score = score, moment = moment, jacobian = jacobian)
}

# One Newton update for binary CBPS balancing moments.
.teffects_cbps_binary_update <- function(coefficients, data) {
  Q <- data$Q
  z <- data$z
  stat <- data$stat
  n <- data$n
  pi_treated <- stats::plogis(as.numeric(Q %*% coefficients))
  weights <- .teffects_cbps_binary_weights(pi_treated, z, stat)
  moment <- .teffects_weighted_colmeans(Q, weights$score)
  jacobian <- Matrix::crossprod(
    Q,
    .teffects_row_multiply(Q, weights$derivative)
  ) /
    n
  list(step = -as.numeric(Matrix::solve(jacobian, moment)))
}

# One Newton update for multivalued CBPS balancing moments.
.teffects_cbps_multivariate_update <- function(coefficients, data) {
  Q <- data$Q
  indicator <- data$indicator
  groups <- data$groups
  pi_k <- .teffects_multinomial_components(
    Q,
    coefficients,
    groups,
    derivatives = FALSE
  )$probability
  equations <- .teffects_cbps_multinomial_system(
    Q,
    indicator,
    pi_k,
    keep_score = FALSE
  )
  list(step = -as.numeric(solve(equations$jacobian, equations$moment)))
}

.teffects_fit_cbps <- function(Q, d, groups, stat) {
  n <- nrow(Q)
  K <- length(groups)

  if (K == 2L) {
    # Binary CBPS, solved by Newton iterations on the balancing moments.
    z <- as.numeric(d == groups[[2L]])
    coefficients <- .teffects_glm(
      Q,
      z,
      "logit",
      coefficients_only = TRUE
    )
    coefficients <- .teffects_newton(
      coefficients,
      .teffects_cbps_binary_update,
      list(Q = Q, z = z, stat = stat, n = n),
      what = "CBPS balancing moments"
    )$parameters

    pi_treated <- stats::plogis(as.numeric(Q %*% coefficients))
    pi_k <- cbind(1 - pi_treated, pi_treated)
    colnames(pi_k) <- as.character(groups)
    weights <- .teffects_cbps_binary_weights(pi_treated, z, stat)
    score <- as.matrix(.teffects_row_multiply(Q, weights$score))
    jacobian <- as.matrix(Matrix::crossprod(
      Q,
      .teffects_row_multiply(Q, weights$derivative)
    )) /
      n
    # l_k(X_i) = d log pi_k(X_i; gamma) / d gamma for the two arms.
    l_k <- list(
      -as.matrix(.teffects_row_multiply(Q, pi_treated)),
      as.matrix(.teffects_row_multiply(Q, 1 - pi_treated))
    )
    names(coefficients) <- colnames(score) <- colnames(Q)
    return(.teffects_treatment_fit(
      coefficients,
      pi_k,
      score,
      jacobian,
      l_k,
      "cbps"
    ))
  }

  # Multivalued CBPS, starting from the multinomial-logit fit.
  coefficients <- .teffects_multinomial_logit(
    Q,
    d,
    groups,
    coefficients_only = TRUE
  )
  indicator <- vapply(
    groups,
    function(group) as.numeric(d == group),
    numeric(n)
  )
  coefficients <- .teffects_newton(
    coefficients,
    .teffects_cbps_multivariate_update,
    list(Q = Q, indicator = indicator, groups = groups),
    what = "multivalued CBPS balancing moments"
  )$parameters

  components <- .teffects_multinomial_components(Q, coefficients, groups)
  pi_k <- components$probability
  l_k <- components$dlog_probability
  equations <- .teffects_cbps_multinomial_system(Q, indicator, pi_k)
  score <- equations$score
  jacobian <- equations$jacobian
  colnames(score) <- names(coefficients)
  .teffects_treatment_fit(
    coefficients,
    pi_k,
    score,
    jacobian,
    l_k,
    "cbps"
  )
}

## Inverse-probability tilting ----

# One Newton update on the binary tilting moments. The moment is
# P_n[Q {Z/p - 1}] and the Jacobian is -P_n[Q Q' Z(1 - p)/p].
.teffects_ipt_target_update <- function(coefficients, data) {
  Q <- data$Q
  comparison <- data$comparison
  n <- data$n
  probability <- stats::plogis(as.numeric(Q %*% coefficients))
  score_weight <- comparison / probability - 1
  derivative_weight <- comparison * (1 - probability) / probability
  jacobian <- -as.matrix(Matrix::crossprod(
    Q,
    .teffects_row_multiply(Q, derivative_weight)
  )) /
    n
  list(
    step = -as.numeric(Matrix::solve(
      jacobian,
      .teffects_weighted_colmeans(Q, score_weight)
    ))
  )
}

# Inverse-probability tilting for a binary ATT or ATC. Fit the comparison arm's
# tilt so that
#
#   P_n[ Q_i { Z_{ki} / pi_k(X_i; gamma_k) - 1 } ] = 0,
#
# where k is the comparison arm (0 for ATT, 1 for ATC). The tilt enters
# through a logistic link pi_k = plogis(Q_i' gamma_k).
.teffects_fit_ipt_target <- function(Q, d, groups, stat) {
  n <- nrow(Q)
  comparison_position <- if (identical(stat, "att")) {
    1L
  } else if (identical(stat, "atc")) {
    2L
  }
  comparison <- as.numeric(d == groups[[comparison_position]])
  coefficients <- .teffects_glm(
    Q,
    comparison,
    "logit",
    coefficients_only = TRUE
  )

  coefficients <- .teffects_newton(
    coefficients,
    .teffects_ipt_target_update,
    list(Q = Q, comparison = comparison, n = n),
    what = "IPT balancing moments"
  )$parameters

  probability <- stats::plogis(as.numeric(Q %*% coefficients))
  score_weight <- comparison / probability - 1
  derivative_weight <- comparison * (1 - probability) / probability
  pi_k <- if (comparison_position == 1L) {
    cbind(probability, 1 - probability)
  } else {
    cbind(1 - probability, probability)
  }
  colnames(pi_k) <- as.character(groups)
  score <- as.matrix(.teffects_row_multiply(Q, score_weight))
  jacobian <- -as.matrix(Matrix::crossprod(
    Q,
    .teffects_row_multiply(Q, derivative_weight)
  )) /
    n
  # l_k(X_i) = d log pi_k(X_i; gamma) / d gamma for the two arms.
  dlog_comparison <- as.matrix(.teffects_row_multiply(
    Q,
    1 - probability
  ))
  dlog_other <- -as.matrix(.teffects_row_multiply(Q, probability))
  l_k <- if (comparison_position == 1L) {
    list(dlog_comparison, dlog_other)
  } else {
    list(dlog_other, dlog_comparison)
  }
  names(coefficients) <- colnames(score) <- sprintf(
    "IPT[%s]::%s",
    toupper(stat),
    colnames(Q)
  )

  .teffects_treatment_fit(
    coefficients,
    pi_k,
    score,
    jacobian,
    l_k,
    "ipt"
  )
}


# One Newton update for multivalued tilting moments, updating every arm
# parameter vector gamma_k in a single pass.
.teffects_ipt_update <- function(coefficients, data) {
  Q <- data$Q
  indicator <- data$indicator
  n <- data$n
  p <- data$p
  K <- data$K
  step <- matrix(0, p, K)
  for (k in seq_len(K)) {
    p_k <- stats::plogis(as.numeric(Q %*% coefficients[, k]))
    score_weight <- indicator[, k] / p_k - 1
    derivative_weight <- indicator[, k] * (1 - p_k) / p_k
    moment <- .teffects_weighted_colmeans(Q, score_weight)
    jacobian <- -Matrix::crossprod(
      Q,
      .teffects_row_multiply(Q, derivative_weight)
    ) /
      n
    step[, k] <- -as.numeric(Matrix::solve(jacobian, moment))
  }
  list(step = step)
}

# Multivalued inverse-probability tilting. Each arm has its own tilt and
# parameter vector gamma_k.
#
#   P_n[ Q_i { Z_{ki} / pi_k(X_i; gamma_k) - 1 } ] = 0,  k = 1, ..., K.
.teffects_fit_ipt <- function(Q, d, groups) {
  n <- nrow(Q)
  p <- ncol(Q)
  K <- length(groups)
  coefficients <- matrix(
    0,
    p,
    K,
    dimnames = list(colnames(Q), as.character(groups))
  )
  pi_k <- matrix(
    NA_real_,
    n,
    K,
    dimnames = list(NULL, as.character(groups))
  )
  indicator <- vapply(
    groups,
    function(group) as.numeric(d == group),
    numeric(n)
  )

  # MLE values are stable starting values for the tilting moments.
  for (k in seq_len(K)) {
    coefficients[, k] <- .teffects_glm(
      Q,
      indicator[, k],
      "logit",
      coefficients_only = TRUE
    )
  }

  coefficients <- .teffects_newton(
    coefficients,
    .teffects_ipt_update,
    list(Q = Q, indicator = indicator, n = n, p = p, K = K),
    what = "IPT balancing moments"
  )$parameters

  score <- matrix(0, n, K * p)
  jacobian <- matrix(0, K * p, K * p)
  l_k <- vector("list", K)
  for (k in seq_len(K)) {
    positions <- (k - 1L) * p + seq_len(p)
    p_k <- stats::plogis(as.numeric(Q %*% coefficients[, k]))
    score_weight <- indicator[, k] / p_k - 1
    derivative_weight <- indicator[, k] * (1 - p_k) / p_k
    pi_k[, k] <- p_k
    score[, positions] <- as.matrix(.teffects_row_multiply(
      Q,
      score_weight
    ))
    jacobian[positions, positions] <- -as.matrix(Matrix::crossprod(
      Q,
      .teffects_row_multiply(Q, derivative_weight)
    )) /
      n
    derivative <- matrix(0, n, K * p)
    derivative[, positions] <- as.matrix(.teffects_row_multiply(Q, 1 - p_k))
    l_k[[k]] <- derivative
  }
  colnames(score) <- unlist(lapply(seq_len(K), function(k) {
    sprintf("IPT[%s]::%s", groups[[k]], colnames(Q))
  }))
  coefficients_vector <- as.numeric(coefficients)
  names(coefficients_vector) <- colnames(score)

  .teffects_treatment_fit(
    coefficients_vector,
    pi_k,
    score,
    jacobian,
    l_k,
    "ipt"
  )
}

.teffects_fit_psmethod <- function(
  design,
  treatment,
  groups,
  psmethod,
  stat
) {
  if (psmethod %in% c("logit", "probit")) {
    return(.teffects_fit_generalized_treatment_model(
      design = design,
      treatment = treatment,
      groups = groups,
      tmodel = psmethod
    ))
  }

  fit <- if (identical(psmethod, "ipt")) {
    if (stat %in% c("att", "atc")) {
      if (length(groups) > 2L) {
        stop(
          paste(
            "Multivalued IPT does not support `stat = \"att\"`",
            "or `stat = \"atc\"`."
          ),
          call. = FALSE
        )
      }
      .teffects_fit_ipt_target(design, treatment, groups, stat)
    } else if (stat %in% c("ate", "pomeans")) {
      .teffects_fit_ipt(design, treatment, groups)
    }
  } else if (identical(psmethod, "cbps")) {
    .teffects_fit_cbps(design, treatment, groups, stat)
  }

  list(
    model = fit,
    coefficients = fit$coefficients,
    probability = fit$probability,
    score = fit$score,
    jacobian = fit$jacobian,
    dlog_probability = fit$dlog_probability
  )
}

# Generalized propensity score for binary and multivalued treatments. For
# binary logit or probit, the score is
#
#   s_pi(O_i; gamma) = Q_i {D_i - pi(X_i; gamma)} g_i
#                      / [pi(X_i; gamma) {1 - pi(X_i; gamma)}],
#
# with g_i = G'(Q_i' gamma). For logit this simplifies to Q_i {D_i - pi(X_i)},
# and multivalued treatments use the usual multinomial-logit score.
.teffects_fit_generalized_treatment_model <- function(
  design,
  treatment,
  groups,
  tmodel,
  prediction_only = FALSE
) {
  # `Q` and `d` are the vignette's Q_i and D_i.
  Q <- design
  d <- treatment
  if (length(groups) > 2L) {
    if (!identical(tmodel, "logit")) {
      stop(
        "Multivalued treatment models require `tmodel = \"logit\"`.",
        call. = FALSE
      )
    }
    model <- .teffects_multinomial_logit(
      Q,
      d,
      groups,
      diagnostics = !prediction_only
    )
    if (prediction_only) {
      return(list(probability = model$probability))
    }
    return(list(
      model = model,
      coefficients = model$coefficients,
      probability = model$probability,
      score = model$score,
      jacobian = model$jacobian,
      dlog_probability = model$dlog_probability
    ))
  }

  z <- as.numeric(d == groups[[2L]])
  model <- .teffects_glm(
    Q,
    z,
    tmodel,
    compute_vcov = !prediction_only
  )
  pi_treated <- model$fitted.values
  pi_k <- cbind(1 - pi_treated, pi_treated)
  colnames(pi_k) <- as.character(groups)
  if (prediction_only) {
    return(list(probability = pi_k))
  }

  g_i <- model$family$mu.eta(model$linear.predictors)
  observed_sign <- ifelse(z == 1, 1, -1)
  pi_observed <- ifelse(
    z == 1,
    pi_treated,
    1 - pi_treated
  )
  score <- as.matrix(.teffects_row_multiply(
    Q,
    observed_sign * g_i / pi_observed
  ))
  g_derivative <- if (identical(tmodel, "probit")) {
    -model$linear.predictors * g_i
  } else if (identical(tmodel, "logit")) {
    g_i * (1 - 2 * pi_treated)
  }
  score_derivative <- observed_sign *
    g_derivative /
    pi_observed -
    g_i^2 / pi_observed^2
  jacobian <- as.matrix(Matrix::crossprod(
    Q,
    .teffects_row_multiply(Q, score_derivative)
  )) /
    nrow(Q)
  # l_k(X_i) = d log pi_k(X_i; gamma) / d gamma for the two arms.
  l_k <- list(
    -as.matrix(.teffects_row_multiply(
      Q,
      g_i / (1 - pi_treated)
    )),
    as.matrix(.teffects_row_multiply(
      Q,
      g_i / pi_treated
    ))
  )
  colnames(score) <- names(model$coefficients)

  list(
    model = model,
    coefficients = model$coefficients,
    probability = pi_k,
    score = score,
    jacobian = jacobian,
    dlog_probability = l_k
  )
}


#' @export
coef.teffects_nuisance <- function(object, ...) object$coefficients

#' @export
vcov.teffects_nuisance <- function(object, ...) {
  if (!is.null(object$vcov)) {
    return(object$vcov)
  }
  if (!is.null(object$factor)) {
    return(as.matrix(
      Matrix::solve(
        object$factor,
        Matrix::Diagonal(length(object$coefficients))
      )
    ))
  }
  if (!is.null(object$bread)) {
    return(.teffects_solve(
      object$bread,
      diag(length(object$coefficients))
    ))
  }
  stop(
    "This nuisance fit does not contain covariance information.",
    call. = FALSE
  )
}
