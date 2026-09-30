#' Regression adjustment
#'
#' Estimates potential-outcome means and treatment effects by fitting a
#' separate linear outcome regression in each treatment arm and averaging its
#' predictions over the covariate distribution of the target population.
#'
#' @inheritParams teffects_estimators
#' @return A `teffects_ra` object.
#' @references Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025).
#'   "Covariate Balancing and the Equivalence of Weighting and Doubly Robust
#'   Estimators of Average Treatment Effects." arXiv:2310.18563.
#' @export
ra <- function(
  fml,
  data,
  omodel = c("linear", "logit", "probit", "poisson"),
  stat = c("ate", "att", "atc", "pomeans"),
  weights = NULL,
  vce = c("robust", "bootstrap", "jackknife"),
  ...
) {
  if ("tmodel" %in% ...names()) {
    stop("`tmodel` is not supported by `ra()`.", call. = FALSE)
  }
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_set_arg(omodel, "strict match")
  stat <- tolower(stat)
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_set_arg(vce, "strict match")
  dreamerr::check_set_value(
    omodel,
    "strict match(linear)",
    .arg_name = "omodel",
    .message = "Only `omodel = \"linear\"` is implemented for `ra()` so far."
  )
  dreamerr::check_set_value(
    vce,
    "strict match(robust)",
    .arg_name = "vce",
    .message = "Only `vce = \"robust\"` is implemented for `ra()` so far."
  )
  .teffects_ra_linear(fml, data, stat, weights, match.call())
}

# Implements ordinary regression adjustment from the vignette "Regression
# adjustment".
.teffects_ra_linear <- function(
  fml,
  data,
  stat,
  weights,
  call
) {
  ## Data and targets ----
  prepared <- .teffects_prepare(
    fml,
    data,
    treatment_design = "none",
    weights = weights
  )
  y <- prepared$outcome
  R <- prepared$outcome_design
  n <- prepared$n
  treated_groups <- prepared$treated_levels
  control_level <- prepared$control_level
  groups <- c(control_level, treated_groups)
  if (length(groups) < 2L) {
    stop("`ra()` requires at least two treatment levels.", call. = FALSE)
  }
  treated <- match(treated_groups, groups)
  d <- prepared$treatment

  ## Outcome regressions ----
  # Each arm solves a least-squares problem, giving
  # mu_k(X_i) = R_i' beta_k. The helper retains the residuals and the
  # coefficient factorization needed for phi_{beta_k,i}.
  outcome_fit <- .teffects_fit_outcome_models(R, y, d, groups, prepared$weights)
  mu_k <- outcome_fit$predictions

  ## Influence function ----
  # For a target population t on the rows `target_rows`,
  #
  #   theta_{k|t} = P_n{a_{ti} mu_k(X_i)} / q_t,
  #   phi_{k|t,i} = a_{ti} {mu_k(X_i) - theta_{k|t}} / q_t
  #                 + b_t' phi_{beta_k,i}.
  if (stat %in% c("ate", "pomeans")) {
    # t = star: the target is the full population, so a_{ti} = 1 and q_t = 1.
    target_rows <- rep(TRUE, n)
    theta <- vapply(
      seq_len(ncol(mu_k)),
      function(k) .teffects_observation_mean(mu_k[, k], prepared$weights),
      numeric(1)
    )
    theta_if <- matrix(NA_real_, n, length(groups))
    for (k in seq_along(groups)) {
      theta_if[, k] <- .teffects_ra_influence(
        k,
        target_rows,
        theta[[k]],
        outcome_fit,
        R,
        mu_k,
        n
      )
    }
    formatted <- .teffects_format_population(
      theta,
      theta_if,
      groups,
      control_level,
      stat
    )
  } else if (stat %in% c("att", "atc")) {
    # t = g: the target is the treated group for ATT and the control group for
    # ATC. Each contrast reports theta_{j|t} and theta_{0|t} in pair order.
    number_of_means <- 2L * length(treated)
    theta <- numeric(number_of_means)
    theta_if <- matrix(NA_real_, n, number_of_means)
    mean_position <- 0L
    for (j_group in treated_groups) {
      j <- match(j_group, groups)
      target <- if (identical(stat, "att")) {
        j_group
      } else if (identical(stat, "atc")) {
        control_level
      }
      target_rows <- d == target
      for (h_group in c(j_group, control_level)) {
        h <- match(h_group, groups)
        mean_position <- mean_position + 1L
        theta[[mean_position]] <- .teffects_observation_mean(
          mu_k[target_rows, h],
          prepared$weights[target_rows]
        )
        theta_if[, mean_position] <- .teffects_ra_influence(
          h,
          target_rows,
          theta[[mean_position]],
          outcome_fit,
          R,
          mu_k,
          n
        )
      }
    }
    formatted <- .teffects_format_pairs(
      theta,
      theta_if,
      groups,
      control_level,
      treated,
      stat
    )
  }
  .teffects_result(
    formatted$estimate,
    formatted$influence,
    "ra",
    list(
      stat = stat,
      control = prepared$control_level,
      treated = if (stat %in% c("att", "atc")) {
        treated_groups
      } else if (stat %in% c("ate", "pomeans")) {
        NULL
      },
      target = names(formatted$estimate),
      sample = which(prepared$complete),
      osample = .teffects_osample(
        prepared$n_original,
        prepared$complete
      ),
      formula = fml,
      call = call,
      outcome_models = outcome_fit$models,
      auxiliary = list(
        coefficients = unlist(lapply(
          outcome_fit$models,
          `[[`,
          "coefficients"
        ))
      )
    )
  )
}

# Influence function phi_{k|t,i} for a target population on `target_rows`.
.teffects_ra_influence <- function(
  k,
  target_rows,
  theta_kt,
  outcome_fit,
  R,
  mu_k,
  n
) {
  a_ti <- numeric(n)
  a_ti[target_rows] <- 1
  q_t <- mean(target_rows)
  b_t <- Matrix::colMeans(R[target_rows, , drop = FALSE])
  (a_ti / q_t) *
    (mu_k[, k] - theta_kt) +
    .teffects_outcome_influence(outcome_fit, R, k, b_t)
}

#' Inverse-probability-weighted regression adjustment
#'
#' Estimates potential-outcome means and treatment effects by fitting a
#' propensity-score-weighted outcome regression in each treatment arm.
#'
#' @inheritParams teffects_estimators
#' @return A `teffects_ipwra` object.
#' @name ipwra
#' @export
ipwra <- function(
  fml,
  data,
  treatment_fml = NULL,
  omodel = c("linear", "logit", "probit", "poisson"),
  psmethod = c("logit", "probit", "ipt", "cbps"),
  stat = c("ate", "att", "atc", "pomeans"),
  weights = NULL,
  vce = c("robust", "bootstrap", "jackknife"),
  pstolerance = 1e-5,
  ...
) {
  if ("tmodel" %in% ...names()) {
    stop("`tmodel` is not supported; use `psmethod`.", call. = FALSE)
  }
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_arg(treatment_fml, "NULL os formula")
  dreamerr::check_value(
    pstolerance,
    "numeric scalar GE{0} LT{0.5}",
    .arg_name = "pstolerance",
  )
  dreamerr::check_set_arg(omodel, "strict match")
  dreamerr::check_set_arg(psmethod, "strict match")
  stat <- tolower(stat)
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_set_arg(vce, "strict match")
  dreamerr::check_set_value(
    omodel,
    "strict match(linear)",
    .arg_name = "omodel",
    .message = "Only `omodel = \"linear\"` is implemented for `ipwra()` so far."
  )
  dreamerr::check_set_value(
    vce,
    "strict match(robust)",
    .arg_name = "vce",
    .message = "Only `vce = \"robust\"` is implemented for `ipwra()` so far."
  )
  .teffects_ipwra(
    fml = fml,
    data = data,
    treatment_fml = treatment_fml,
    psmethod = psmethod,
    stat = stat,
    weights = weights,
    pstolerance = pstolerance,
    call = match.call()
  )
}

# Implements inverse-probability-weighted regression adjustment from the
# vignette "Inverse-probability-weighted regression adjustment".
.teffects_ipwra <- function(
  fml,
  data,
  treatment_fml,
  psmethod,
  stat,
  weights,
  pstolerance,
  call
) {
  ## Data and targets ----
  prepared <- .teffects_prepare(fml, data, treatment_fml, weights = weights)
  y <- prepared$outcome
  R <- prepared$outcome_design
  d <- prepared$treatment
  n <- prepared$n
  treated_groups <- prepared$treated_levels
  control_level <- prepared$control_level
  groups <- c(control_level, treated_groups)
  if (length(groups) < 2L) {
    stop("`ipwra()` requires at least two treatment levels.", call. = FALSE)
  }
  control_position <- match(control_level, groups)
  treated_positions <- match(treated_groups, groups)
  K <- length(groups)

  ## Propensity score ----
  treatment_fit <- .teffects_fit_psmethod(
    design = prepared$treatment_design,
    treatment = d,
    groups = groups,
    psmethod = psmethod,
    stat = stat
  )
  pi_k <- treatment_fit$probability
  l_k <- treatment_fit$dlog_probability
  .teffects_overlap_check(pi_k, pstolerance)
  # Compute one propensity-score influence row per observation.
  phi_gamma <- -t(.teffects_solve(
    treatment_fit$jacobian,
    t(treatment_fit$score)
  ))

  ## Reported means ----
  # Every reported coefficient is a potential-outcome mean theta_{k|t} for an
  # outcome level k and target population t. The target is the full population
  # (t = star, coded NULL) for ATE and POMs, t = j for ATT, and t = 0 for ATC.
  if (stat %in% c("ate", "pomeans")) {
    number_of_means <- K
    levels <- seq_len(K)
    targets <- rep(list(NULL), K)
  } else if (stat %in% c("att", "atc")) {
    number_of_means <- 2L * length(treated_positions)
    levels <- integer(number_of_means)
    targets <- vector("list", number_of_means)
    position <- 0L
    for (j in treated_positions) {
      target <- if (identical(stat, "att")) {
        j
      } else if (identical(stat, "atc")) {
        control_position
      }
      for (k in c(j, control_position)) {
        position <- position + 1L
        levels[[position]] <- k
        targets[[position]] <- target
      }
    }
  }

  ## Weighted regressions and influence ----
  # Each reported mean uses its own transported weight w_{k|t} and hence its own
  # weighted regression. The coefficient influence is computed below.
  #
  #   phi_{beta_{k|t},i}^w = (A_{k|t}^w)^{-1} [ w_{k|t}(O_i) R_i
  #     {Y_i - mu_k(X_i)} + C_{k|t} phi_{gamma,i} ],
  #
  # with C_{k|t} = E[w_{k|t}(O_i) R_i {Y_i - mu_k(X_i)}
  #                   {l_t(X_i) - l_k(X_i)}'].
  theta <- numeric(number_of_means)
  phi <- matrix(0, n, number_of_means)
  outcome_models <- vector("list", number_of_means)
  beta_influence <- vector("list", number_of_means)
  for (position in seq_len(number_of_means)) {
    k <- levels[[position]]
    target <- targets[[position]]
    rows <- d == groups[[k]]
    target_rows <- if (is.null(target)) {
      rep(TRUE, n)
    } else {
      d == groups[[target]]
    }
    a_ti <- numeric(n)
    a_ti[target_rows] <- 1
    q_t <- mean(a_ti)
    # Compute the target-population regressor mean.
    b_t <- Matrix::colMeans(.teffects_row_multiply(R, a_ti)) / q_t

    r_t <- if (is.null(target)) rep(1, n) else pi_k[, target]
    w_kt <- numeric(n)
    w_kt[rows] <- r_t[rows] / pi_k[rows, k]
    model <- .teffects_lm(R[rows, , drop = FALSE], y[rows], w_kt[rows])
    mu_k <- as.numeric(R %*% model$coefficients)
    residual <- y - mu_k
    beta_score <- .teffects_row_multiply(R, w_kt * residual)
    l_t <- if (is.null(target)) 0 else l_k[[target]]
    C_kt <- as.matrix(Matrix::crossprod(beta_score, l_t - l_k[[k]])) / n
    bread_inverse <- .teffects_lm_solve(model, diag(ncol(R)))
    beta_influence[[position]] <- n *
      as.matrix(
        (beta_score + phi_gamma %*% t(C_kt)) %*% bread_inverse
      )

    # Compute the outcome mean and its influence.
    theta[[position]] <- mean(mu_k[target_rows])
    phi[, position] <- (a_ti / q_t) *
      (mu_k - theta[[position]]) +
      as.numeric(beta_influence[[position]] %*% b_t)
    outcome_models[[position]] <- model
  }

  names(outcome_models) <- if (stat %in% c("ate", "pomeans")) {
    as.character(groups)
  } else if (stat %in% c("att", "atc")) {
    unlist(lapply(treated_positions, function(j) {
      target <- if (identical(stat, "att")) {
        j
      } else if (identical(stat, "atc")) {
        control_position
      }
      c(
        sprintf("%s|%s", groups[[j]], groups[[target]]),
        sprintf("%s|%s", control_level, groups[[target]])
      )
    }))
  }

  ## Auxiliary coefficient covariance ----
  # This also determines the influence of every outcome coefficient
  # and of gamma. Stacking them in the order
  # [beta_0, ..., beta_{K-1}, gamma] gives the joint covariance reported in
  # `x$teffects$auxiliary`.
  auxiliary_order <- order(
    levels,
    vapply(
      targets,
      function(target) {
        if (is.null(target)) 0L else target
      },
      integer(1)
    )
  )
  auxiliary_influence <- do.call(
    cbind,
    c(beta_influence[auxiliary_order], list(phi_gamma))
  )
  auxiliary_vcov <- as.matrix(Matrix::crossprod(auxiliary_influence))
  auxiliary_vcov <- (auxiliary_vcov + t(auxiliary_vcov)) / (2 * n^2)
  auxiliary_coef <- c(
    unlist(lapply(auxiliary_order, function(position) {
      k <- levels[[position]]
      target <- targets[[position]]
      prefix <- if (is.null(target)) {
        sprintf("OME[%s]", groups[[k]])
      } else {
        sprintf("OME[%s|%s]", groups[[k]], groups[[target]])
      }
      value <- stats::coef(outcome_models[[position]])
      names(value) <- sprintf("%s::%s", prefix, names(value))
      value
    })),
    stats::setNames(
      treatment_fit$coefficients,
      sprintf("TME::%s", names(treatment_fit$coefficients))
    )
  )
  dimnames(auxiliary_vcov) <- list(names(auxiliary_coef), names(auxiliary_coef))

  formatted <- if (stat %in% c("ate", "pomeans")) {
    .teffects_format_population(theta, phi, groups, control_level, stat)
  } else if (stat %in% c("att", "atc")) {
    .teffects_format_pairs(
      theta,
      phi,
      groups,
      control_level,
      treated_positions,
      stat
    )
  }
  .teffects_result(
    formatted$estimate,
    formatted$influence,
    "ipwra",
    list(
      stat = stat,
      psmethod = psmethod,
      outcome_model_name = "linear",
      treatment_model_name = psmethod,
      control = control_level,
      treated = if (stat %in% c("ate", "att", "atc")) {
        treated_groups
      } else if (stat == "pomeans") {
        NULL
      },
      target = names(formatted$estimate),
      sample = which(prepared$complete),
      osample = .teffects_osample(
        prepared$n_original,
        prepared$complete
      ),
      formula = fml,
      treatment_formula = treatment_fml,
      call = call,
      propensity_score = pi_k,
      auxiliary = list(
        coefficients = auxiliary_coef,
        vcov = auxiliary_vcov
      ),
      outcome_models = outcome_models,
      treatment_model = treatment_fit$model
    )
  )
}


#' @export
print.teffects_ra <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Treatment-effects estimation\n")
  cat("Estimator: regression adjustment\n")
  cat("Outcome model: linear\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}

#' @export
print.teffects_ipwra <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Treatment-effects estimation\n")
  cat("Estimator: inverse-probability-weighted regression adjustment\n")
  cat("Propensity-score method:", x$teffects$psmethod, "\n")
  cat("Outcome model: linear\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
