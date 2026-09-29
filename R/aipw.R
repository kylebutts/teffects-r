## Main ----

#' Augmented inverse-probability weighting
#'
#' Estimates potential-outcome means and treatment effects by combining
#' arm-specific linear outcome regressions with inverse-probability-weighted
#' residual corrections.
#'
#' @inheritParams teffects_estimators
#' @return A `teffects_aipw` object.
#' @references Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025).
#'   "Covariate Balancing and the Equivalence of Weighting and Doubly Robust
#'   Estimators of Average Treatment Effects." arXiv:2310.18563.
#' @export
aipw <- function(
  fml,
  data,
  treatment_fml = NULL,
  omodel = c("linear", "logit", "probit", "poisson"),
  psmethod = c("logit", "probit", "ipt", "cbps"),
  stat = c("ate", "att", "atc", "pomeans"),
  cme = c("ml", "nls", "wnls"),
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
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_set_arg(cme, "strict match")
  dreamerr::check_set_arg(vce, "strict match")
  dreamerr::check_set_value(
    omodel,
    "strict match(linear)",
    .arg_name = "omodel",
    .message = "Only `omodel = \"linear\"` is implemented for `aipw()` yet."
  )
  dreamerr::check_set_value(
    cme,
    "strict match(ml)",
    .arg_name = "cme",
    .message = "Only `cme = \"ml\"` is implemented for `aipw()` yet."
  )
  dreamerr::check_set_value(
    vce,
    "strict match(robust)",
    .arg_name = "vce",
    .message = "Only `vce = \"robust\"` is implemented for `aipw()` so far."
  )
  .teffects_aipw(
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


## Estimator ----

# Implements the parametric augmented IPW estimator from the vignette
# "Augmented inverse-probability weighting".
.teffects_aipw <- function(
  fml,
  data,
  treatment_fml,
  psmethod,
  stat,
  pstolerance,
  weights,
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
    stop("`aipw()` requires at least two treatment levels.", call. = FALSE)
  }
  treated_positions <- match(treated_groups, groups)

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

  ## Outcome regressions ----
  outcome_fit <- .teffects_fit_outcome_models(R, y, d, groups, prepared$weights)
  mu_k <- outcome_fit$predictions

  ## AIPW means and influence ----
  # The augmented pseudo-outcome is
  #
  #   h_{k|t}(O_i) = a_{ti} mu_k(X_i) + w_{k|t}(O_i) {Y_i - mu_k(X_i)},
  #
  # so that theta_{k|t}^{AIPW} = P_n h_{k|t} / P_n a_{ti}.
  # The two nuisance models contribute to the influence function through
  #
  #   c_{beta,k|t}^{AIPW}  = q_t^{-1} E[{a_{ti} - w_{k|t}(O_i)} R_i],
  #   c_{gamma,k|t}^{AIPW} = q_t^{-1} E[w_{k|t}(O_i) {Y_i - mu_k(X_i)}
  #                                     {l_t(X_i) - l_k(X_i)}].
  #
  # `target = NULL` marks the full population, where a_{ti} = 1, q_t = 1, and
  # l_t(X_i) = 0.
  if (stat %in% c("ate", "pomeans")) {
    # t = star: one potential-outcome mean per treatment level.
    number_of_means <- length(groups)
    targets <- rep(list(NULL), number_of_means)
    outcome_levels <- groups
  } else if (stat %in% c("att", "atc")) {
    # ATT uses t = j and ATC uses t = 0; each reports theta_{j|t} and
    # theta_{0|t} for every non-control level j.
    number_of_means <- 2L * length(treated_positions)
    targets <- vector("list", number_of_means)
    outcome_levels <- groups[rep(NA_integer_, number_of_means)]
    position <- 0L
    for (j in treated_groups) {
      target <- if (identical(stat, "att")) {
        j
      } else if (identical(stat, "atc")) {
        control_level
      }
      for (k in c(j, control_level)) {
        position <- position + 1L
        targets[[position]] <- target
        outcome_levels[[position]] <- k
      }
    }
  }

  theta <- numeric(number_of_means)
  phi <- matrix(0, n, number_of_means)
  pseudo_outcome <- matrix(0, n, number_of_means)
  c_beta <- matrix(0, ncol(R), number_of_means)
  c_gamma <- matrix(NA_real_, ncol(treatment_fit$score), number_of_means)
  for (position in seq_len(number_of_means)) {
    fit_k <- .teffects_aipw_potential_outcome(
      outcome_levels[[position]],
      targets[[position]],
      d,
      groups,
      n,
      pi_k,
      y,
      mu_k,
      l_k,
      R,
      prepared$weights
    )
    theta[[position]] <- fit_k$theta
    phi[, position] <- fit_k$phi
    pseudo_outcome[, position] <- fit_k$h_kt
    c_beta[, position] <- fit_k$c_beta
    c_gamma[, position] <- fit_k$c_gamma
  }

  # Add the outcome-model and propensity-score corrections. All
  # targets share one factorization of J_gamma.
  for (position in seq_len(number_of_means)) {
    if (any(c_beta[, position] != 0)) {
      phi[, position] <- phi[, position] +
        .teffects_outcome_influence(
          outcome_fit,
          R,
          match(outcome_levels[[position]], groups),
          c_beta[, position]
        )
    }
  }
  phi <- phi -
    as.matrix(
      treatment_fit$score %*%
        .teffects_solve(
          t(treatment_fit$jacobian),
          c_gamma
        )
    )

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
  colnames(pseudo_outcome) <- as.character(outcome_levels)
  .teffects_result(
    formatted$estimate,
    formatted$influence,
    "aipw",
    list(
      stat = stat,
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
      pseudo_outcome = pseudo_outcome,
      outcome_models = outcome_fit$models,
      treatment_model = treatment_fit$model,
      auxiliary = list(
        outcome = lapply(outcome_fit$models, stats::coef),
        treatment = stats::coef(treatment_fit$model)
      )
    )
  )
}

# AIPW pseudo-outcome, mean, and influence for outcome level `k` and target
# population `target` (`NULL` marks the full population).
.teffects_aipw_potential_outcome <- function(
  k,
  target,
  d,
  groups,
  n,
  pi_k,
  y,
  mu_k,
  l_k,
  R,
  observation_weights = NULL
) {
  k_position <- match(k, groups)
  target_position <- if (is.null(target)) NULL else match(target, groups)
  rows <- d == k
  target_rows <- if (is.null(target)) rep(TRUE, n) else d == target
  a_ti <- numeric(n)
  a_ti[target_rows] <- 1
  q_t <- mean(target_rows)
  r_t <- if (is.null(target)) rep(1, n) else pi_k[, target_position]
  w_kt <- numeric(n)
  w_kt[rows] <- r_t[rows] / pi_k[rows, k_position]
  if (!is.null(observation_weights)) {
    w_kt[rows] <- w_kt[rows] * observation_weights[rows]
  }
  residual <- y - mu_k[, k_position]
  h_kt <- a_ti * mu_k[, k_position] + w_kt * residual
  theta <- .teffects_observation_mean(h_kt, observation_weights) / q_t
  l_t <- if (is.null(target)) {
    0
  } else {
    l_k[[target_position]][rows, , drop = FALSE]
  }
  list(
    theta = theta,
    h_kt = h_kt,
    phi = (a_ti * (mu_k[, k_position] - theta) + w_kt * residual) / q_t,
    c_beta = (Matrix::colMeans(.teffects_row_multiply(R, a_ti)) -
      Matrix::colMeans(.teffects_row_multiply(R, w_kt))) /
      q_t,
    c_gamma = as.numeric(Matrix::crossprod(
      l_t - l_k[[k_position]][rows, , drop = FALSE],
      w_kt[rows] * residual[rows]
    )) /
      (n * q_t),
    rows = rows
  )
}


## Methods ----
#' @export
print.teffects_aipw <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Treatment-effects estimation\n")
  cat("Estimator: augmented inverse-probability weighting\n")
  cat("Outcome model: linear\n")
  cat("Propensity-score method:", x$teffects$treatment_model_name, "\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
