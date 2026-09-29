## Main ----

#' Inverse-probability weighting
#'
#' Estimates potential-outcome means and treatment effects from normalized
#' inverse-probability-weighted outcome means. Propensity scores may be fit by
#' maximum-likelihood logit or probit, inverse-probability tilting, or
#' covariate-balancing propensity scores.
#'
#' @inheritParams teffects_estimators
#' @return A `teffects_ipw` object.
#' @references Graham, B. S., C. C. D. X. Pinto, and D. Egel (2012).
#'   "Inverse Probability Tilting for Moment Condition Models with Missing
#'   Data." *Review of Economic Studies*, 79(3), 1053--1079.
#' @references Imai, K., and M. Ratkovic (2014). "Covariate Balancing
#'   Propensity Score." *Journal of the Royal Statistical Society: Series B*,
#'   76(1), 243--263.
#' @references Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025).
#'   "Covariate Balancing and the Equivalence of Weighting and Doubly Robust
#'   Estimators of Average Treatment Effects." arXiv:2310.18563.
#' @export
ipw <- function(
  fml,
  data,
  psmethod = c("logit", "probit", "ipt", "cbps"),
  stat = c("ate", "att", "atc", "pomeans"),
  weights = NULL,
  vce = c("robust", "bootstrap", "jackknife"),
  pstolerance = 1e-5,
  psaction = c("error", "clamp", "trim"),
  ...
) {
  if ("tmodel" %in% ...names()) {
    stop("`tmodel` is not supported; use `psmethod`.", call. = FALSE)
  }
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_value(
    pstolerance,
    "numeric scalar GE{0} LT{0.5}",
    .arg_name = "pstolerance",
  )
  dreamerr::check_set_arg(psmethod, "strict match")
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_set_arg(psaction, "strict match")

  dreamerr::check_set_arg(vce, "strict match")
  dreamerr::check_set_value(
    vce,
    "strict match(robust)",
    .arg_name = "vce",
    .message = "Only `vce = \"robust\"` is implemented for `ipw()` so far."
  )
  .teffects_ipw(
    fml = fml,
    data = data,
    psmethod = psmethod,
    stat = stat,
    weights = weights,
    pstolerance = pstolerance,
    psaction = psaction,
    call = match.call()
  )
}


## Estimator ----

.teffects_ipw <- function(
  fml,
  data,
  psmethod,
  stat,
  pstolerance,
  psaction,
  weights,
  call
) {
  prepared <- .teffects_prepare(
    fml = fml,
    data = data,
    outcome_design = "none",
    treatment_design = "full",
    weights = weights
  )
  y <- prepared$outcome
  d <- prepared$treatment
  treatment_candidate <- prepared$treatment_design
  complete <- prepared$complete
  n <- prepared$n

  treated_groups <- prepared$treated_levels
  control_level <- prepared$control_level
  groups <- c(control_level, treated_groups)
  if (length(groups) < 2L) {
    stop("`ipw()` requires at least two treatment levels.", call. = FALSE)
  }
  # Fit the generalized propensity score P(D = j | X).
  repeat {
    treatment_design <- .teffects_independent_columns(treatment_candidate)
    treatment_fit <- .teffects_fit_psmethod(
      design = treatment_design,
      treatment = d,
      groups = groups,
      psmethod = psmethod,
      stat = stat
    )
    raw_probability <- treatment_fit$probability
    if (!identical(psaction, "trim")) {
      break
    }
    overlap_failure <- .teffects_overlap_check(
      raw_probability,
      pstolerance,
      error = FALSE
    )
    if (!any(overlap_failure)) {
      break
    }
    keep <- !overlap_failure
    if (length(unique(d[keep])) != length(groups)) {
      stop(
        "Trimming the propensity score removed an entire treatment level.",
        call. = FALSE
      )
    }
    complete[complete] <- keep
    treatment_candidate <- treatment_candidate[keep, , drop = FALSE]
    y <- y[keep]
    d <- d[keep]
    prepared$weights <- prepared$weights[keep]
    n <- length(y)
  }
  probability <- raw_probability
  if (identical(psaction, "error")) {
    .teffects_overlap_check(probability, pstolerance)
  } else if (identical(psaction, "clamp")) {
    probability <- .teffects_clamp_probability(probability, pstolerance)
    .teffects_overlap_check(probability, 0)
    if (pstolerance > 0) {
      treatment_fit$dlog_probability <- .teffects_clamp_dlog_probability(
        raw_probability,
        probability,
        treatment_fit$dlog_probability,
        pstolerance
      )
    }
  }

  treated_positions <- match(treated_groups, groups)

  ## Influence function ----
  # The transported weight for outcome level k and target t is
  #
  #   w_{k|t}(O_i) = Z_{ki} r_t(X_i) / pi_k(X_i),
  #
  # with r_t(X_i) = 1 and l_t(X_i) = 0 when t is the full population.
  # The Hájek estimator normalizes w by its own sample mean, the H_{k|t} of
  # `target = NULL` marks t = star.
  pi_k <- probability
  l_k <- treatment_fit$dlog_probability
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

  potential_outcome_mean <- numeric(number_of_means)
  potential_outcome_influence <- matrix(0, n, number_of_means)
  nuisance_loading <- matrix(
    NA_real_,
    nrow = ncol(treatment_fit$score),
    ncol = number_of_means
  )
  for (position in seq_len(number_of_means)) {
    mean_k <- .teffects_ipw_potential_outcome(
      outcome_levels[[position]],
      targets[[position]],
      d,
      groups,
      n,
      pi_k,
      y,
      l_k,
      prepared$weights
    )
    potential_outcome_mean[[position]] <- mean_k$theta
    potential_outcome_influence[mean_k$rows, position] <- mean_k$phi
    nuisance_loading[, position] <- mean_k$c_kt
  }

  # Add the propensity-score correction to every target. All targets share one
  # factorization of the score Jacobian J_gamma.
  potential_outcome_influence[] <- potential_outcome_influence -
    as.matrix(
      treatment_fit$score %*%
        .teffects_solve(
          t(treatment_fit$jacobian),
          nuisance_loading
        )
    )

  formatted <- if (stat %in% c("ate", "pomeans")) {
    .teffects_format_population(
      potential_outcome_mean,
      potential_outcome_influence,
      groups,
      control_level,
      stat
    )
  } else if (stat %in% c("att", "atc")) {
    .teffects_format_pairs(
      potential_outcome_mean,
      potential_outcome_influence,
      groups,
      control_level,
      treated_positions,
      stat
    )
  }
  .teffects_result(
    formatted$estimate,
    formatted$influence,
    "ipw",
    list(
      stat = stat,
      treatment_model_name = psmethod,
      control = control_level,
      treated = if (stat %in% c("ate", "att", "atc")) {
        treated_groups
      } else if (identical(stat, "pomeans")) {
        NULL
      },
      target = names(formatted$estimate),
      sample = complete,
      formula = fml,
      call = call,
      propensity_score = probability,
      propensity_score_raw = raw_probability,
      psaction = psaction,
      auxiliary = list(
        coefficients = treatment_fit$coefficients,
        vcov = stats::vcov(treatment_fit$model)
      ),
      treatment_model = treatment_fit$model
    )
  )
}

# Hajek potential-outcome mean and influence for outcome level `k` and target
# population `target` (`NULL` marks the full population).
.teffects_ipw_potential_outcome <- function(
  k,
  target,
  d,
  groups,
  n,
  pi_k,
  y,
  l_k,
  observation_weights = NULL
) {
  k_position <- match(k, groups)
  target_position <- if (is.null(target)) NULL else match(target, groups)
  rows <- d == k
  r_t <- if (is.null(target)) rep(1, n) else pi_k[, target_position]
  w <- r_t[rows] / pi_k[rows, k_position]
  if (!is.null(observation_weights)) {
    w <- w * observation_weights[rows]
  }
  w <- w / sum(w)
  theta <- sum(w * y[rows])
  w_residual <- w * (y[rows] - theta)
  l_t <- if (is.null(target)) {
    0
  } else {
    l_k[[target_position]][rows, , drop = FALSE]
  }
  list(
    # Known-propensity term on the package's scale.
    theta = theta,
    phi = n * w_residual,
    # c_{k|t}^{IPW} multiplying phi_{gamma,i}.
    c_kt = as.numeric(Matrix::crossprod(
      l_t - l_k[[k_position]][rows, , drop = FALSE],
      w_residual
    )),
    rows = rows
  )
}

#' @export
print.teffects_ipw <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Treatment-effects estimation\n")
  cat("Estimator: inverse-probability weighting\n")
  cat("Treatment model:", x$teffects$treatment_model_name, "\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
