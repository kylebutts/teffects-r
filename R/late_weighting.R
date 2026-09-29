.teffects_late_resolve_instrument <- function(
  instrument,
  expression,
  data
) {
  if (
    is.character(instrument) &&
      length(instrument) == 1L &&
      instrument %in% names(data)
  ) {
    name <- instrument
    instrument <- data[[instrument]]
  } else {
    name <- paste(deparse(expression), collapse = "")
  }
  list(value = instrument, name = name)
}

.teffects_late_complete <- function(instrument) {
  complete <- !is.na(instrument)
  if (is.numeric(instrument) || is.complex(instrument)) {
    complete <- complete & is.finite(instrument)
  }
  complete
}

.teffects_late_sample <- function(
  prepared,
  instrument,
  instrument_ref,
  estimator
) {
  if (is.factor(instrument)) {
    instrument <- as.character(instrument)
  }
  if (!is.null(instrument_ref) && is.factor(instrument_ref)) {
    instrument_ref <- as.character(instrument_ref)
  }

  treatment <- prepared$treatment
  treated_levels <- prepared$treated_levels
  treatment_levels <- c(prepared$control_level, treated_levels)
  if (length(treatment_levels) != 2L) {
    stop(
      sprintf("`%s()` requires a binary actual treatment.", estimator),
      call. = FALSE
    )
  }
  treatment_ref <- prepared$control_level
  if (!(treatment_ref %in% treatment_levels)) {
    stop(
      "The treatment `ref` is not observed in the estimation sample.",
      call. = FALSE
    )
  }
  treatment_exposed <- treated_levels[[1L]]

  instrument <- instrument[prepared$complete]
  instrument_levels <- unique(instrument)
  if (length(instrument_levels) != 2L) {
    stop("`instrument` must have exactly two observed levels.", call. = FALSE)
  }
  if (is.null(instrument_ref)) {
    default_ref <- if (is.logical(instrument)) {
      FALSE
    } else if (is.numeric(instrument)) {
      0
    } else {
      instrument_levels[[1L]]
    }
    instrument_ref <- if (default_ref %in% instrument_levels) {
      default_ref
    } else {
      instrument_levels[[1L]]
    }
  }
  if (
    length(instrument_ref) != 1L ||
      is.na(instrument_ref) ||
      !(instrument_ref %in% instrument_levels)
  ) {
    stop(
      "`instrument_ref` must identify one observed instrument level.",
      call. = FALSE
    )
  }
  instrument_exposed <- instrument_levels[instrument_levels != instrument_ref][[
    1L
  ]]

  list(
    treatment = as.numeric(treatment == treatment_exposed),
    instrument = as.numeric(instrument == instrument_exposed),
    instrument_levels = c(instrument_ref, instrument_exposed)
  )
}

.teffects_late_effect <- function(reduced_form, first_stage) {
  if (
    !is.finite(first_stage) || abs(first_stage) <= sqrt(.Machine$double.eps)
  ) {
    stop(
      "The estimated instrument first stage is zero or numerically negligible.",
      call. = FALSE
    )
  }
  reduced_form / first_stage
}

.teffects_late_weighting <- function(
  estimator,
  fml,
  data,
  instrument,
  instrument_name,
  instrument_fml,
  instrument_ref,
  imodel,
  pstolerance,
  weights,
  call
) {
  prepared <- .teffects_prepare(
    fml = fml,
    data = data,
    treatment_fml = instrument_fml,
    outcome_design = "none",
    additional_complete = .teffects_late_complete(instrument),
    weights = weights
  )
  y <- prepared$outcome
  keep <- prepared$complete
  sample <- .teffects_late_sample(
    prepared,
    instrument,
    instrument_ref,
    paste0("late_", estimator)
  )
  d01 <- sample$treatment
  z <- sample$instrument
  treatment_control <- prepared$control_level
  treatment_treated <- prepared$treated_levels[[1L]]
  instrument_ref <- sample$instrument_levels[[1L]]
  instrument_exposed <- sample$instrument_levels[[2L]]

  Q <- prepared$treatment_design
  gamma <- .teffects_glm(Q, z, imodel, coefficients_only = TRUE)
  is_kappa <- identical(estimator, "kappa")
  if (!is_kappa) {
    gamma <- .teffects_fit_late_balancing(Q, z, gamma, imodel)
  }
  family <- stats::binomial(imodel)
  linear_predictor <- as.numeric(Q %*% gamma)
  propensity <- family$linkinv(linear_predictor)
  .teffects_overlap_check(propensity, pstolerance)

  w0 <- (1 - z) / (1 - propensity)
  w1 <- z / propensity
  sum_w0 <- sum(w0)
  sum_w1 <- sum(w1)
  mean_y0 <- sum(w0 * y) / sum_w0
  mean_y1 <- sum(w1 * y) / sum_w1
  mean_d0 <- sum(w0 * d01) / sum_w0
  mean_d1 <- sum(w1 * d01) / sum_w1
  reduced_form <- mean_y1 - mean_y0
  first_stage <- mean_d1 - mean_d0
  theta <- c(
    LATE = .teffects_late_effect(reduced_form, first_stage),
    `ATE[Y|Z]` = reduced_form,
    `ATE[D|Z]` = first_stage,
    `Y[0]` = mean_y0,
    `Y[1]` = mean_y1,
    `D[0]` = mean_d0,
    `D[1]` = mean_d1,
    gamma
  )
  names(theta)[-(seq_len(7L))] <- sprintf("Z::%s", colnames(Q))
  q <- ncol(Q)

  moment_data <- list(
    is_kappa = is_kappa,
    estimator = estimator,
    q = q,
    Q = Q,
    family = family,
    z = z,
    y = y,
    d01 = d01
  )
  moment <- .teffects_late_weighting_moment(theta, moment_data)
  jacobian <- .teffects_numerical_moment_jacobian(
    function(parameters) {
      .teffects_late_weighting_moment(parameters, moment_data, mean_only = TRUE)
    },
    theta
  )
  influence_weights <- .teffects_influence_weights(
    jacobian,
    seq_along(theta)
  )
  influence <- .teffects_apply_influence_weights(
    moment,
    influence_weights[, seq_len(3L), drop = FALSE]
  )
  auxiliary_vcov <- .teffects_influence_vcov(
    moment,
    influence_weights[, -seq_len(3L), drop = FALSE]
  )
  estimate <- theta[seq_len(3L)]
  colnames(influence) <- names(estimate)

  instrument_model <- structure(
    list(
      coefficients = gamma,
      linear.predictors = linear_predictor,
      fitted.values = propensity,
      family = family,
      vcov = auxiliary_vcov[4L + seq_len(q), 4L + seq_len(q), drop = FALSE]
    ),
    class = "teffects_nuisance"
  )

  .teffects_result(
    estimate,
    influence,
    paste0("late_", estimator),
    list(
      stat = "late",
      instrument_model_name = if (is_kappa) {
        imodel
      } else if (estimator == "balancing") {
        paste("balancing", imodel)
      },
      control = treatment_control,
      treated = treatment_treated,
      instrument_levels = c(instrument_ref, instrument_exposed),
      instrument = instrument_name,
      target = names(estimate),
      sample = which(keep),
      osample = .teffects_osample(prepared$n_original, keep),
      formula = fml,
      instrument_formula = instrument_fml,
      call = call,
      propensity_score = propensity,
      normalized_weights = cbind(
        control = w0 / sum_w0,
        exposed = w1 / sum_w1
      ),
      instrument_model = instrument_model,
      auxiliary = list(
        coefficients = theta[-seq_len(3L)],
        vcov = auxiliary_vcov
      )
    )
  )
}

# Moment function for the normalized weighting estimators. `moment_data`
# bundles the estimation-sample objects used by the function.
.teffects_late_weighting_moment <- function(
  parameters,
  moment_data,
  mean_only = FALSE
) {
  is_kappa <- moment_data$is_kappa
  estimator <- moment_data$estimator
  q <- moment_data$q
  Q <- moment_data$Q
  family <- moment_data$family
  z <- moment_data$z
  y <- moment_data$y
  d01 <- moment_data$d01

  tau <- parameters[[1L]]
  delta_y <- parameters[[2L]]
  delta_d <- parameters[[3L]]
  mu_y0 <- parameters[[4L]]
  mu_y1 <- parameters[[5L]]
  mu_d0 <- parameters[[6L]]
  mu_d1 <- parameters[[7L]]
  g <- parameters[7L + seq_len(q)]
  eta <- as.numeric(Q %*% g)
  e <- pmin(pmax(family$linkinv(eta), 1e-10), 1 - 1e-10)
  weight0 <- (1 - z) / (1 - e)
  weight1 <- z / e
  instrument_weight <- if (is_kappa) {
    (z - e) * family$mu.eta(eta) / (e * (1 - e))
  } else if (estimator == "balancing") {
    z / e - (1 - z) / (1 - e)
  }
  if (mean_only) {
    return(c(
      delta_y - tau * delta_d,
      mu_y1 - mu_y0 - delta_y,
      mu_d1 - mu_d0 - delta_d,
      mean(weight0 * (y - mu_y0)),
      mean(weight1 * (y - mu_y1)),
      mean(weight0 * (d01 - mu_d0)),
      mean(weight1 * (d01 - mu_d1)),
      .teffects_weighted_colmeans(Q, instrument_weight)
    ))
  }
  as.matrix(cbind(
    delta_y - tau * delta_d,
    mu_y1 - mu_y0 - delta_y,
    mu_d1 - mu_d0 - delta_d,
    weight0 * (y - mu_y0),
    weight1 * (y - mu_y1),
    weight0 * (d01 - mu_d0),
    weight1 * (d01 - mu_d1),
    .teffects_row_multiply(Q, instrument_weight)
  ))
}

# One Newton update for instrument balancing moments.
.teffects_late_balancing_update <- function(coefficients, data) {
  Q <- data$Q
  z <- data$z
  family <- data$family
  eta <- as.numeric(Q %*% coefficients)
  propensity <- pmin(pmax(family$linkinv(eta), 1e-10), 1 - 1e-10)
  density <- family$mu.eta(eta)
  weight <- z / propensity - (1 - z) / (1 - propensity)
  derivative_weight <- -density *
    (z / propensity^2 + (1 - z) / (1 - propensity)^2)
  moment <- .teffects_weighted_colmeans(Q, weight)
  jacobian <- Matrix::crossprod(
    Q,
    .teffects_row_multiply(Q, derivative_weight)
  ) /
    nrow(Q)
  list(step = -as.numeric(Matrix::solve(jacobian, moment)))
}

.teffects_fit_late_balancing <- function(Q, z, coefficients, link) {
  family <- stats::binomial(link)
  coefficients <- .teffects_newton(
    coefficients,
    .teffects_late_balancing_update,
    list(Q = Q, z = z, family = family),
    what = "instrument balancing moments"
  )$parameters
  names(coefficients) <- colnames(Q)
  coefficients
}

.print_late_weighting <- function(x, label, digits) {
  cat("Local treatment-effects estimation\n")
  cat("Estimator:", label, "\n")
  cat("Statistic: LATE\n")
  cat("Outcome model: weighted mean\n")
  cat("Treatment model: weighted mean\n")
  cat("Instrument model:", x$teffects$instrument_model_name, "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}

#' @export
print.teffects_late_kappa <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  .print_late_weighting(x, "normalized kappa", digits)
}

#' @export
print.teffects_late_balancing <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  .print_late_weighting(x, "normalized covariate balancing", digits)
}
