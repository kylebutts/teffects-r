.teffects_late_ipwra <- function(
  fml,
  data,
  instrument,
  instrument_name,
  instrument_fml,
  instrument_ref,
  stat,
  imodel,
  pstolerance,
  weights,
  call
) {
  prepared <- .teffects_prepare(
    fml = fml,
    data = data,
    treatment_fml = instrument_fml,
    additional_complete = .teffects_late_complete(instrument),
    weights = weights
  )
  y <- prepared$outcome
  X <- prepared$outcome_design
  keep <- prepared$complete
  n <- length(y)
  sample <- .teffects_late_sample(
    prepared,
    instrument,
    instrument_ref,
    "late_ipwra"
  )
  d01 <- sample$treatment
  z <- sample$instrument
  treatment_control <- prepared$control_level
  treatment_treated <- prepared$treated_levels[[1L]]
  instrument_ref <- sample$instrument_levels[[1L]]
  instrument_exposed <- sample$instrument_levels[[2L]]

  Q <- prepared$treatment_design
  instrument_model <- .teffects_glm(Q, z, imodel)
  gamma <- instrument_model$coefficients
  propensity <- instrument_model$fitted.values
  .teffects_overlap_check(propensity, pstolerance)

  is_late <- identical(stat, "late")
  if (is_late) {
    weights0 <- (1 - z) / (1 - propensity)
    weights1 <- z / propensity
    # Fit each nuisance model on its instrument arm. Predictions and moments
    # below still use the compact full-sample design.
    outcome_models <- list(
      `0` = .teffects_lm(
        X[z == 0, , drop = FALSE],
        y[z == 0],
        weights0[z == 0]
      ),
      `1` = .teffects_lm(
        X[z == 1, , drop = FALSE],
        y[z == 1],
        weights1[z == 1]
      )
    )
    treatment_models <- list(
      `0` = .teffects_late_treatment_arm_model(
        z == 0,
        weights0,
        d01,
        X,
        n
      ),
      `1` = .teffects_late_treatment_arm_model(
        z == 1,
        weights1,
        d01,
        X,
        n
      )
    )
    beta_y <- c(
      outcome_models[[1L]]$coefficients,
      outcome_models[[2L]]$coefficients
    )
    beta_d <- c(
      treatment_models[[1L]]$coefficients,
      treatment_models[[2L]]$coefficients
    )
    mu0 <- as.numeric(X %*% outcome_models[[1L]]$coefficients)
    mu1 <- as.numeric(X %*% outcome_models[[2L]]$coefficients)
    rho0 <- treatment_models[[1L]]$fitted.values
    rho1 <- treatment_models[[2L]]$fitted.values
    reduced_form <- mean(mu1 - mu0)
    first_stage <- mean(rho1 - rho0)
  } else if (stat == "latt") {
    odds_weight <- (1 - z) * propensity / (1 - propensity)
    rows <- z == 0
    # The transported regression is identified by the reference arm only;
    # its predictions and estimating moments remain full-sample objects.
    outcome_models <- list(
      `0` = .teffects_lm(
        X[rows, , drop = FALSE],
        y[rows],
        odds_weight[rows]
      )
    )
    treatment_models <- list(
      `0` = .teffects_late_treatment_arm_model(
        z == 0,
        odds_weight,
        d01,
        X,
        n
      )
    )
    beta_y <- outcome_models[[1L]]$coefficients
    beta_d <- treatment_models[[1L]]$coefficients
    target <- z == 1
    mu0 <- as.numeric(X %*% beta_y)
    rho0 <- treatment_models[[1L]]$fitted.values
    reduced_form <- mean(y[target] - mu0[target])
    first_stage <- mean(d01[target] - rho0[target])
  }
  local_effect <- .teffects_late_effect(reduced_form, first_stage)

  p <- ncol(X)
  q <- ncol(Q)
  d0_degenerate <- !is.null(treatment_models[[1L]]$degenerate_value)
  d1_degenerate <- is_late &&
    !is.null(treatment_models[[2L]]$degenerate_value)
  family <- stats::binomial(imodel)
  theta <- c(local_effect, reduced_form, first_stage, beta_y, beta_d, gamma)
  target_names <- c(
    toupper(stat),
    if (is_late) {
      "ATE[Y|Z]"
    } else if (stat == "latt") {
      "ATT[Y|Z]"
    },
    if (is_late) {
      "ATE[D|Z]"
    } else if (stat == "latt") {
      "ATT[D|Z]"
    }
  )
  nuisance_names <- if (is_late) {
    c(
      sprintf(
        "Y[%s]::%s",
        c(rep(instrument_ref, p), rep(instrument_exposed, p)),
        rep(colnames(X), 2L)
      ),
      if (!d0_degenerate) sprintf("D[%s]::%s", instrument_ref, colnames(X)),
      if (!d1_degenerate) sprintf("D[%s]::%s", instrument_exposed, colnames(X)),
      sprintf("Z::%s", colnames(Q))
    )
  } else if (stat == "latt") {
    c(
      sprintf("Y[%s]::%s", instrument_ref, colnames(X)),
      if (!d0_degenerate) sprintf("D[%s]::%s", instrument_ref, colnames(X)),
      sprintf("Z::%s", colnames(Q))
    )
  }
  names(theta) <- c(target_names, nuisance_names)

  moment_data <- list(
    is_late = is_late,
    stat = stat,
    p = p,
    q = q,
    d0_degenerate = d0_degenerate,
    d1_degenerate = d1_degenerate,
    treatment_models = treatment_models,
    Q = Q,
    family = family,
    z = z,
    X = X,
    y = y,
    d01 = d01
  )
  moment <- .teffects_late_ipwra_moment(theta, moment_data)
  jacobian <- .teffects_numerical_moment_jacobian(
    function(parameters) {
      .teffects_late_ipwra_moment(parameters, moment_data, mean_only = TRUE)
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
  estimate <- theta[seq_len(3L)]
  colnames(influence) <- names(estimate)

  .teffects_result(
    estimate,
    influence,
    "late_ipwra",
    list(
      stat = stat,
      outcome_model_name = "linear",
      treatment_model_name = "logit",
      instrument_model_name = imodel,
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
      outcome_models = outcome_models,
      treatment_models = treatment_models,
      instrument_model = instrument_model,
      auxiliary = list(
        coefficients = theta[-seq_len(3L)],
        vcov = .teffects_influence_vcov(
          moment,
          influence_weights[, -seq_len(3L), drop = FALSE]
        )
      )
    )
  )
}

# Fit the arm-specific logistic treatment regression, predicting on the full
# sample. Returns a degenerate model when the arm has no treatment variation.
.teffects_late_treatment_arm_model <- function(rows, weights, d01, X, n) {
  observed <- d01[rows]
  if (all(observed == observed[[1L]])) {
    return(structure(
      list(
        coefficients = numeric(),
        fitted.values = rep(observed[[1L]], n),
        degenerate_value = observed[[1L]],
        vcov = matrix(numeric(), 0L, 0L)
      ),
      class = "teffects_nuisance"
    ))
  }
  model <- .teffects_glm(
    X[rows, , drop = FALSE],
    d01[rows],
    "logit",
    weights = weights[rows]
  )
  model$linear.predictors <- as.numeric(X %*% model$coefficients)
  model$fitted.values <- stats::plogis(model$linear.predictors)
  model$degenerate_value <- NULL
  model
}

# Moment function for the doubly robust LATE/LATT estimator. `moment_data`
# bundles the estimation-sample objects used by the function.
.teffects_late_ipwra_moment <- function(
  parameters,
  moment_data,
  mean_only = FALSE
) {
  is_late <- moment_data$is_late
  stat <- moment_data$stat
  p <- moment_data$p
  q <- moment_data$q
  d0_degenerate <- moment_data$d0_degenerate
  d1_degenerate <- moment_data$d1_degenerate
  treatment_models <- moment_data$treatment_models
  Q <- moment_data$Q
  family <- moment_data$family
  z <- moment_data$z
  X <- moment_data$X
  y <- moment_data$y
  d01 <- moment_data$d01

  tau <- parameters[[1L]]
  delta_y <- parameters[[2L]]
  delta_d <- parameters[[3L]]
  cursor <- 3L
  if (is_late) {
    by0 <- parameters[cursor + seq_len(p)]
    cursor <- cursor + p
    by1 <- parameters[cursor + seq_len(p)]
    cursor <- cursor + p
    if (!d0_degenerate) {
      bd0 <- parameters[cursor + seq_len(p)]
      cursor <- cursor + p
    }
    if (!d1_degenerate) {
      bd1 <- parameters[cursor + seq_len(p)]
      cursor <- cursor + p
    }
  } else if (stat == "latt") {
    by0 <- parameters[cursor + seq_len(p)]
    cursor <- cursor + p
    if (!d0_degenerate) {
      bd0 <- parameters[cursor + seq_len(p)]
      cursor <- cursor + p
    }
  }
  g <- parameters[cursor + seq_len(q)]
  eta <- as.numeric(Q %*% g)
  e <- pmin(pmax(family$linkinv(eta), 1e-10), 1 - 1e-10)
  density <- family$mu.eta(eta)
  instrument_weight <- (z - e) * density / (e * (1 - e))
  fitted_y0 <- as.numeric(X %*% by0)
  fitted_d0 <- if (d0_degenerate) {
    treatment_models[[1L]]$fitted.values
  } else {
    stats::plogis(as.numeric(X %*% bd0))
  }

  if (is_late) {
    fitted_y1 <- as.numeric(X %*% by1)
    fitted_d1 <- if (d1_degenerate) {
      treatment_models[[2L]]$fitted.values
    } else {
      stats::plogis(as.numeric(X %*% bd1))
    }
    w0 <- (1 - z) / (1 - e)
    w1 <- z / e
    if (mean_only) {
      return(c(
        delta_y - tau * delta_d,
        mean(fitted_y1 - fitted_y0 - delta_y),
        mean(fitted_d1 - fitted_d0 - delta_d),
        .teffects_weighted_colmeans(X, w0 * (y - fitted_y0)),
        .teffects_weighted_colmeans(X, w1 * (y - fitted_y1)),
        if (!d0_degenerate) {
          .teffects_weighted_colmeans(X, w0 * (d01 - fitted_d0))
        },
        if (!d1_degenerate) {
          .teffects_weighted_colmeans(X, w1 * (d01 - fitted_d1))
        },
        .teffects_weighted_colmeans(Q, instrument_weight)
      ))
    }
    moment_parts <- list(
      delta_y - tau * delta_d,
      fitted_y1 - fitted_y0 - delta_y,
      fitted_d1 - fitted_d0 - delta_d,
      .teffects_row_multiply(X, w0 * (y - fitted_y0)),
      .teffects_row_multiply(X, w1 * (y - fitted_y1))
    )
    if (!d0_degenerate) {
      moment_parts <- c(
        moment_parts,
        list(
          .teffects_row_multiply(X, w0 * (d01 - fitted_d0))
        )
      )
    }
    if (!d1_degenerate) {
      moment_parts <- c(
        moment_parts,
        list(
          .teffects_row_multiply(X, w1 * (d01 - fitted_d1))
        )
      )
    }
    moments <- do.call(
      cbind,
      c(
        moment_parts,
        list(.teffects_row_multiply(Q, instrument_weight))
      )
    )
  } else if (stat == "latt") {
    w0 <- (1 - z) * e / (1 - e)
    if (mean_only) {
      return(c(
        delta_y - tau * delta_d,
        mean(z * (y - fitted_y0 - delta_y)),
        mean(z * (d01 - fitted_d0 - delta_d)),
        .teffects_weighted_colmeans(X, w0 * (y - fitted_y0)),
        if (!d0_degenerate) {
          .teffects_weighted_colmeans(X, w0 * (d01 - fitted_d0))
        },
        .teffects_weighted_colmeans(Q, instrument_weight)
      ))
    }
    moment_parts <- list(
      delta_y - tau * delta_d,
      z * (y - fitted_y0 - delta_y),
      z * (d01 - fitted_d0 - delta_d),
      .teffects_row_multiply(X, w0 * (y - fitted_y0))
    )
    if (!d0_degenerate) {
      moment_parts <- c(
        moment_parts,
        list(
          .teffects_row_multiply(X, w0 * (d01 - fitted_d0))
        )
      )
    }
    moments <- do.call(
      cbind,
      c(
        moment_parts,
        list(.teffects_row_multiply(Q, instrument_weight))
      )
    )
  }
  as.matrix(moments)
}

.teffects_numerical_moment_jacobian <- function(mean_moment_function, theta) {
  q <- length(theta)
  jacobian <- matrix(NA_real_, q, q)
  step_scale <- .Machine$double.eps^(1 / 3)
  for (j in seq_len(q)) {
    step <- step_scale * max(1, abs(theta[[j]]))
    upper <- lower <- theta
    upper[[j]] <- upper[[j]] + step
    lower[[j]] <- lower[[j]] - step
    jacobian[, j] <- (mean_moment_function(upper) -
      mean_moment_function(lower)) /
      (2 * step)
  }
  jacobian
}

#' @export
print.teffects_late_ipwra <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Local treatment-effects estimation\n")
  cat("Estimator: doubly robust IPW regression adjustment\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n")
  cat("Outcome model: linear\n")
  cat("Treatment model: logit\n")
  cat("Instrument model:", x$teffects$instrument_model_name, "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
