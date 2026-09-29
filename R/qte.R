## Main ----

#' Marginal quantile treatment effects
#'
#' Estimates quantiles of marginal potential-outcome distributions and their
#' contrasts under conditional independence. The implementation follows the
#' inverse-probability-weighted (IPW) and efficient-influence-function (EIF)
#' estimators in Cattaneo (2010). The EIF estimator uses linear series
#' regressions for the conditional distribution functions.
#'
#' @inheritParams teffects_estimators
#' @param probs Numeric vector of quantile indexes strictly between zero and
#'   one.
#' @param method Estimation method. `"eif"` is the augmented, efficient-
#'   influence-function estimator and `"ipw"` is normalized inverse-
#'   probability weighting.
#' @param vce Variance estimator. Currently only `"bootstrap"` is supported.
#' @param reps Number of successful nonparametric bootstrap replications.
#'
#' @details
#' The reported QTE at index `tau` is the difference between quantiles of two
#' marginal potential-outcome distributions. It is not, without additional
#' assumptions, the quantile of unit-level treatment effects. Inference is
#' bootstrap-based because analytic standard errors require density estimation
#' at every estimated quantile.
#'
#' This implementation targets population quantiles for continuous outcomes.
#' It does not estimate quantile treatment effects on the treated.
#'
#' @return A `teffects_qte` object. The `potential_quantiles` element in
#'   `x$teffects` contains the estimated quantile for every treatment level.
#' @references
#' Cattaneo, M. D. (2010). Efficient semiparametric estimation of multi-valued
#' treatment effects under ignorability. *Journal of Econometrics*, 155(2),
#' 138--154.
#'
#' Cattaneo, M. D., Drukker, D. M., and Holland, A. D. (2013). Estimation of
#' multivalued treatment effects under conditional independence. *The Stata
#' Journal*, 13(3), 407--450.
#' @export
qte <- function(
  fml,
  data,
  probs = c(0.25, 0.5, 0.75),
  treatment_fml = NULL,
  method = c("eif", "ipw"),
  tmodel = c("logit", "probit"),
  weights = NULL,
  vce = "bootstrap",
  reps = 999L,
  pstolerance = 1e-5,
  ...
) {
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_arg(treatment_fml, "NULL os formula")
  dreamerr::check_set_arg(
    vce,
    "strict match(bootstrap)",
    .message = "Only `vce = \"bootstrap\"` is implemented for `qte()`."
  )
  dreamerr::check_arg(probs, "class(numeric, integer)")
  dreamerr::check_arg(
    probs,
    "numeric vector no na len(1,) GT{0} LT{1}"
  )
  dreamerr::check_value(reps, "integer scalar GE{0}")
  if (reps < 2L) {
    stop("Argument `reps` needs to be 2 or larger.")
  }
  dreamerr::check_value(
    pstolerance,
    "numeric scalar GE{0} LT{0.5}",
    .arg_name = "pstolerance",
  )
  dreamerr::check_set_arg(method, "strict match")
  dreamerr::check_set_arg(tmodel, "strict match")
  probs <- sort(unique(as.numeric(probs)))

  .teffects_qte(
    fml = fml,
    data = data,
    probs = probs,
    treatment_fml = treatment_fml,
    weights = weights,
    method = method,
    tmodel = tmodel,
    reps = reps,
    pstolerance = pstolerance,
    call = match.call()
  )
}


## Estimator ----

.teffects_qte <- function(
  fml,
  data,
  probs,
  treatment_fml,
  weights,
  method,
  tmodel,
  reps,
  pstolerance,
  call
) {
  use_eif <- identical(method, "eif")
  prepared <- .teffects_prepare(
    fml = fml,
    data = data,
    treatment_fml = treatment_fml,
    outcome_design = if (use_eif) {
      "reduced"
    } else if (identical(method, "ipw")) {
      "none"
    },
    treatment_design = "full",
    weights = weights
  )
  y <- prepared$outcome
  X <- prepared$outcome_design
  d <- prepared$treatment
  treatment_candidate <- prepared$treatment_design
  keep <- prepared$complete
  n <- prepared$n

  treated_groups <- prepared$treated_levels
  groups <- c(prepared$control_level, treated_groups)
  if (length(groups) < 2L) {
    stop("`qte()` requires at least two treatment levels.", call. = FALSE)
  }
  control <- prepared$control_level
  if (use_eif && !("(Intercept)" %in% colnames(X))) {
    stop(
      "`method = \"eif\"` requires an intercept in the outcome formula.",
      call. = FALSE
    )
  }
  fit <- .teffects_qte_core(
    y = y,
    X = X,
    d = d,
    groups = groups,
    treated_groups = treated_groups,
    control = control,
    treatment_candidate = treatment_candidate,
    probs = probs,
    method = method,
    tmodel = tmodel,
    pstolerance = pstolerance
  )

  bootstrap <- matrix(NA_real_, reps, length(fit$estimate))
  colnames(bootstrap) <- names(fit$estimate)
  completed <- 0L
  attempts <- 0L
  maximum_attempts <- max(5L * reps, reps + 20L)
  while (completed < reps && attempts < maximum_attempts) {
    attempts <- attempts + 1L
    index <- sample.int(n, n, replace = TRUE)
    replicate_fit <- tryCatch(
      .teffects_qte_core(
        y = y[index],
        X = if (use_eif) {
          X[index, , drop = FALSE]
        } else if (identical(method, "ipw")) {
          NULL
        },
        d = d[index],
        groups = groups,
        treated_groups = treated_groups,
        control = control,
        treatment_candidate = treatment_candidate[index, , drop = FALSE],
        probs = probs,
        method = method,
        tmodel = tmodel,
        pstolerance = pstolerance,
        diagnostics = FALSE
      ),
      error = function(e) NULL
    )
    if (
      !is.null(replicate_fit) &&
        length(replicate_fit$estimate) == ncol(bootstrap) &&
        all(is.finite(replicate_fit$estimate))
    ) {
      completed <- completed + 1L
      bootstrap[completed, ] <- replicate_fit$estimate
    }
  }
  if (completed < reps) {
    stop(
      sprintf(
        "Only %d of %d requested bootstrap replications succeeded after %d attempts.",
        completed,
        reps,
        attempts
      ),
      call. = FALSE
    )
  }
  covariance <- stats::cov(bootstrap)
  dimnames(covariance) <- list(names(fit$estimate), names(fit$estimate))

  .teffects_result(
    fit$estimate,
    estimator = "qte",
    details = list(
      stat = "qte",
      method = method,
      probs = probs,
      treatment_model_name = tmodel,
      control = control,
      treated = treated_groups,
      target = names(fit$estimate),
      sample = which(keep),
      osample = .teffects_osample(prepared$n_original, keep),
      formula = fml,
      treatment_formula = treatment_fml,
      call = call,
      propensity_score = fit$probability,
      treatment_model = fit$treatment_model,
      potential_quantiles = fit$potential_quantiles,
      cdf_weights = fit$cdf_weights,
      bootstrap = bootstrap,
      bootstrap_reps = reps,
      bootstrap_attempts = attempts
    ),
    covariance = covariance
  )
}

.teffects_qte_core <- function(
  y,
  X,
  d,
  groups,
  treated_groups = NULL,
  control = NULL,
  treatment_candidate,
  probs,
  method,
  tmodel,
  pstolerance,
  diagnostics = TRUE
) {
  if (is.null(control)) {
    control <- groups[[1L]]
  }
  if (is.null(treated_groups)) {
    treated_groups <- groups[-1L]
  }
  if (!diagnostics && !all(groups %in% d)) {
    stop("A bootstrap sample omitted a treatment level.", call. = FALSE)
  }
  n <- length(y)
  K <- length(groups)
  treatment_design <- .teffects_independent_columns(treatment_candidate)
  if (diagnostics && K == 2L) {
    treatment_model <- .teffects_glm(
      treatment_design,
      as.numeric(d == groups[[2L]]),
      tmodel
    )
    probability <- cbind(
      1 - treatment_model$fitted.values,
      treatment_model$fitted.values
    )
    colnames(probability) <- as.character(groups)
  } else {
    treatment_fit <- .teffects_fit_generalized_treatment_model(
      design = treatment_design,
      treatment = d,
      groups = groups,
      tmodel = tmodel,
      prediction_only = !diagnostics
    )
    probability <- treatment_fit$probability
    if (diagnostics) {
      treatment_model <- treatment_fit$model
    }
  }
  .teffects_overlap_check(probability, pstolerance)

  potential_quantiles <- matrix(
    NA_real_,
    nrow = length(probs),
    ncol = K,
    dimnames = if (diagnostics) {
      list(format(probs, trim = TRUE), as.character(groups))
    }
  )
  cdf_weights <- if (diagnostics) {
    stats::setNames(vector("list", K), as.character(groups))
  }
  mean_X <- if (identical(method, "eif")) {
    Matrix::colMeans(X)
  } else if (identical(method, "ipw")) {
    NULL
  }
  for (j in seq_len(K)) {
    rows <- d == groups[[j]]
    inverse_probability <- 1 / probability[rows, j]

    if (identical(method, "ipw")) {
      arm_weights <- inverse_probability / sum(inverse_probability)
    } else if (method == "eif") {
      X_group <- X[rows, , drop = FALSE]
      bread <- Matrix::crossprod(X_group)
      c_term <- mean_X -
        as.numeric(Matrix::crossprod(X_group, inverse_probability)) / n
      adjustment <- as.numeric(
        X_group %*% .teffects_solve(bread, c_term)
      )
      arm_weights <- inverse_probability / n + adjustment
    }
    if (diagnostics) {
      cdf_weights[[j]] <- arm_weights
    }
    potential_quantiles[, j] <- .teffects_cdf_quantile(
      y[rows],
      arm_weights,
      probs
    )
  }

  estimate <- .teffects_qte_contrasts(
    potential_quantiles,
    probs,
    groups,
    treated_groups,
    control
  )
  if (!diagnostics) {
    return(list(estimate = estimate))
  }

  list(
    estimate = estimate,
    potential_quantiles = potential_quantiles,
    cdf_weights = cdf_weights,
    probability = probability,
    treatment_model = treatment_model
  )
}

.teffects_cdf_quantile <- function(y, weights, probs) {
  order_y <- order(y)
  y <- y[order_y]
  weights <- weights[order_y]
  runs <- rle(y)
  values <- runs$values
  cumulative <- cumsum(weights)[cumsum(runs$lengths)]
  tolerance <- sqrt(.Machine$double.eps)
  monotone <- all(weights >= -tolerance) &&
    all(diff(cumulative) >= -tolerance)
  if (monotone) {
    crossing <- findInterval(probs, cumulative, left.open = TRUE) + 1L
    return(values[pmin(crossing, length(values))])
  }

  vapply(
    probs,
    function(prob) {
      # The EIF series estimator can produce a nonmonotone finite-sample CDF.
      # Cattaneo's estimator minimizes the absolute sample moment in this case.
      values[[which.min(abs(cumulative - prob))]]
    },
    numeric(1)
  )
}

.teffects_qte_contrasts <- function(
  potential_quantiles,
  probs,
  groups,
  treated_groups,
  control
) {
  treated <- match(treated_groups, groups)
  control_position <- match(control, groups)
  pieces <- vector("list", length(probs))
  piece_names <- vector("list", length(probs))
  for (r in seq_along(probs)) {
    pieces[[r]] <- c(
      potential_quantiles[r, treated] -
        potential_quantiles[r, control_position],
      potential_quantiles[r, control_position]
    )
    tau <- format(probs[[r]], trim = TRUE)
    piece_names[[r]] <- c(
      sprintf("QTE[%s vs %s]@%s", treated_groups, control, tau),
      sprintf("POquantile[%s]@%s", control, tau)
    )
  }
  estimate <- unlist(pieces, use.names = FALSE)
  names(estimate) <- unlist(piece_names, use.names = FALSE)
  estimate
}


## Methods ----

#' @export
print.teffects_qte <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Quantile treatment-effects estimation\n")
  cat(
    "Estimator:",
    if (identical(x$teffects$method, "eif")) {
      "efficient influence function"
    } else if (identical(x$teffects$method, "ipw")) {
      "inverse-probability weighting"
    },
    "\n"
  )
  cat("Treatment model:", x$teffects$treatment_model_name, "\n")
  cat("Bootstrap replications:", x$teffects$bootstrap_reps, "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
