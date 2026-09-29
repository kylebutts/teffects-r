## Output ----
# Assemble the fitted object
.teffects_result <- function(
  coefficients,
  influence = NULL,
  estimator,
  details,
  covariance = NULL
) {
  if (is.null(covariance)) {
    n <- nrow(influence)
    covariance <- as.matrix(Matrix::crossprod(influence)) / n^2
  }

  dimnames(covariance) <- list(names(coefficients), names(coefficients))
  standard_error <- sqrt(diag(covariance))
  statistic <- coefficients / standard_error
  table <- cbind(
    Estimate = coefficients,
    `Std. Error` = standard_error,
    `z value` = statistic,
    `Pr(>|z|)` = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  )

  structure(
    list(
      coefficients = coefficients,
      vcov = covariance,
      se = standard_error,
      coeftable = table,
      influence = influence,
      teffects = c(list(estimator = estimator), details)
    ),
    class = c(paste0("teffects_", estimator), "teffects")
  )
}

### Methods ----
# Report potential-outcome means, and for a treatment-effect statistic the
# contrast of the corresponding influence values.
.teffects_format_population <- function(
  pomean,
  influence,
  groups,
  control_level,
  stat
) {
  control_position <- match(control_level, groups)
  if (identical(stat, "pomeans")) {
    estimate <- pomean
    names(estimate) <- sprintf("POmean[%s]", groups)
  } else if (identical(stat, "ate")) {
    control_group <- groups[[control_position]]
    treated_groups <- groups[groups != control_group]
    treated <- match(treated_groups, groups)
    estimate <- c(
      pomean[treated] - pomean[[control_position]],
      pomean[[control_position]]
    )
    influence <- cbind(
      influence[, treated, drop = FALSE] - influence[, control_position],
      influence[, control_position]
    )
    names(estimate) <- c(
      sprintf("ATE[%s vs %s]", groups[treated], control_level),
      sprintf("POmean[%s]", control_level)
    )
  }
  colnames(influence) <- names(estimate)
  list(estimate = estimate, influence = influence)
}

# Report treatment-effect contrasts from means already stacked in pair order
# (treated level, control level).
.teffects_format_pairs <- function(
  pair_mean,
  pair_influence,
  groups,
  control_level,
  treated,
  stat
) {
  treated_mean <- seq.int(1L, length(pair_mean), by = 2L)
  control_mean <- treated_mean + 1L
  effect <- pair_mean[treated_mean] - pair_mean[control_mean]
  effect_influence <- pair_influence[, treated_mean, drop = FALSE] -
    pair_influence[, control_mean, drop = FALSE]

  if (identical(stat, "atc")) {
    estimate <- c(effect, pair_mean[[control_mean[[1L]]]])
    influence <- cbind(
      effect_influence,
      pair_influence[, control_mean[[1L]], drop = FALSE]
    )
    control_names <- sprintf("POmean[%s]", control_level)
  } else if (identical(stat, "att")) {
    estimate <- c(effect, pair_mean[control_mean])
    influence <- cbind(
      effect_influence,
      pair_influence[, control_mean, drop = FALSE]
    )
    control_names <- if (length(treated) == 1L) {
      sprintf("POmean[%s]", control_level)
    } else {
      sprintf("POmean[%s|%s]", control_level, groups[treated])
    }
  }
  names(estimate) <- c(
    sprintf(
      "%s[%s vs %s]",
      toupper(stat),
      groups[treated],
      control_level
    ),
    control_names
  )
  colnames(influence) <- names(estimate)
  list(estimate = estimate, influence = influence)
}

#' @export
coef.teffects <- function(object, ...) object$coefficients

#' @export
vcov.teffects <- function(object, ...) object$vcov

#' Confidence intervals for treatment-effect estimates
#'
#' Computes normal-approximation confidence intervals from the estimated
#' coefficients and covariance matrix of a treatment-effects model.
#'
#' @param object A fitted `teffects` object.
#' @param parm Optional numeric or character vector selecting coefficients.
#' @param level Confidence level as a proportion in `(0, 1)`.
#' @param ... Additional arguments passed to [stats::confint.default()].
#'
#' @return A matrix with one row per selected coefficient and columns giving
#'   the lower and upper confidence limits.
#' @export
confint.teffects <- function(object, parm, level = 0.95, ...) {
  message <- "`level` must be one number strictly between zero and one."
  dreamerr::check_value(
    level,
    "class(numeric, integer)",
    .arg_name = "level",
    .message = message
  )
  dreamerr::check_value(
    level,
    "numeric scalar GT{0} LT{1}",
    .arg_name = "level",
    .message = message
  )
  if (missing(parm)) {
    stats::confint.default(object, level = level, ...)
  } else {
    stats::confint.default(object, parm = parm, level = level, ...)
  }
}
