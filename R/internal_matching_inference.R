## Match predictions ----

.teffects_predict_matches <- function(
  y,
  treatment,
  matching,
  bias_predictions = NULL
) {
  po0 <- po1 <- rep(NA_real_, length(y))
  po0[treatment == 0L] <- y[treatment == 0L]
  po1[treatment == 1L] <- y[treatment == 1L]
  focal_treated <- matching$focal[treatment[matching$focal] == 1L]
  focal_control <- matching$focal[treatment[matching$focal] == 0L]
  plain_donor_outcome <- if (
    (length(focal_treated) && is.null(bias_predictions$control)) ||
      (length(focal_control) && is.null(bias_predictions$treated))
  ) {
    as.numeric(matching$operator %*% y)
  } else {
    NULL
  }

  if (length(focal_treated)) {
    donor_outcome <- if (is.null(bias_predictions$control)) {
      plain_donor_outcome
    } else {
      as.numeric(
        matching$operator %*% (y - bias_predictions$control)
      ) +
        bias_predictions$control
    }
    po0[focal_treated] <- donor_outcome[focal_treated]
  }
  if (length(focal_control)) {
    donor_outcome <- if (is.null(bias_predictions$treated)) {
      plain_donor_outcome
    } else {
      as.numeric(
        matching$operator %*% (y - bias_predictions$treated)
      ) +
        bias_predictions$treated
    }
    po1[focal_control] <- donor_outcome[focal_control]
  }
  list(po0 = po0, po1 = po1, effect = po1 - po0)
}


## Local moments ----

.teffects_local_moments <- function(plan, y, z = NULL) {
  operator <- plan$operator
  degree <- plan$degree
  focal <- plan$focal
  mean_y <- as.numeric(operator %*% y)
  mean_y2 <- as.numeric(operator %*% (y^2))
  mean_y[focal[degree == 0L]] <- NA_real_
  mean_y2[focal[degree == 0L]] <- NA_real_
  correction <- degree / (degree - 1)
  variance <- rep(NA_real_, length(y))
  valid <- degree > 1L
  variance[focal[valid]] <- pmax(
    correction[valid] *
      (mean_y2[focal[valid]] - mean_y[focal[valid]]^2),
    0
  )

  result <- list(mean = mean_y, variance = variance)
  if (is.null(z)) {
    return(result)
  }

  mean_z <- as.matrix(operator %*% z)
  mean_zy <- as.matrix(
    operator %*% .teffects_row_multiply(z, y)
  )
  covariance <- matrix(0, nrow(z), ncol(z))
  covariance[focal[valid], ] <-
    correction[valid] *
    (mean_zy[focal[valid], , drop = FALSE] -
      mean_z[focal[valid], , drop = FALSE] * mean_y[focal[valid]])
  colnames(covariance) <- colnames(z)
  result$covariance <- covariance
  result
}


## Matching variance ----

.teffects_matching_variance <- function(
  y,
  treatment,
  matching,
  predictions,
  estimate,
  stat,
  vce,
  vce_nneighbor,
  index,
  same_neighbor_offset = 1L,
  same_plan = NULL
) {
  n <- length(y)
  target <- matching$target
  reuse <- matching$reuse

  if (identical(vce, "robust")) {
    if (is.null(same_plan)) {
      variance_rows <- if (stat %in% c("att", "atc")) {
        which(target | reuse$k > 0)
      } else {
        seq_len(n)
      }
      same_plan <- .teffects_build_neighborhood_plan(
        index,
        vce_nneighbor + same_neighbor_offset,
        same = TRUE,
        include_self = TRUE,
        rows = variance_rows
      )
    }
    conditional_variance <- .teffects_local_moments(
      same_plan,
      y
    )$variance
    if (anyNA(conditional_variance[target])) {
      stop(
        paste(
          "Too few within-treatment matches to estimate the robust",
          "matching variance."
        ),
        call. = FALSE
      )
    }
  } else if (identical(vce, "iid")) {
    operator <- matching$operator
    edge_i <- operator@i + 1L
    edge_j <- rep.int(seq_len(ncol(operator)), diff(operator@p))
    pair_effect <-
      (2 * treatment[edge_i] - 1) *
      (y[edge_i] - y[edge_j])
    squared_residual <- sum(operator@x * (pair_effect - estimate)^2) /
      length(matching$focal)
    conditional_variance <- rep(squared_residual / 2, n)
  }

  if (stat %in% c("att", "atc")) {
    contribution <- numeric(n)
    contribution[target] <- (predictions$effect[target] - estimate)^2
    donor <- !target & reuse$k > 0
    contribution[donor] <- conditional_variance[donor] *
      (reuse$k[donor]^2 - reuse$k_prime[donor])
    variance <- sum(contribution) / sum(target)^2
  } else if (identical(stat, "ate")) {
    contribution <-
      (predictions$effect - estimate)^2 +
      conditional_variance *
        (reuse$k^2 + 2 * reuse$k - reuse$k_prime)
    variance <- sum(contribution) / n^2
  }

  list(
    estimate = estimate,
    variance = variance,
    conditional_variance = conditional_variance,
    reuse = reuse,
    contribution = contribution
  )
}


## Propensity-score adjustment ----

.teffects_psmatch_adjustment <- function(
  y,
  treatment,
  propensity,
  treatment_model,
  design,
  stat,
  estimate,
  vce_nneighbor,
  nneighbor,
  index,
  same_plan,
  clamped = rep(FALSE, length(propensity))
) {
  if (identical(stat, "atc")) {
    treatment <- 1L - treatment
    propensity <- 1 - propensity
    stat <- "att"
    estimate <- -estimate
  }
  z <- design
  gamma_vcov <- stats::vcov(treatment_model)
  density <- treatment_model$family$mu.eta(treatment_model$linear.predictors)
  density[clamped] <- 0

  opposite_plan <- .teffects_build_neighborhood_plan(
    index,
    vce_nneighbor,
    same = FALSE,
    include_self = FALSE
  )
  covariance_same <- .teffects_local_moments(
    same_plan,
    y,
    z
  )$covariance
  opposite_moments <- .teffects_local_moments(
    opposite_plan,
    y,
    z
  )
  covariance_opposite <- opposite_moments$covariance
  covariance_1 <- covariance_0 <- matrix(0, nrow(z), ncol(z))
  covariance_1[treatment == 1L, ] <- covariance_same[treatment == 1L, ]
  covariance_1[treatment == 0L, ] <- covariance_opposite[treatment == 0L, ]
  covariance_0[treatment == 0L, ] <- covariance_same[treatment == 0L, ]
  covariance_0[treatment == 1L, ] <- covariance_opposite[treatment == 1L, ]

  if (identical(stat, "ate")) {
    c_tau <- colMeans(
      density * (covariance_1 / propensity + covariance_0 / (1 - propensity))
    )
    adjustment <- -drop(t(c_tau) %*% gamma_vcov %*% c_tau)
    return(list(adjustment = adjustment, c = c_tau))
  }

  same_minus_i <- .teffects_build_neighborhood_plan(
    index,
    vce_nneighbor,
    same = TRUE,
    include_self = FALSE
  )
  y_same <- as.numeric(same_minus_i$operator %*% y)
  y_opposite <- opposite_moments$mean
  y1_tilde <- ifelse(treatment == 1L, y_same, y_opposite)
  y0_tilde <- ifelse(treatment == 0L, y_same, y_opposite)
  n_treated <- sum(treatment)
  c1 <- as.numeric(Matrix::crossprod(
    z,
    density * (y1_tilde - y0_tilde - estimate)
  )) /
    n_treated
  c2 <- as.numeric(Matrix::crossprod(
    density,
    covariance_1 + propensity / (1 - propensity) * covariance_0
  )) /
    n_treated
  c_delta <- c1 + c2

  z_no_intercept <- z[, colnames(z) != "(Intercept)", drop = FALSE]
  dense_z <- as.matrix(z_no_intercept)
  z_scaling <- tryCatch(
    solve(stats::cov(dense_z)),
    error = function(error) {
      stop(
        "The propensity-score covariates have a singular covariance matrix.",
        call. = FALSE
      )
    }
  )
  covariate_distance <- .teffects_make_quadratic_distance(
    dense_z,
    z_scaling
  )
  covariate_index <- .teffects_matching_index(
    covariate_distance,
    treatment
  )
  omega_z <- .teffects_build_neighborhood_plan(
    covariate_index,
    nneighbor,
    same = FALSE,
    include_self = FALSE
  )
  y_omega_z <- as.numeric(omega_z$operator %*% y)
  derivative <- as.numeric(Matrix::crossprod(
    z,
    density * ((2 * treatment - 1) * (y - y_omega_z) - estimate)
  )) /
    n_treated
  adjustment <-
    -drop(t(c_delta) %*% gamma_vcov %*% c_delta) +
    drop(t(derivative) %*% gamma_vcov %*% derivative)
  list(adjustment = adjustment, c = c_delta, derivative = derivative)
}


## Result ----

.teffects_matching_result <- function(
  y,
  treatment,
  matching,
  stat,
  level_names,
  vce,
  vce_nneighbor,
  index,
  ps_adjustment = NULL,
  same_neighbor_offset = 1L,
  predictions,
  estimate,
  same_plan = NULL,
  metadata
) {
  inference <- .teffects_matching_variance(
    y = y,
    treatment = treatment,
    matching = matching,
    predictions = predictions,
    estimate = estimate,
    stat = stat,
    vce = vce,
    vce_nneighbor = vce_nneighbor,
    index = index,
    same_neighbor_offset = same_neighbor_offset,
    same_plan = same_plan
  )
  target_variance <- inference$variance +
    if (is.null(ps_adjustment)) 0 else ps_adjustment$adjustment
  if (!is.finite(target_variance) || target_variance < 0) {
    stop("The matching variance estimate is not positive.", call. = FALSE)
  }
  coefficient_name <- sprintf(
    "%s[%s vs %s]",
    toupper(stat),
    level_names[[2L]],
    level_names[[1L]]
  )
  covariance <- matrix(
    target_variance,
    dimnames = list(coefficient_name, coefficient_name)
  )
  .teffects_result(
    stats::setNames(estimate, coefficient_name),
    estimator = metadata$estimator,
    details = c(
      metadata[setdiff(names(metadata), "estimator")],
      list(
        stat = stat,
        target = coefficient_name,
        matching_weights = matching$weights,
        matching_degree = matching$degree,
        matching_edge_count = matching$edge_count,
        matching_inference = inference,
        propensity_adjustment = ps_adjustment,
        potential_outcomes = predictions
      )
    ),
    covariance = covariance
  )
}
