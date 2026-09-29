## Propensity score tools ----
# Flag observations that violate the propensity-score overlap tolerance and,
# when `error = TRUE`, fail loudly. Callers that trim or clamp request the
# flags with `error = FALSE`.
.teffects_overlap_check <- function(propensity, pstolerance, error = TRUE) {
  failure <- if (pstolerance == 0) {
    !is.finite(propensity) | propensity <= 0 | propensity >= 1
  } else {
    !is.finite(propensity) |
      propensity < pstolerance |
      propensity > 1 - pstolerance
  }
  if (!is.null(dim(propensity))) {
    failure <- rowSums(failure) > 0L
  }
  number_failed <- sum(failure)
  if (error && number_failed) {
    stop(
      sprintf(
        "%d observation%s fail the propensity-score overlap tolerance.",
        number_failed,
        if (number_failed == 1L) "" else "s"
      ),
      call. = FALSE
    )
  }
  failure
}


# Project each generalized propensity-score vector onto the probability
# simplex with every component in [tolerance, 1 - tolerance].  The projection
# preserves the adding-up constraint, unlike componentwise clipping when the
# treatment has more than two levels.
.teffects_clamp_probability <- function(probability, tolerance) {
  if (tolerance == 0) {
    return(probability)
  }
  vector_input <- is.null(dim(probability))
  if (vector_input) {
    return(pmin(pmax(probability, tolerance), 1 - tolerance))
  }
  probability <- as.matrix(probability)
  K <- ncol(probability)
  if (K * tolerance > 1) {
    stop(
      sprintf(
        "`pstolerance` must not exceed 1 / %d when `psaction = \"clamp\"`.",
        K
      ),
      call. = FALSE
    )
  }
  lower <- tolerance
  upper <- 1 - tolerance

  if (K == 2L) {
    # Binary treatment, just clamp
    treated <- pmin(pmax(probability[, 2L], lower), upper)
    clamped <- cbind(1 - treated, treated)
  } else {
    # Multivalued treatment,
    # The projection preserves the adding-up constraint, unlike componentwise
    # clamping when the treatment has more than two levels.
    clamped <- probability
    outside <- rowSums(probability < lower | probability > upper) > 0L
    remaining <- 1 - K * lower
    if (remaining == 0) {
      clamped[,] <- lower
    } else {
      for (i in which(outside)) {
        shifted <- probability[i, ] - lower
        ordered <- sort(shifted, decreasing = TRUE)
        excess <- (cumsum(ordered) - remaining) / seq_len(K)
        active <- max(which(ordered > excess))
        clamped[i, ] <- pmax(shifted - excess[[active]], 0) + lower
      }
    }
  }
  dimnames(clamped) <- dimnames(probability)
  clamped
}

.teffects_clamp_dlog_probability <- function(
  probability,
  clamped_probability,
  dlog_probability,
  tolerance
) {
  K <- ncol(probability)
  result <- lapply(dlog_probability, as.matrix)
  boundary_tolerance <- sqrt(.Machine$double.eps)
  boundary <- clamped_probability <= tolerance + boundary_tolerance |
    clamped_probability >= 1 - tolerance - boundary_tolerance
  for (i in which(rowSums(boundary) > 0L)) {
    free <- clamped_probability[i, ] > tolerance + boundary_tolerance &
      clamped_probability[i, ] < 1 - tolerance - boundary_tolerance
    if (!any(free)) {
      for (j in seq_len(K)) {
        result[[j]][i, ] <- 0
      }
      next
    }
    dp <- do.call(
      rbind,
      lapply(seq_len(K), function(j) {
        probability[i, j] * dlog_probability[[j]][i, ]
      })
    )
    dq <- matrix(0, K, ncol(dp))
    dq[free, ] <- sweep(
      dp[free, , drop = FALSE],
      2L,
      colMeans(dp[free, , drop = FALSE]),
      "-"
    )
    for (j in seq_len(K)) {
      result[[j]][i, ] <- dq[j, ] / clamped_probability[i, j]
    }
  }
  result
}
