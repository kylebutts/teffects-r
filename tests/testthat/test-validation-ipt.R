test_that("IPT reproduces Graham's iptATE calculation", {
  data <- read_stata(stata_fixture("RPSPolicyEvalData.dta"))
  y <- data$log_calories_pc2002
  d <- data$RPS
  X <- cbind(
    log_calories_pc2000 = data$log_calories_pc2000,
    log_real_exp_pc2000 = data$log_real_exp_pc2000,
    HHSize2000 = data$HHSize2000,
    `(Intercept)` = 1
  )
  n <- nrow(X)
  threshold <- log(1 / (n - 1))
  a <- -(n - 1) * (1 + threshold + 0.5 * threshold^2)
  b <- n + (n - 1) * threshold
  c <- -(n - 1)

  phi <- function(v) {
    ifelse(v > threshold, v - exp(-v), a + b * v + 0.5 * c * v^2)
  }
  phi1 <- function(v) {
    ifelse(v > threshold, 1 + exp(-v), b + c * v)
  }
  phi2 <- function(v) {
    ifelse(v > threshold, -exp(-v), c)
  }
  fit_tilt <- function(selected) {
    objective <- function(beta) {
      eta <- drop(X %*% beta)
      -sum(selected * phi(eta) - eta)
    }
    gradient <- function(beta) {
      eta <- drop(X %*% beta)
      -drop(crossprod(X, selected * phi1(eta) - 1))
    }
    stats::optim(
      stats::lm.fit(X, selected)$coefficients,
      objective,
      gradient,
      method = "BFGS",
      control = list(reltol = 1e-13, maxit = 10000)
    )$par
  }

  treated_beta <- fit_tilt(d)
  control_beta <- fit_tilt(1 - d)
  treated_eta <- drop(X %*% treated_beta)
  control_fit_eta <- drop(X %*% control_beta)
  control_ado_eta <- -control_fit_eta
  treated_probability <- stats::plogis(treated_eta)
  control_probability <- stats::plogis(control_fit_eta)
  reference_ate <- mean(
    d / treated_probability * y - (1 - d) / control_probability * y
  )

  treated_score <- (d * phi1(treated_eta) - 1) * X
  control_score <- ((1 - d) * phi1(control_fit_eta) - 1) * X
  effect_score <- d / treated_probability * y -
    (1 - d) / control_probability * (y + reference_ate)
  treated_hessian <- crossprod(X, d * phi2(treated_eta) * X)
  control_hessian <- crossprod(X, (1 - d) * phi2(control_fit_eta) * X)
  effect_treated_derivative <- colSums(d * -exp(-treated_eta) * y * X)
  effect_control_derivative <- colSums(
    (1 - d) * exp(control_ado_eta) * (y + reference_ate) * X
  )
  p <- ncol(X)
  derivative <- rbind(
    cbind(treated_hessian, matrix(0, p, p), 0),
    cbind(matrix(0, p, p), control_hessian, 0),
    c(effect_treated_derivative, effect_control_derivative, -n)
  )
  score <- cbind(treated_score, control_score, effect_score)
  reference_vcov <- solve(derivative) %*%
    crossprod(score) %*%
    t(solve(derivative))

  fit <- ipw(
    log_calories_pc2002 ~
      treat(RPS, ref = 0) +
      log_calories_pc2000 +
      log_real_exp_pc2000 +
      HHSize2000,
    data,
    pstolerance = 0,
    psmethod = "ipt"
  )

  expect_lt(
    abs(unname(stats::coef(fit)[[1L]]) - reference_ate),
    1e-6
  )
  expect_lt(
    abs(
      sqrt(stats::vcov(fit)[1L, 1L]) -
        sqrt(reference_vcov[2L * p + 1L, 2L * p + 1L])
    ),
    1e-8
  )
})
