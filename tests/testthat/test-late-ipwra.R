test_that("late(estimator = 'ipwra') estimates LATE with logical first-stage predictions", {
  set.seed(8192)
  n <- 1200
  data <- data.frame(x = rnorm(n))
  data$z <- rbinom(n, 1, plogis(-0.1 + 0.5 * data$x))
  latent <- -0.4 + 1.4 * data$z + 0.5 * data$x + rlogis(n)
  data$d <- as.numeric(latent > 0)
  data$y <- 1 + 2 * data$d + 0.8 * data$x + rnorm(n)

  fit <- teffects:::late(
    y ~ treat(d, ref = 0) + x,
    data,
    instrument = z,
    estimator = "ipwra"
  )
  expect_named(coef(fit), c("LATE", "ATE[Y|Z]", "ATE[D|Z]"))
  expect_equal(
    unname(coef(fit)[[1L]]),
    unname(coef(fit)[[2L]] / coef(fit)[[3L]]),
    tolerance = 1e-10
  )
  fitted_treatment <- unlist(lapply(
    fit$teffects$treatment_models,
    function(model) model$fitted.values
  ))
  expect_true(all(fitted_treatment > 0 & fitted_treatment < 1))
  expect_equal(fit$vcov, crossprod(fit$influence) / nrow(fit$influence)^2)
  expect_lt(abs(unname(coef(fit)[[1L]]) - 2), 0.35)
})

test_that("late(estimator = 'ipwra') estimates LATT and supports a separate instrument model", {
  set.seed(16384)
  n <- 1400
  data <- data.frame(x = rnorm(n), q = rnorm(n))
  data$z <- rbinom(n, 1, plogis(0.2 + 0.6 * data$x - 0.3 * data$q))
  latent <- -0.6 + 1.6 * data$z + 0.4 * data$x + rlogis(n)
  data$d <- as.numeric(latent > 0)
  data$y <- 0.5 + 1.5 * data$d + data$x + rnorm(n)

  fit <- teffects:::late(
    y ~ treat(d, ref = 0) + x,
    data,
    instrument = "z",
    instrument_fml = ~ x + q,
    estimator = "ipwra",
    stat = "latt",
    imodel = "probit"
  )
  expect_named(coef(fit), c("LATT", "ATT[Y|Z]", "ATT[D|Z]"))
  expect_equal(
    unname(coef(fit)[[1L]]),
    unname(coef(fit)[[2L]] / coef(fit)[[3L]]),
    tolerance = 1e-10
  )
  expect_equal(fit$teffects$instrument, "z")
  expect_lt(abs(unname(coef(fit)[[1L]]) - 1.5), 0.4)
})

test_that("late(estimator = 'ipwra') validates its binary treatment and instrument", {
  data <- data.frame(
    y = 1:6,
    d = c(0, 1, 2, 0, 1, 2),
    z = c(0, 1, 0, 1, 0, 1),
    x = c(-1, 0, 1, -1, 0, 1)
  )
  expect_error(
    teffects:::late(y ~ treat(d) + x, data, z, estimator = "ipwra"),
    "binary actual treatment"
  )
  data$d <- c(0, 1, 0, 1, 0, 1)
  data$z <- c(0, 1, 2, 0, 1, 2)
  expect_error(
    teffects:::late(y ~ treat(d) + x, data, z, estimator = "ipwra"),
    "exactly two observed levels"
  )
})

test_that("late(estimator = 'ipwra') supports one-sided noncompliance", {
  set.seed(32768)
  n <- 1000
  data <- data.frame(x = rnorm(n))
  data$z <- rbinom(n, 1, plogis(0.2 + 0.4 * data$x))
  data$d <- data$z * rbinom(n, 1, plogis(0.5 + 0.3 * data$x))
  data$y <- 1 + 1.25 * data$d + data$x + rnorm(n)

  for (target in c("late", "latt")) {
    fit <- teffects:::late(
      y ~ treat(d, ref = 0) + x,
      data,
      instrument = z,
      stat = target,
      estimator = "ipwra"
    )
    expect_equal(fit$teffects$treatment_models[[1L]]$degenerate_value, 0)
    expect_length(coef(fit$teffects$treatment_models[[1L]]), 0L)
    expect_lt(abs(unname(coef(fit)[[1L]]) - 1.25), 0.35)
  }
})
