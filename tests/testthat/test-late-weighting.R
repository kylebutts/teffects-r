test_that("late(estimator = 'kappa') computes the normalized weighting estimator", {
  set.seed(65536)
  n <- 1200
  data <- data.frame(x = rnorm(n))
  data$z <- rbinom(n, 1, plogis(-0.2 + 0.6 * data$x))
  data$d <- rbinom(n, 1, plogis(-0.5 + 1.5 * data$z + 0.3 * data$x))
  data$y <- 1 + 1.75 * data$d + 0.8 * data$x + rnorm(n)

  fit <- teffects:::late(
    y ~ treat(d, ref = 0),
    data,
    instrument = z,
    instrument_fml = ~x,
    estimator = "kappa"
  )
  expect_named(coef(fit), c("LATE", "ATE[Y|Z]", "ATE[D|Z]"))
  expect_equal(
    unname(coef(fit)[[1L]]),
    unname(coef(fit)[[2L]] / coef(fit)[[3L]]),
    tolerance = 1e-10
  )

  e <- fit$teffects$propensity_score
  w0 <- (1 - data$z) / (1 - e)
  w1 <- data$z / e
  reduced_form <- weighted.mean(data$y, w1) - weighted.mean(data$y, w0)
  first_stage <- weighted.mean(data$d, w1) - weighted.mean(data$d, w0)
  expect_equal(
    unname(coef(fit)),
    c(
      reduced_form / first_stage,
      reduced_form,
      first_stage
    ),
    tolerance = 1e-9
  )
  expect_equal(fit$vcov, crossprod(fit$influence) / nrow(fit$influence)^2)
})

test_that("late(estimator = 'kappa') reproduces the Stata manual example", {
  skip_if_not_installed("haven")
  card <- read_stata(
    stata_fixture("card95.dta")
  )
  fit <- teffects:::late(
    lwage ~ treat(somecol),
    card,
    instrument = nearc4
  )
  expect_equal(unname(coef(fit)[["LATE"]]), 1.278672, tolerance = 5e-7)
  expect_equal(
    unname(sqrt(vcov(fit)["LATE", "LATE"])),
    0.2203624,
    tolerance = 5e-7
  )
})

test_that("late(estimator = 'balancing') exactly balances instrument covariates", {
  set.seed(131072)
  n <- 1000
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  data$z <- rbinom(n, 1, plogis(0.1 + 0.7 * data$x1 - 0.4 * data$x2))
  data$d <- rbinom(n, 1, plogis(-0.6 + 1.7 * data$z + 0.4 * data$x1))
  data$y <- 0.5 + 2 * data$d + data$x1 - 0.5 * data$x2 + rnorm(n)

  fit <- teffects:::late(
    y ~ treat(d, ref = 0),
    data,
    instrument = "z",
    instrument_fml = ~ x1 + x2,
    estimator = "balancing"
  )
  X <- cbind(1, data$x1, data$x2)
  e <- fit$teffects$propensity_score
  expect_equal(
    colMeans((data$z / e - (1 - data$z) / (1 - e)) * X),
    c(0, 0, 0),
    tolerance = 1e-7
  )
  expect_equal(unname(colSums(fit$teffects$normalized_weights)), c(1, 1))
})

test_that("late weighting estimators are translation and scale invariant", {
  set.seed(262144)
  n <- 900
  data <- data.frame(x = rnorm(n))
  data$z <- rbinom(n, 1, plogis(0.2 + 0.5 * data$x))
  data$d <- rbinom(n, 1, plogis(-0.8 + 1.8 * data$z + 0.2 * data$x))
  data$y <- 2 + 1.5 * data$d + data$x + rnorm(n)

  for (estimator in c("kappa", "balancing")) {
    original <- teffects:::late(
      y ~ treat(d),
      data,
      z,
      instrument_fml = ~x,
      estimator = estimator
    )
    shifted <- transform(data, y = 10 + 3 * y)
    transformed <- teffects:::late(
      y ~ treat(d),
      shifted,
      z,
      instrument_fml = ~x,
      estimator = estimator
    )
    expect_equal(
      unname(coef(transformed)[[1L]]),
      3 * unname(coef(original)[[1L]]),
      tolerance = 1e-8
    )
  }
})

test_that("late weighting estimators validate binary inputs", {
  data <- data.frame(
    y = 1:6,
    d = c(0, 1, 2, 0, 1, 2),
    z = c(0, 1, 0, 1, 0, 1),
    x = c(-1, 0, 1, -1, 0, 1)
  )
  expect_error(
    teffects:::late(
      y ~ treat(d),
      data,
      z,
      instrument_fml = ~x,
      estimator = "kappa"
    ),
    "binary actual treatment"
  )
  data$d <- c(0, 1, 0, 1, 0, 1)
  data$z <- c(0, 1, 2, 0, 1, 2)
  expect_error(
    teffects:::late(
      y ~ treat(d),
      data,
      z,
      instrument_fml = ~x,
      estimator = "balancing"
    ),
    "exactly two observed levels"
  )
})

test_that("late estimators exclude non-finite numeric instruments", {
  expect_identical(
    teffects:::.teffects_late_complete(c(0, 1, Inf, -Inf, NA_real_)),
    c(TRUE, TRUE, FALSE, FALSE, FALSE)
  )
})
