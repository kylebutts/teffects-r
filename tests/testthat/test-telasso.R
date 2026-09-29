test_that("telasso estimates ATEs with plugin-selected nuisance models", {
  set.seed(20260918)
  n <- 800
  x <- matrix(rnorm(n * 12L), n)
  data <- data.frame(x)
  names(data) <- paste0("x", seq_len(ncol(x)))
  data$d <- rbinom(n, 1, plogis(0.7 * data$x1 - 0.4 * data$x2))
  data$y <- 2 * data$d + data$x1 + 0.5 * data$x3 + rnorm(n)

  fit <- telasso(
    y ~ treat(d, ref = 0) + x..,
    data = data
  )

  expect_named(coef(fit), c("ATE[1 vs 0]", "POmean[0]"))
  expect_equal(unname(coef(fit)[[1L]]), 2, tolerance = 0.2)
  expect_equal(fit$vcov, crossprod(fit$influence) / n^2)
  expect_true(all(fit$teffects$propensity_score > 0))
  expect_true(all(fit$teffects$propensity_score < 1))
})

test_that("telasso reproduces Stata's cattaneo2 reference result", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  fit <- telasso(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    data = data,
    treatment_fml = ~ mmarried + mage + I(mage^2) + fbaby + medu
  )

  expect_stata_result(fit, "telasso", "ate", c("effect", "control_pom"))
})

test_that("telasso cross-fitting is reproducible and preserves RNG state", {
  set.seed(99)
  n <- 500
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  data$d <- rbinom(n, 1, plogis(data$x1))
  data$y <- data$d + data$x1 + rnorm(n)
  formula <- y ~ treat(d) + x1 + x2 + x3
  state <- .Random.seed

  first <- telasso(formula, data, xfolds = 3, resamples = 2, seed = 42)
  expect_identical(.Random.seed, state)
  second <- telasso(formula, data, xfolds = 3, resamples = 2, seed = 42)

  expect_equal(coef(first), coef(second))
  expect_equal(vcov(first), vcov(second))
  expect_length(first$teffects$nuisance, 2L)
})

test_that("telasso validates its currently supported scope", {
  data <- data.frame(
    y = rnorm(30),
    d = rep(0:2, 10),
    x = rnorm(30)
  )
  expect_error(
    telasso(y ~ treat(d) + x, data),
    "binary treatment"
  )
  data$d <- rep(0:1, 15)
  expect_error(
    telasso(y ~ treat(d) + x, data, omodel = "poisson"),
    "Only `omodel"
  )
  expect_error(
    telasso(y ~ treat(d) + x, data, resamples = 2),
    "requires `xfolds > 1`"
  )
})
