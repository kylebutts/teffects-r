test_that("probit IPW reproduces Stata's cattaneo2 ATE", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))

  est <- ipw(
    bweight ~ treat(mbsmoke, ref = 0) +
      mmarried +
      mage +
      I(mage^2) +
      fbaby +
      medu,
    data = data,
    psmethod = "probit"
  )
  reference <- stata_reference("ipw", "ate", c("effect", "control_pom"))

  expect_equal(
    unname(coef(est)),
    reference$estimate,
    tolerance = 5e-7
  )
  expect_equal(
    unname(sqrt(diag(vcov(est)))),
    reference$std_error,
    tolerance = 5e-7
  )
  expect_equal(est$vcov, crossprod(est$influence) / nrow(est$influence)^2)
})

test_that("probit IPW reproduces Stata's cattaneo2 ATT", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))

  est <- ipw(
    bweight ~ treat(mbsmoke, ref = 0) +
      mmarried +
      mage +
      I(mage^2) +
      fbaby +
      medu,
    data = data,
    psmethod = "probit",
    stat = "att"
  )
  reference <- stata_reference("ipw", "att", c("effect", "control_pom"))

  expect_equal(
    unname(coef(est)),
    reference$estimate,
    tolerance = 5e-7
  )
  expect_equal(
    unname(sqrt(diag(vcov(est)))),
    reference$std_error,
    tolerance = 5e-7
  )
})

test_that("IPW reproduces the examples printed in the manual", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  formula <- bweight ~ treat(mbsmoke, ref = 0) +
    mmarried +
    mage +
    I(mage^2) +
    fbaby +
    medu

  ate <- ipw(formula, data = data, psmethod = "probit")
  att <- ipw(
    formula,
    data = data,
    psmethod = "probit",
    stat = "att"
  )

  for (case in "manual_ipw_3") {
    expect_stata_result(
      ate,
      "ipw",
      "ate",
      c("effect", "control_pom"),
      case
    )
  }
})

test_that("multivalued propensity scores form a probability simplex", {
  set.seed(8192)
  n <- 1200
  data <- data.frame(x = rnorm(n))
  latent <- cbind(0, 0.2 + 0.4 * data$x, -0.3 - 0.5 * data$x)
  probability <- exp(latent) / rowSums(exp(latent))
  data$d <- apply(probability, 1L, function(p) sample(0:2, 1L, prob = p))
  data$y <- 1 + data$d + data$x + rnorm(n)

  fit <- ipw(y ~ treat(d, ref = 0) + x, data)
  expect_equal(dim(fit$teffects$propensity_score), c(n, 3L))
  expect_equal(
    rowSums(fit$teffects$propensity_score),
    rep(1, n),
    tolerance = 1e-10
  )
})
