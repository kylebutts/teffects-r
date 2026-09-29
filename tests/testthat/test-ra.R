test_that("linear RA reproduces Stata's cattaneo2 ATE", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))

  est <- ra(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    data = data
  )
  reference <- stata_reference("ra", "ate", c("effect", "control_pom"))

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

test_that("linear RA reproduces Stata's cattaneo2 ATT", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))

  est <- ra(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    data = data,
    stat = "att"
  )
  reference <- stata_reference("ra", "att", c("effect", "control_pom"))

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

test_that("linear RA reproduces Stata's cattaneo2 POMs", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))

  est <- ra(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    data = data,
    stat = "pomeans"
  )
  reference <- stata_reference(
    "ra",
    "pomeans",
    c("control_pom", "treated_pom")
  )

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
  expect_length(est$teffects$auxiliary$coefficients, 10L)
  expect_equal(
    unname(est$teffects$auxiliary$coefficients),
    unname(unlist(
      lapply(
        est$teffects$outcome_models,
        stats::coef
      ),
      use.names = FALSE
    ))
  )
})

test_that("treat(ref=) selects the RA reference level", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))

  est <- ra(
    bweight ~ treat(mbsmoke, ref = 1) + prenatal1 + mmarried + mage + fbaby,
    data = data
  )
  reference <- stata_reference(
    "ra",
    "pomeans",
    c("treated_pom", "control_pom")
  )

  expect_named(coef(est), c("ATE[0 vs 1]", "POmean[1]"))
  expect_equal(
    unname(coef(est)),
    c(
      reference$estimate[[2L]] - reference$estimate[[1L]],
      reference$estimate[[1L]]
    ),
    tolerance = 5e-7
  )
})

test_that("RA reproduces the supported examples printed in the manual", {
  cattaneo2 <- read_stata(stata_fixture("cattaneo2.dta"))
  formula <- bweight ~
    treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby

  ate <- ra(formula, data = cattaneo2)
  att <- ra(formula, data = cattaneo2, stat = "att")
  pomeans <- ra(formula, data = cattaneo2, stat = "pomeans")

  for (case in "manual_ra_4") {
    expect_stata_result(
      ate,
      "ra",
      "ate",
      c("effect", "control_pom"),
      case
    )
  }
  expect_stata_result(
    pomeans,
    "ra",
    "pomeans",
    c("control_pom", "treated_pom"),
    "manual_ra_3"
  )

  cattaneo3 <- read_stata(stata_fixture("cattaneo3.dta"))
  comparison <- ra(
    bweight ~ treat(mbsmoke, ref = 0) +
      factor(mmarried) * factor(magecat),
    data = cattaneo3
  )
  expect_stata_result(
    comparison,
    "ra",
    "ate",
    c("effect", "control_pom"),
    "manual_nnmatch_4_ra_compare"
  )
})

test_that("ra performs ordinary regression adjustment", {
  set.seed(77)
  data <- data.frame(x = rnorm(200))
  data$d <- rbinom(200, 1, plogis(data$x))
  data$y <- 1 + data$d + data$x + rnorm(200)

  default <- ra(y ~ treat(d, ref = 0) + x, data)
  expect_identical(default$teffects$estimator, "ra")
})

test_that("RA rejects outcome models that are not implemented yet", {
  expect_error(
    ra(mpg ~ treat(am, ref = 0) + wt, mtcars, omodel = "poisson"),
    'Only `omodel = "linear"`'
  )
  expect_error(
    ra(mpg ~ treat(am, ref = 0) + wt, mtcars, omodel = "logit"),
    'Only `omodel = "linear"`'
  )
})

test_that("RA defaults to zero and otherwise the first observed treatment", {
  numeric_data <- data.frame(
    y = c(4, 1, 6, 2),
    d = c(1, 0, 1, 0)
  )
  numeric_est <- ra(y ~ treat(d), numeric_data)
  expect_named(coef(numeric_est), c("ATE[1 vs 0]", "POmean[0]"))

  character_data <- data.frame(
    y = c(4, 1, 6, 2),
    d = c("treated", "control", "treated", "control")
  )
  character_est <- ra(y ~ treat(d), character_data)
  expect_named(
    coef(character_est),
    c("ATE[control vs treated]", "POmean[treated]")
  )
})
