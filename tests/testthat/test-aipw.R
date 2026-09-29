test_that("probit AIPW reproduces Stata's cattaneo2 ATE", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- aipw(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    data = data,
    treatment_fml = ~ mmarried + mage + I(mage^2) + fbaby + medu,
    psmethod = "probit"
  )
  reference <- stata_reference("aipw", "ate", c("effect", "control_pom"))
  expect_equal(unname(coef(est)), reference$estimate, tolerance = 5e-7)
  expect_equal(
    unname(sqrt(diag(vcov(est)))),
    reference$std_error,
    tolerance = 5e-7
  )
  expect_equal(est$vcov, crossprod(est$influence) / nrow(est$influence)^2)
})

test_that("probit AIPW reproduces Stata's cattaneo2 POMs", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- aipw(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    data = data,
    treatment_fml = ~ mmarried + mage + I(mage^2) + fbaby + medu,
    psmethod = "probit",
    stat = "pomeans"
  )
  reference <- stata_reference(
    "aipw",
    "pomeans",
    c("control_pom", "treated_pom")
  )
  expect_equal(unname(coef(est)), reference$estimate, tolerance = 5e-7)
  expect_equal(
    unname(sqrt(diag(vcov(est)))),
    reference$std_error,
    tolerance = 5e-7
  )
  expect_named(coef(est), c("POmean[0]", "POmean[1]"))
})

test_that("probit AIPW reproduces Stata's cattaneo2 ATT", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- aipw(
    bweight ~ treat(mbsmoke, ref = 0) + fbaby + mage + mmarried + prenatal1,
    data = data,
    treatment_fml = ~ fbaby + foreign + medu + mmarried,
    psmethod = "probit",
    stat = "att"
  )
  reference <- stata_reference("aipw", "att", c("effect", "control_pom"))
  expect_equal(unname(coef(est)), reference$estimate, tolerance = 5e-7)
  expect_equal(
    unname(sqrt(diag(vcov(est)))),
    reference$std_error,
    tolerance = 5e-7
  )
})

test_that("AIPW reproduces the supported examples printed in the manual", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  formula <- bweight ~
    treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby
  treatment_formula <- ~ mmarried + mage + I(mage^2) + fbaby + medu

  ate <- aipw(
    formula,
    data = data,
    treatment_fml = treatment_formula,
    psmethod = "probit"
  )
  pomeans <- aipw(
    formula,
    data = data,
    treatment_fml = treatment_formula,
    psmethod = "probit",
    stat = "pomeans"
  )
  att <- aipw(
    bweight ~ treat(mbsmoke, ref = 0) +
      fbaby +
      mage +
      mmarried +
      prenatal1,
    data = data,
    treatment_fml = ~ fbaby + foreign + medu + mmarried,
    psmethod = "probit",
    stat = "att"
  )

  expect_stata_result(
    pomeans,
    "aipw",
    "pomeans",
    c("control_pom", "treated_pom"),
    "manual_aipw_2"
  )
  expect_stata_result(
    att,
    "aipw",
    "att",
    c("effect", "control_pom"),
    "manual_aipw_5"
  )
})
