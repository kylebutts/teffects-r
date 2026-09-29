test_that("probit-weighted RA reproduces Stata's cattaneo2 ATE", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- ipwra(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    treatment_fml = ~ mmarried + mage + I(mage^2) + fbaby + medu,
    data = data,
    psmethod = "probit"
  )
  reference <- stata_reference("ipwra", "ate", c("effect", "control_pom"))
  expect_equal(unname(coef(est)), reference$estimate, tolerance = 5e-7)
  expect_equal(
    unname(sqrt(diag(vcov(est)))),
    reference$std_error,
    tolerance = 5e-7
  )
  expect_equal(est$vcov, crossprod(est$influence) / nrow(est$influence)^2)
})

test_that("weighted RA pomeans reproduce Stata and expose auxiliary models", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- ipwra(
    bweight ~ treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby,
    treatment_fml = ~ mmarried + mage + I(mage^2) + fbaby + medu,
    data = data,
    psmethod = "probit",
    stat = "pomeans"
  )
  reference <- stata_reference(
    "ipwra",
    "pomeans",
    c("control_pom", "treated_pom")
  )
  expect_equal(unname(coef(est)), reference$estimate, tolerance = 5e-7)
  expect_equal(
    unname(sqrt(diag(vcov(est)))),
    reference$std_error,
    tolerance = 5e-7
  )
})

test_that("weighted RA reproduces the supported examples printed in the manual", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  formula <- bweight ~
    treat(mbsmoke, ref = 0) + prenatal1 + mmarried + mage + fbaby
  treatment_formula <- ~ mmarried + mage + I(mage^2) + fbaby + medu

  ate <- ipwra(
    formula,
    treatment_fml = treatment_formula,
    data = data,
    psmethod = "probit"
  )
  pomeans <- ipwra(
    formula,
    treatment_fml = treatment_formula,
    data = data,
    psmethod = "probit",
    stat = "pomeans"
  )

  expect_stata_result(
    pomeans,
    "ipwra",
    "pomeans",
    c("control_pom", "treated_pom"),
    "manual_ipwra_2"
  )
})

test_that("weighted RA supports multivalued treatments", {
  set.seed(32768)
  n <- 1000
  data <- data.frame(x = rnorm(n))
  latent <- cbind(0, 0.25 + 0.45 * data$x, -0.25 - 0.35 * data$x)
  probability <- exp(latent) / rowSums(exp(latent))
  data$d <- apply(probability, 1L, function(p) sample(0:2, 1L, prob = p))
  data$y <- 1 + data$d + data$x + rnorm(n)
  formula <- y ~ treat(d, ref = 0) + x

  ate <- ipwra(formula, data, psmethod = "logit")
  expect_named(coef(ate), c("ATE[1 vs 0]", "ATE[2 vs 0]", "POmean[0]"))

  att <- ipwra(
    formula,
    data,
    stat = "att",
    psmethod = "logit"
  )
  expect_named(
    coef(att),
    c("ATT[1 vs 0]", "ATT[2 vs 0]", "POmean[0|1]", "POmean[0|2]")
  )
})
