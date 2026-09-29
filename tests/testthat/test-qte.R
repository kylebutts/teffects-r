test_that("QTE recovers shifted empirical quantiles without covariates", {
  base <- stats::qnorm((1:200 - 0.5) / 200)
  data <- data.frame(
    d = rep(0:1, each = 200),
    y = c(base, base + 2)
  )

  set.seed(1048576)
  eif <- qte(
    y ~ treat(d, ref = 0),
    data,
    probs = c(0.25, 0.5, 0.75),
    method = "eif",
    reps = 20
  )
  set.seed(1048576)
  ipw <- qte(
    y ~ treat(d, ref = 0),
    data,
    probs = c(0.25, 0.5, 0.75),
    method = "ipw",
    reps = 20
  )

  expect_equal(unname(coef(eif)[c(1, 3, 5)]), rep(2, 3), tolerance = 1e-12)
  expect_equal(coef(eif), coef(ipw), tolerance = 1e-12)
  expect_equal(dim(vcov(eif)), c(6L, 6L))
  expect_null(eif$influence)
  expect_equal(dim(eif$teffects$potential_quantiles), c(3L, 2L))
})

test_that("QTE validates its quantiles and EIF specification", {
  data <- data.frame(y = 1:20, d = rep(0:1, each = 10), x = rep(1:10, 2))

  expect_error(
    qte(y ~ treat(d), data, probs = c(0, 0.5), reps = 2),
    "`probs`"
  )
  expect_error(
    qte(y ~ treat(d), data, probs = 0.5, vce = "robust", reps = 2),
    "Only `vce"
  )
  expect_error(
    qte(y ~ 0 + treat(d) + x, data, probs = 0.5, reps = 2),
    "requires an intercept"
  )
})
