test_that("CBPS matches the CRAN reference implementation", {
  skip_if_not_installed("CBPS")

  set.seed(20260918)
  n <- 1000
  binary <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  binary$d <- rbinom(
    n,
    1,
    stats::plogis(0.2 + 0.5 * binary$x1 - 0.3 * binary$x2)
  )
  binary$y <- 1 + 2 * binary$d + binary$x1 + rnorm(n)

  for (stat in c("ate", "att")) {
    ours <- ipw(
      y ~ treat(d, ref = 0) + x1 + x2,
      binary,
      stat = stat,
      psmethod = "cbps"
    )
    reference <- CBPS::CBPS(
      d ~ x1 + x2,
      binary,
      ATT = if (stat == "ate") 0 else 1,
      method = "exact",
      standardize = TRUE
    )
    ours_probability <- ours$teffects$propensity_score[, 2L]
    reference_probability <- stats::fitted(reference)
    ours_weight <- ifelse(
      binary$d == 1,
      1 / ours_probability,
      if (stat == "ate") {
        1 / (1 - ours_probability)
      } else {
        ours_probability / (1 - ours_probability)
      }
    )
    ours_weight <- ave(
      ours_weight,
      binary$d,
      FUN = function(x) x / sum(x)
    )

    expect_lt(
      max(abs(
        as.numeric(stats::coef(ours$teffects$treatment_model)) -
          as.numeric(stats::coef(reference))
      )),
      3e-4
    )
    expect_equal(ours_probability, reference_probability, tolerance = 2e-4)
    expect_equal(
      as.numeric(ours_weight),
      as.numeric(reference$weights),
      tolerance = 5e-3
    )
    expect_lt(
      max(abs(
        sqrt(diag(stats::vcov(ours$teffects$treatment_model))) /
          sqrt(diag(reference$var)) -
          1
      )),
      0.25
    )
  }
})

test_that("multivalued CBPS matches the CRAN reference implementation", {
  skip_if_not_installed("CBPS")

  set.seed(20260919)
  n <- 1200
  multiple <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  linear_predictor <- cbind(
    0,
    0.2 + 0.4 * multiple$x1,
    -0.3 + 0.3 * multiple$x2
  )
  probability <- exp(linear_predictor) / rowSums(exp(linear_predictor))
  multiple$d <- apply(
    probability,
    1L,
    function(p) sample(0:2, 1L, prob = p)
  )
  multiple$y <- 1 + multiple$d + multiple$x1 + rnorm(n)

  ours <- ipw(
    y ~ treat(d, ref = 0) + x1 + x2,
    multiple,
    psmethod = "cbps"
  )
  reference <- CBPS::CBPS(
    factor(d) ~ x1 + x2,
    multiple,
    ATT = 0,
    method = "exact",
    standardize = TRUE
  )

  expect_equal(
    ours$teffects$propensity_score,
    stats::fitted(reference),
    tolerance = 5e-5
  )
})
