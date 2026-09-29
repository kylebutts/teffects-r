test_that("psmatch uses internal matching weights for ATE", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- psmatch(
    bweight ~ treat(mbsmoke, ref = 0) +
      mmarried +
      mage +
      I(mage^2) +
      fbaby +
      medu,
    data = data
  )

  keep <- est$teffects$sample
  treatment <- data$mbsmoke[keep] == 1
  outcome <- data$bweight[keep]
  matching_weight <- est$teffects$matching_weights
  expected <- stats::weighted.mean(
    outcome[treatment],
    matching_weight[treatment]
  ) -
    stats::weighted.mean(
      outcome[!treatment],
      matching_weight[!treatment]
    )
  expect_equal(unname(coef(est)), expected)
})

test_that("psmatch reproduces the remaining examples printed in the manual", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  formula <- bweight ~ treat(mbsmoke, ref = 0) +
    mmarried +
    mage +
    I(mage^2) +
    fbaby +
    medu

  caliper <- psmatch(formula, data = data, caliper = 0.1)
  att <- psmatch(
    formula,
    data = data,
    stat = "att",
    vce = "iid",
    caliper = 0.03
  )
  four_neighbors <- psmatch(formula, data = data, nneighbor = 4L)

  expect_stata_result(
    caliper,
    "psmatch",
    "ate",
    "effect",
    "manual_psmatch_2_caliper"
  )
  expect_stata_result(
    att,
    "psmatch",
    "att",
    "effect",
    "manual_psmatch_2_att"
  )
  expect_stata_result(
    four_neighbors,
    "psmatch",
    "ate",
    "effect",
    "manual_psmatch_3"
  )
  expect_error(
    psmatch(formula, data = data, caliper = 0.03),
    "could not be matched"
  )
})

test_that("psmatch reports generated matches in original row coordinates", {
  set.seed(902)
  data <- data.frame(x = stats::rnorm(80))
  data$d <- stats::rbinom(80, 1, stats::plogis(0.4 * data$x))
  data$y <- 1 + data$d + data$x + stats::rnorm(80)
  data$x[[1L]] <- NA_real_

  fit <- psmatch(
    y ~ treat(d) + x,
    data,
    stat = "att",
    vce = "iid",
    generate = "neighbor"
  )
  matches <- fit$teffects$generated_matches$att
  focal <- as.integer(rownames(matches))
  for (i in seq_len(nrow(matches))) {
    matched <- as.integer(stats::na.omit(matches[i, ]))
    expect_true(all(data$d[matched] != data$d[focal[[i]]]))
  }
  expect_true(fit$teffects$osample[[1L]])
})

test_that("psmatch accepts ATC", {
  set.seed(1902)
  data <- data.frame(x = stats::rnorm(120))
  data$d <- stats::rbinom(120, 1, stats::plogis(data$x))
  data$y <- 2 + 1.5 * data$d + data$x + stats::rnorm(120)

  fit <- psmatch(
    y ~ treat(d, ref = 0) + x,
    data,
    stat = "atc"
  )

  expect_s3_class(fit, "teffects_psmatch")
  expect_identical(fit$teffects$stat, "atc")
  expect_named(coef(fit), "ATC[1 vs 0]")
  expect_true(is.finite(unname(coef(fit))))
  expect_true(is.finite(unname(vcov(fit))))
})

test_that("psmatch rejects unsupported matching controls", {
  data <- data.frame(y = 1:8, d = rep(0:1, 4L), x = 1:8)

  expect_error(psmatch(y ~ treat(d) + x, data, nneighbor = 0L))
  expect_error(psmatch(y ~ treat(d) + x, data, vce_nneighbor = 0L))
  expect_error(
    psmatch(y ~ treat(d) + x, data, weights = rep(1, nrow(data))),
    "not currently supported"
  )
})

test_that("psmatch uses psmethod and rejects unimplemented methods", {
  data <- data.frame(y = 1:8, d = rep(0:1, 4L), x = 1:8)

  expect_silent(psmatch(y ~ treat(d) + x, data, psmethod = "probit"))
  expect_error(
    psmatch(y ~ treat(d) + x, data, psmethod = "ipt"),
    "not implemented yet"
  )
})
