test_that("IPT balances every treatment arm against the full sample", {
  set.seed(1024)
  n <- 600
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  latent <- cbind(0, 0.4 + data$x1, -0.3 + data$x2)
  probability <- exp(latent) / rowSums(exp(latent))
  data$d <- apply(probability, 1L, function(p) sample(0:2, 1L, prob = p))
  data$y <- 1 + data$d + data$x1 - data$x2 + rnorm(n)

  fit <- ipw(
    y ~ treat(d, ref = 0) + x1 + x2,
    data,
    psmethod = "ipt"
  )
  expect_named(coef(fit), c("ATE[1 vs 0]", "ATE[2 vs 0]", "POmean[0]"))

  X <- cbind(1, data$x1, data$x2)
  p <- fit$teffects$propensity_score
  for (j in seq_len(3L)) {
    weighted <- colMeans((data$d == j - 1L) / p[, j] * X)
    expect_equal(weighted, colMeans(X), tolerance = 1e-7)
  }
})

test_that("binary CBPS satisfies its ATE and ATT balance equations", {
  set.seed(2048)
  n <- 800
  data <- data.frame(x = rnorm(n))
  data$d <- rbinom(n, 1, plogis(0.2 + 0.7 * data$x))
  data$y <- 2 * data$d + data$x + rnorm(n)
  X <- cbind(1, data$x)

  ate <- ipw(
    y ~ treat(d, ref = 0) + x,
    data,
    stat = "ate",
    psmethod = "cbps"
  )
  p_ate <- ate$teffects$propensity_score[, 2L]
  expect_equal(
    colMeans((data$d - p_ate) / (p_ate * (1 - p_ate)) * X),
    c(0, 0),
    tolerance = 1e-7
  )

  att <- ipw(
    y ~ treat(d, ref = 0) + x,
    data,
    stat = "att",
    psmethod = "cbps"
  )
  p_att <- att$teffects$propensity_score[, 2L]
  expect_equal(
    colMeans((data$d - p_att) / (1 - p_att) * X),
    c(0, 0),
    tolerance = 1e-7
  )
  expect_named(coef(att), c("ATT[1 vs 0]", "POmean[0]"))
})

test_that("multivalued CBPS balances adjacent treatment arms", {
  set.seed(4096)
  n <- 900
  data <- data.frame(x = rnorm(n))
  latent <- cbind(0, 0.3 + 0.5 * data$x, -0.2 - 0.4 * data$x)
  probability <- exp(latent) / rowSums(exp(latent))
  data$d <- apply(probability, 1L, function(p) sample(0:2, 1L, prob = p))
  data$y <- data$d + data$x + rnorm(n)

  fit <- ipw(y ~ treat(d, ref = 0) + x, data, psmethod = "cbps")
  p <- fit$teffects$propensity_score
  X <- cbind(1, data$x)
  for (j in 2:3) {
    balance <- colMeans(
      (data$d == j - 1L) / p[, j] * X - (data$d == j - 2L) / p[, j - 1L] * X
    )
    expect_equal(balance, c(0, 0), tolerance = 1e-7)
  }
  expect_named(coef(fit), c("ATE[1 vs 0]", "ATE[2 vs 0]", "POmean[0]"))
})

# Equivalences from: Słoczyński, Uysal, and Wooldridge
test_that("IPT ATT uses normalized odds weights and equals linear IPWRA", {
  set.seed(1536)
  n <- 700
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  data$d <- rbinom(n, 1, plogis(0.1 + 0.6 * data$x1 - 0.4 * data$x2))
  data$y <- 1 + 2 * data$d + data$x1 + data$d * data$x2 + rnorm(n)

  fit <- ipw(
    y ~ treat(d, ref = 0) + x1 + x2,
    data,
    stat = "att",
    psmethod = "ipt"
  )
  p0 <- fit$teffects$propensity_score[, 1L]
  odds <- (1 - p0) / p0
  X <- cbind(1, data$x1, data$x2)
  treated_mean <- colMeans(X[data$d == 1, , drop = FALSE])
  weighted_control_mean <- colSums(
    (data$d == 0) * odds * X
  ) /
    sum((data$d == 0) * odds)
  expect_equal(weighted_control_mean, treated_mean, tolerance = 1e-7)
  expect_equal(sum((data$d == 0) * odds), sum(data$d), tolerance = 1e-7)

  control_regression <- lm(y ~ x1 + x2, data, weights = (data$d == 0) * odds)
  ipwra_control_mean <- mean(predict(control_regression, data[data$d == 1, ]))
  ipw_effect <- mean(data$y[data$d == 1]) -
    weighted.mean(data$y[data$d == 0], odds[data$d == 0])
  ipwra_effect <- mean(data$y[data$d == 1]) - ipwra_control_mean
  expect_equal(ipw_effect, ipwra_effect, tolerance = 1e-8)
  expect_equal(unname(coef(fit)[[1L]]), ipw_effect, tolerance = 1e-8)
  expect_named(coef(fit), c("ATT[1 vs 0]", "POmean[0]"))
})


test_that("linear adjustment is numerically identical for IPT ATE", {
  set.seed(8192)
  n <- 700
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  data$d <- rbinom(n, 1, plogis(0.2 + 0.6 * data$x1 - 0.4 * data$x2))
  data$y <- 1 + 2 * data$d + data$x1 + data$d * data$x2 + rnorm(n)
  fml <- y ~ treat(d, ref = 0) + x1 + x2

  weighting <- ipw(fml, data, stat = "ate", psmethod = "ipt")
  augmented <- aipw(fml, data, stat = "ate", psmethod = "ipt")
  adjusted <- ipwra(fml, data, stat = "ate", psmethod = "ipt")

  expect_equal(coef(augmented), coef(weighting), tolerance = 1e-10)
  expect_equal(coef(adjusted), coef(weighting), tolerance = 1e-10)
  expect_equal(vcov(augmented), vcov(weighting), tolerance = 1e-10)
  expect_equal(vcov(adjusted), vcov(weighting), tolerance = 1e-10)
})

test_that("linear adjustment is numerically identical for CBPS ATT", {
  set.seed(16384)
  n <- 700
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  data$d <- rbinom(n, 1, plogis(-0.1 + 0.5 * data$x1 + 0.3 * data$x2))
  data$y <- 1 + 2 * data$d + data$x1 - data$d * data$x2 + rnorm(n)
  fml <- y ~ treat(d, ref = 0) + x1 + x2

  weighting <- ipw(fml, data, stat = "att", psmethod = "cbps")
  augmented <- aipw(fml, data, stat = "att", psmethod = "cbps")
  adjusted <- ipwra(fml, data, stat = "att", psmethod = "cbps")

  expect_equal(coef(augmented), coef(weighting), tolerance = 1e-10)
  expect_equal(coef(adjusted), coef(weighting), tolerance = 1e-10)
  expect_equal(vcov(augmented), vcov(weighting), tolerance = 1e-10)
  expect_equal(vcov(adjusted), vcov(weighting), tolerance = 1e-10)
})
