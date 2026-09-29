atc_test_data <- function(n = 500L, seed = 90210L) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  propensity <- stats::plogis(-0.2 + 0.7 * x1 - 0.4 * x2)
  d <- stats::rbinom(n, 1, propensity)
  y <- 1 + (1.5 + 0.4 * x1) * d + x1 - 0.5 * x2 + stats::rnorm(n)
  data.frame(y, d, x1, x2)
}

test_that("ATC is the sign-reversed ATT after reversing treatment labels", {
  data <- atc_test_data(350L, 8128L)
  data$d_reverse <- 1 - data$d
  formula <- y ~ treat(d, ref = 0) + x1 + x2
  reverse_formula <- y ~ treat(d_reverse, ref = 0) + x1 + x2
  estimators <- list(
    ra = function(fml, stat) ra(fml, data, stat = stat),
    ipw = function(fml, stat) ipw(fml, data, stat = stat),
    aipw = function(fml, stat) aipw(fml, data, stat = stat),
    ipwra = function(fml, stat) {
      ipwra(fml, data, stat = stat, psmethod = "logit")
    },
    psmatch = function(fml, stat) {
      psmatch(fml, data, stat = stat, vce = "iid")
    },
    nnmatch = function(fml, stat) {
      nnmatch(fml, data, stat = stat, vce = "iid")
    }
  )

  for (estimator in estimators) {
    atc <- estimator(formula, "atc")
    reversed_att <- estimator(reverse_formula, "att")
    expect_equal(unname(coef(atc)[[1L]]), -unname(coef(reversed_att)[[1L]]))
    expect_equal(unname(vcov(atc)[1L, 1L]), unname(vcov(reversed_att)[1L, 1L]))
  }
})
