test_that("nnmatch uses internal nearest-neighbor weights for ATE", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- nnmatch(
    bweight ~ treat(mbsmoke, ref = 0) +
      mage +
      prenatal1 +
      mmarried +
      fbaby,
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

test_that("nnmatch reproduces the remaining examples printed in the manual", {
  data <- read_stata(stata_fixture("cattaneo3.dta"))

  exact_continuous <- nnmatch(
    bweight ~ treat(mbsmoke, ref = 0) + mage,
    data = data,
    ematch = ~ prenatal1 + mmarried + fbaby,
    metric = "euclidean"
  )
  exact_only <- nnmatch(
    bweight ~ treat(mbsmoke, ref = 0),
    data = data,
    ematch = ~ factor(mmarried) + factor(magecat)
  )

  expect_stata_result(
    exact_continuous,
    "nnmatch",
    "ate",
    "effect",
    "manual_nnmatch_2"
  )
  expect_stata_result(
    exact_only,
    "nnmatch",
    "ate",
    "effect",
    "manual_nnmatch_4"
  )
})

test_that("nnmatch applies exact matching and multiple neighbors", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  est <- nnmatch(
    bweight ~ treat(mbsmoke, ref = 0) + mage,
    data = data,
    stat = "att",
    nneighbor = 2,
    ematch = ~ prenatal1 + mmarried + fbaby,
    metric = "euclidean",
    generate = "neighbor"
  )

  expect_named(coef(est), "ATT[1 vs 0]")
  expect_true(all(est$teffects$matching_degree >= 2L))
  expect_named(est$teffects$generated_matches, "att")
  generated_names <- names(as.data.frame(
    est$teffects$generated_matches$att
  ))
  expect_gte(length(generated_names), 2L)
  expect_identical(
    generated_names,
    paste0("neighbor", seq_along(generated_names))
  )
})

test_that("nnmatch bias adjustment uses sparse outcome models", {
  data <- read_stata(stata_fixture("cattaneo2.dta"))
  formula <- bweight ~ treat(mbsmoke, ref = 0) + mage + fage
  unadjusted <- nnmatch(
    formula,
    data = data,
    ematch = ~ prenatal1 + mmarried + fbaby
  )
  adjusted <- nnmatch(
    formula,
    data = data,
    ematch = ~ prenatal1 + mmarried + fbaby,
    biasadj = ~ mage + fage
  )

  expect_equal(
    unname(coef(adjusted) - coef(unadjusted)),
    adjusted$teffects$bias_adjustment
  )
  expect_stata_result(
    adjusted,
    "nnmatch",
    "ate",
    "effect",
    "manual_nnmatch_3"
  )
})

test_that("nnmatch aligns generated matches with the original sample", {
  data <- data.frame(
    y = c(2, 4, 3, 7, 5, 9, 6, 11, 8, 13, 10, 15),
    d = rep(0:1, 6L),
    x = seq_len(12),
    exact = rep(c("a", "a", "b", "b"), 3L)
  )
  data$exact[[1L]] <- NA_character_

  fit <- nnmatch(
    y ~ treat(d) + x,
    data,
    ematch = ~exact,
    biasadj = ~1,
    metric = "euclidean",
    vce = "iid",
    generate = "neighbor"
  )

  expect_true(fit$teffects$osample[[1L]])
  for (matches in fit$teffects$generated_matches) {
    focal <- as.integer(rownames(matches))
    for (i in seq_len(nrow(matches))) {
      matched <- as.integer(stats::na.omit(matches[i, ]))
      expect_true(all(data$d[matched] != data$d[focal[[i]]]))
      expect_true(all(data$exact[matched] == data$exact[focal[[i]]]))
    }
  }
})

test_that("nnmatch validates distance controls", {
  data <- data.frame(y = 1:6, d = rep(0:1, 3L), x = 1:6)

  expect_error(nnmatch(y ~ treat(d) + x, data, dtolerance = NA_real_))
  expect_error(nnmatch(y ~ treat(d) + x, data, caliper = Inf))
  expect_error(nnmatch(y ~ treat(d) + x, data, nneighbor = 0L))
  expect_error(nnmatch(y ~ treat(d) + x, data, vce_nneighbor = 0L))
  expect_error(
    nnmatch(y ~ treat(d) + x, data, weights = rep(1, nrow(data))),
    "not currently supported"
  )
  expect_error(
    nnmatch(
      y ~ treat(d) + x,
      data,
      metric = "matrix",
      metric_matrix = matrix(-1)
    ),
    "positive semidefinite"
  )
})

test_that("nnmatch preserves expanded covariate names", {
  data <- data.frame(
    y = seq_len(8),
    d = rep(0:1, 4L),
    x1 = 1:8,
    x2 = 8:1
  )
  fixest::setFixest_fml(..x = ~ x1 + x2)

  prepared <- teffects:::.teffects_prepare(
    fixest::xpd(y ~ treat(d) + ..x),
    data,
    outcome_design = "full",
    treatment_design = "none"
  )

  expect_identical(
    colnames(prepared$outcome_design),
    c("(Intercept)", "x1", "x2")
  )
})

test_that("nnmatch follows fixest naming for logical covariates", {
  data <- data.frame(
    y = seq_len(8),
    d = rep(0:1, 4L),
    unemployed = rep(c(FALSE, TRUE), 4L),
    nodegree = rep(c(TRUE, FALSE), 4L)
  )

  prepared <- teffects:::.teffects_prepare(
    y ~ treat(d) + unemployed + nodegree,
    data,
    outcome_design = "full",
    treatment_design = "none"
  )

  expect_identical(
    colnames(prepared$outcome_design),
    c("(Intercept)", "unemployedTRUE", "nodegreeTRUE")
  )
  expect_equal(
    as.vector(prepared$outcome_design),
    as.vector(stats::model.matrix(~ unemployed + nodegree, data))
  )
})

test_that("nnmatch handles singular Mahalanobis covariance matrices", {
  data <- data.frame(
    y = seq_len(12),
    d = rep(0:1, 6L),
    group = factor(rep(letters[1:3], 4L)),
    x = 1:12
  )

  expect_silent(
    fit <- nnmatch(
      y ~ treat(d) + fixest::i(group) + x,
      data = data
    )
  )
  expect_true(is.finite(unname(coef(fit))))
})

test_that("k-d tree neighbor search preserves exhaustive matching", {
  skip_if_not_installed("FNN")
  set.seed(771)
  n <- 120
  design <- matrix(rnorm(n * 4), n, 4)
  treatment <- rep(0:1, length.out = n)
  exact <- rep(letters[1:3], length.out = n)
  distance <- structure(
    list(
      design = design,
      scaling = crossprod(matrix(rnorm(16), 4, 4)),
      caliper = 2
    ),
    class = "teffects_quadratic_distance"
  )

  for (same in c(TRUE, FALSE)) {
    for (include_self in c(TRUE, FALSE)) {
      expect_identical(
        teffects:::.teffects_neighbor_sets(
          distance,
          treatment,
          number = 3,
          same = same,
          include_self = include_self,
          exact_group = exact,
          use_kd = TRUE
        ),
        teffects:::.teffects_neighbor_sets(
          distance,
          treatment,
          number = 3,
          same = same,
          include_self = include_self,
          exact_group = exact,
          use_kd = FALSE
        )
      )
    }
  }
})
