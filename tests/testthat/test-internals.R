test_that("preparation keeps intercept-only designs compact and reusable", {
  data <- data.frame(
    y = c(4, 1, 6, 2),
    d = c(1, 0, 1, 0)
  )

  prepared <- teffects:::.teffects_prepare(y ~ treat(d), data)

  expect_equal(dim(prepared$outcome_design), c(nrow(data), 1L))
  expect_equal(
    as.matrix(prepared$outcome_design),
    matrix(
      1,
      nrow = nrow(data),
      dimnames = list(NULL, "(Intercept)")
    )
  )
  expect_identical(prepared$treatment_design, prepared$outcome_design)
  expect_identical(prepared$complete, rep(TRUE, nrow(data)))
  expect_identical(prepared$treated_levels, 1)
  expect_identical(c(prepared$control_level, prepared$treated_levels), c(0, 1))
})

test_that("preparation excludes non-finite rows without expanding designs", {
  data <- data.frame(
    y = seq_len(4),
    d = c(0, 1, 0, 1),
    x = c(1, NA, Inf, 4)
  )

  prepared <- teffects:::.teffects_prepare(y ~ treat(d) + x, data)

  expect_equal(dim(prepared$outcome_design), c(2L, 2L))
  expect_identical(prepared$complete, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(prepared$n, 2L)
})

test_that("preparation puts control first and sorts treated levels", {
  data <- data.frame(
    y = seq_len(6),
    d = c(3, 0, 2, 0, 1, 3)
  )

  prepared <- teffects:::.teffects_prepare(y ~ treat(d, ref = 0), data)

  expect_identical(prepared$treated_levels, c(1, 2, 3))
  expect_identical(
    c(prepared$control_level, prepared$treated_levels),
    c(0, 1, 2, 3)
  )
})

test_that("complete-case masks combine vectors, matrices, and frames", {
  sparse_matrix <- Matrix::Matrix(
    matrix(c(1, NA, 3), ncol = 1L),
    sparse = TRUE
  )
  extra <- data.frame(group = factor(c("a", NA, "b")))

  expect_identical(
    teffects:::.teffects_complete_cases(
      c(1, Inf, 3),
      matrix(c(1, 2, 3, 4, 5, 6), nrow = 3L),
      sparse_matrix,
      extra
    ),
    c(TRUE, FALSE, TRUE)
  )
  expect_error(
    teffects:::.teffects_complete_cases(1:3, matrix(1:4, nrow = 2L)),
    "same number of rows"
  )
})

test_that("prepared designs preserve formula semantics and factor levels", {
  data <- data.frame(
    y = 1:6,
    d = rep(0:1, 3L),
    x = 2:7,
    group = factor(c("a", NA, "b", "a", "b", "a"))
  )

  prepared <- teffects:::.teffects_prepare(
    y ~ treat(d) + I(x^2) + log(x) + factor(group),
    data
  )

  expect_identical(prepared$complete, c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))
  expected <- Matrix::sparse.model.matrix(
    ~ I(x^2) + log(x) + factor(group),
    data[prepared$complete, , drop = FALSE]
  )
  rownames(expected) <- NULL
  expect_equal(
    as.matrix(prepared$outcome_design),
    as.matrix(expected)
  )
})

test_that("prepared designs preserve no-intercept formulas and namespaced treat", {
  data <- data.frame(y = 1:4, d = c(0, 1, 0, 1), x = 2:5)

  prepared <- teffects:::.teffects_prepare(
    y ~ 0 + teffects::treat(d, ref = 0) + I(x^2),
    data
  )

  expect_identical(colnames(prepared$outcome_design), "I(x^2)")
  expect_equal(as.numeric(prepared$outcome_design[, 1L]), data$x^2)
})

test_that("preparation aligns the common sample with its complete-case mask", {
  data <- data.frame(
    y = seq_len(5),
    d = c(0, 1, 0, 1, 0),
    x = c(1, 2, NA, 4, 5),
    z = c(1, NA, 3, 4, 5)
  )
  rownames(data) <- paste0("row-", seq_len(nrow(data)))
  prepared <- teffects:::.teffects_prepare(y ~ treat(d) + x, data, ~z)

  expect_identical(prepared$complete, c(TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(prepared$outcome, data$y[prepared$complete])
  expect_equal(prepared$treatment, data$d[prepared$complete])
  expect_equal(nrow(prepared$outcome_design), sum(prepared$complete))
  expect_null(rownames(prepared$outcome_design))
  expect_identical(
    teffects:::.teffects_osample(prepared$n_original, prepared$complete),
    c(FALSE, TRUE, TRUE, FALSE, FALSE)
  )
})

test_that("rank reduction occurs after the common sample is known", {
  data <- data.frame(
    y = 1:5,
    d = c(0, 1, 0, 1, 0),
    group = factor(c("a", "a", "a", "b", "b")),
    z = c(1, 2, 3, NA, NA)
  )
  reduced <- teffects:::.teffects_prepare(y ~ treat(d) + group, data, ~z)
  candidates <- teffects:::.teffects_prepare(
    y ~ treat(d) + group,
    data,
    ~z,
    outcome_design = "full",
    treatment_design = "full"
  )

  expect_identical(reduced$complete, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(Matrix::rankMatrix(reduced$outcome_design)[[1L]], 1L)
  expect_identical(colnames(reduced$outcome_design), "(Intercept)")
  expect_gt(ncol(candidates$outcome_design), 1L)
})

test_that("unreduced designs retain the full-data factor universe", {
  data <- data.frame(
    y = 1:6,
    d = rep(0:1, 3L),
    g = c("a", "b", "c", "a", "b", "c"),
    z = c(1, 2, NA, 4, 5, NA)
  )
  prepared <- teffects:::.teffects_prepare(
    y ~ treat(d) + factor(g),
    data,
    ~z,
    outcome_design = "full"
  )

  expect_true("factor(g)c" %in% colnames(prepared$outcome_design))
  expect_true(all(prepared$outcome_design[, "factor(g)c"] == 0))
})

test_that("control level are selected from the final common sample", {
  data <- data.frame(
    y = 1:4,
    d = c(0, 0, 1, 1),
    x = 1:4,
    z = c(NA, NA, 1, 2)
  )

  implicit <- teffects:::.teffects_prepare(y ~ treat(d) + x, data, ~z)
  expect_error(
    ipwra(y ~ treat(d, ref = 0) + x, data, treatment_fml = ~z, psmethod = "logit"),
    "reference level must be observed in the analysis sample."
  )

  expect_identical(implicit$control_level, 1)
})
