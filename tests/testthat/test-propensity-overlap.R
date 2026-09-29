test_that("multivalued clamping preserves the probability simplex", {
  probability <- rbind(
    c(0.999, 0.0005, 0.0005),
    c(0.8, 0.15, 0.05)
  )
  bounded <- teffects:::.teffects_clamp_probability(probability, 0.1)

  expect_equal(rowSums(bounded), c(1, 1), tolerance = 1e-12)
  expect_gte(min(bounded), 0.1 - 1e-12)
  expect_lte(max(bounded), 0.9 + 1e-12)
  expect_equal(
    teffects:::.teffects_clamp_probability(bounded, 0.1),
    bounded,
    tolerance = 1e-12
  )
  expect_error(
    teffects:::.teffects_clamp_probability(probability, 0.4),
    "must not exceed 1 / 3"
  )
})

test_that("binary clamping is equivalent to scalar clipping", {
  probability <- c(0.001, 0.2, 0.8, 0.999)
  expected <- pmin(pmax(probability, 0.1), 0.9)

  expect_equal(
    teffects:::.teffects_clamp_probability(probability, 0.1),
    expected
  )
  matrix_probability <- cbind(1 - probability, probability)
  bounded <- teffects:::.teffects_clamp_probability(matrix_probability, 0.1)
  expect_equal(bounded[, 2L], expected)
  expect_equal(rowSums(bounded), rep(1, length(probability)))
})
