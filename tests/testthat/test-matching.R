test_that("scalar matching distances reproduce dense neighbor searches", {
  values <- c(0.02, 0.07, 0.07, 0.16, 0.24, 0.28, 0.42, 0.52)
  treatment <- rep(0:1, 4L)
  caliper <- 0.2
  scalar <- structure(
    list(values = values, caliper = caliper),
    class = "teffects_scalar_distance"
  )
  dense <- abs(outer(values, values, "-"))
  dense[dense > caliper] <- Inf

  for (same in c(FALSE, TRUE)) {
    for (include_self in c(FALSE, TRUE)) {
      expected <- teffects:::.teffects_neighbor_sets(
        dense,
        treatment,
        number = 2L,
        same = same,
        include_self = include_self
      )
      actual <- teffects:::.teffects_neighbor_sets(
        scalar,
        treatment,
        number = 2L,
        same = same,
        include_self = include_self
      )
      expect_identical(actual, expected)
    }
  }
})

test_that("scalar neighbor searches preserve multilevel ties and calipers", {
  values <- c(-0.20, -0.10, -0.10, 0, 0, 0.10, 0.10, 0.30, 0.30)
  treatment <- rep(c("control", "low", "high"), 3L)
  rows <- c(9L, 2L, 5L, 1L, 7L)

  for (caliper in list(NULL, 0.11)) {
    scalar <- structure(
      list(values = values, caliper = caliper),
      class = "teffects_scalar_distance"
    )
    dense <- abs(outer(values, values, "-"))
    if (!is.null(caliper)) {
      dense[dense > caliper] <- Inf
    }

    for (same in c(FALSE, TRUE)) {
      for (include_self in c(FALSE, TRUE)) {
        expected <- teffects:::.teffects_neighbor_sets(
          dense,
          treatment,
          number = 2L,
          same = same,
          include_self = include_self,
          rows = rows,
          tolerance = 0
        )
        actual <- teffects:::.teffects_neighbor_sets(
          scalar,
          treatment,
          number = 2L,
          same = same,
          include_self = include_self,
          rows = rows,
          tolerance = 0
        )

        expect_identical(actual, expected)
      }
    }
  }
})
