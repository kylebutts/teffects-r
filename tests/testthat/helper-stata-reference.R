.stata_reference_data <- local({
  reference <- NULL
  function() {
    if (is.null(reference)) {
      reference <<- utils::read.csv(
        system.file(
          "stata/teffects_reference.csv",
          package = "teffects",
          mustWork = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
    reference
  }
})

# Stata value labels are metadata, not data. Strip them on load so the
# matrix builders see plain numeric/factor columns.
read_stata <- function(path) {
  haven::zap_labels(haven::read_dta(path))
}

stata_fixture <- function(file) {
  testthat::test_path("fixtures", file)
}

stata_reference <- function(estimator, statistic, terms, case = "default") {
  reference <- .stata_reference_data()
  rows <- reference$estimator == estimator &
    reference$statistic == statistic &
    reference$case == case
  selected <- reference[rows, , drop = FALSE]
  positions <- match(terms, selected$term)
  if (anyNA(positions)) {
    stop("Requested Stata reference terms are unavailable.", call. = FALSE)
  }
  selected[positions, c("estimate", "std_error"), drop = FALSE]
}

expect_stata_result <- function(
  object,
  estimator,
  statistic,
  terms,
  case = "default",
  tolerance = 5e-4
) {
  reference <- stata_reference(estimator, statistic, terms, case)

  testthat::expect_equal(
    unname(stats::coef(object)),
    reference$estimate,
    tolerance = tolerance
  )
  testthat::expect_equal(
    unname(sqrt(diag(stats::vcov(object)))),
    reference$std_error,
    tolerance = tolerance
  )

  invisible(reference)
}
