## Main ----

#' Nearest-neighbor matching
#'
#' Imputes each missing potential outcome from the outcomes of a small number
#' of nearest neighbors in the opposite treatment group, where proximity is
#' measured by covariate distance. The number of neighbors controls smoothing
#' much like a bandwidth in nonparametric kernel regression. Matching is
#' performed with replacement; exact matching, calipers, and regression bias
#' adjustment can refine the comparison sets.
#'
#' @inheritParams teffects_estimators
#' @param weights Must be `NULL`; observation weights are not currently
#'   supported by matching estimators.
#' @return A `teffects_nnmatch` object.
#' @references Imbens, G. W., and J. M. Wooldridge (2009). "Recent
#'   Developments in the Econometrics of Program Evaluation." *Journal of
#'   Economic Literature*, 47(1), 5--86.
#' @export
nnmatch <- function(
  fml,
  data,
  stat = c("ate", "att", "atc"),
  nneighbor = 1L,
  biasadj = NULL,
  ematch = NULL,
  weights = NULL,
  vce = c("robust", "iid"),
  vce_nneighbor = 2L,
  caliper = NULL,
  dtolerance = sqrt(.Machine$double.eps),
  generate = NULL,
  metric = c("mahalanobis", "ivariance", "euclidean", "matrix"),
  metric_matrix = NULL,
  ...
) {
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_set_arg(vce, "strict match")
  dreamerr::check_set_arg(metric, "strict match")
  dreamerr::check_value(nneighbor, "integer scalar GT{0}")
  dreamerr::check_value(vce_nneighbor, "integer scalar GT{0}")
  dreamerr::check_arg(
    biasadj,
    "NULL | os formula | character vector len(1,) no na"
  )
  dreamerr::check_arg(
    ematch,
    "NULL | os formula | character vector len(1,) no na"
  )
  if (!is.null(caliper)) {
    dreamerr::check_arg(caliper, "class(numeric, integer)")
  }
  dreamerr::check_arg(caliper, "NULL | numeric scalar GT{0} LT{Inf}")
  dreamerr::check_arg(dtolerance, "class(numeric, integer)")
  dreamerr::check_arg(dtolerance, "numeric scalar GE{0} LT{Inf}")
  dreamerr::check_arg(generate, "NULL | character scalar")
  if (!is.null(weights)) {
    stop(
      "`weights` are not currently supported by matching estimators.",
      call. = FALSE
    )
  }
  if (is.character(biasadj)) {
    biasadj <- stats::reformulate(biasadj)
  }
  if (is.character(ematch)) {
    ematch <- stats::reformulate(ematch)
  }
  call <- match.call()
  prepared <- .teffects_prepare(
    fml,
    data,
    outcome_design = "full",
    treatment_design = "none",
    additional_formulas = Filter(
      Negate(is.null),
      list(bias = biasadj, exact = ematch)
    )
  )
  row_map <- which(prepared$complete)
  y <- prepared$outcome
  treated <- prepared$treated_levels
  treatment_levels <- c(prepared$control_level, treated)
  if (length(treatment_levels) != 2L) {
    stop("`nnmatch()` requires exactly two treatment levels.", call. = FALSE)
  }
  control <- prepared$control_level
  group_sizes <- vapply(
    treatment_levels,
    function(level) sum(prepared$treatment == level),
    numeric(1)
  )
  if (nneighbor > min(group_sizes)) {
    stop(
      "`nneighbor` exceeds the size of the smaller treatment group.",
      call. = FALSE
    )
  }

  treatment <- as.integer(prepared$treatment != control)
  design <- as.matrix(prepared$outcome_design)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (!ncol(design)) {
    scaling <- matrix(numeric(), 0L, 0L)
  } else {
    scaling <- switch(
      metric,
      mahalanobis = tryCatch(
        solve(stats::cov(design)),
        error = function(error) {
          stop(
            "The matching covariates have a singular covariance matrix.",
            call. = FALSE
          )
        }
      ),
      ivariance = {
        variance <- apply(design, 2L, stats::var)
        if (any(!is.finite(variance) | variance <= 0)) {
          stop(
            "Inverse-variance matching requires nonconstant covariates.",
            call. = FALSE
          )
        }
        diag(1 / variance)
      },
      euclidean = diag(ncol(design)),
      matrix = {
        if (
          !is.matrix(metric_matrix) ||
            !all(dim(metric_matrix) == ncol(design)) ||
            !is.numeric(metric_matrix) ||
            any(!is.finite(metric_matrix)) ||
            !isSymmetric(metric_matrix)
        ) {
          stop(
            paste(
              "`metric_matrix` must be a finite symmetric square matrix",
              "matching the covariate design."
            ),
            call. = FALSE
          )
        }
        eigenvalues <- eigen(
          metric_matrix,
          symmetric = TRUE,
          only.values = TRUE
        )$values
        if (
          min(eigenvalues) <
            -sqrt(.Machine$double.eps) * max(1, max(abs(eigenvalues)))
        ) {
          stop("`metric_matrix` must be positive semidefinite.", call. = FALSE)
        }
        metric_matrix
      }
    )
  }
  distance <- .teffects_make_quadratic_distance(design, scaling, caliper)
  exact_group <- if (
    is.null(ematch) || !ncol(prepared$additional_frames$exact)
  ) {
    NULL
  } else {
    interaction(prepared$additional_frames$exact, drop = TRUE)
  }
  index <- .teffects_matching_index(
    distance,
    treatment = treatment,
    exact_group = exact_group
  )
  matching <- .teffects_build_match_plan(
    index = index,
    treatment = treatment,
    stat = stat,
    nneighbor = nneighbor,
    tolerance = dtolerance,
    generate = generate,
    row_map = row_map
  )

  bias_adjustment <- 0
  bias_models <- NULL
  bias_predictions <- NULL
  if (!is.null(biasadj)) {
    unadjusted_predictions <- .teffects_predict_matches(
      y,
      treatment,
      matching
    )
    bias_X <- .teffects_design_matrix(
      prepared$additional_terms$bias,
      prepared$additional_frames$bias,
      "full"
    )
    control_rows <- prepared$treatment == control
    treated_rows <- prepared$treatment == treated[[1L]]
    bias_models <- bias_predictions <- list()
    if (stat %in% c("att", "ate")) {
      control_X <- .teffects_independent_columns(
        bias_X[control_rows, , drop = FALSE]
      )
      bias_models$control <- .teffects_lm(
        control_X,
        y[control_rows],
        matching$reuse$k[control_rows]
      )
      bias_predictions$control <- as.numeric(
        bias_X[, colnames(control_X), drop = FALSE] %*%
          bias_models$control$coefficients
      )
    }
    if (stat %in% c("atc", "ate")) {
      treated_X <- .teffects_independent_columns(
        bias_X[treated_rows, , drop = FALSE]
      )
      bias_models$treated <- .teffects_lm(
        treated_X,
        y[treated_rows],
        matching$reuse$k[treated_rows]
      )
      bias_predictions$treated <- as.numeric(
        bias_X[, colnames(treated_X), drop = FALSE] %*%
          bias_models$treated$coefficients
      )
    }
  }
  predictions <- .teffects_predict_matches(
    y,
    treatment,
    matching,
    bias_predictions
  )
  estimate <- mean(predictions$effect[matching$target])
  if (!is.null(biasadj)) {
    bias_adjustment <- estimate -
      mean(
        unadjusted_predictions$effect[matching$target]
      )
  }
  .teffects_matching_result(
    y = y,
    treatment = treatment,
    matching = matching,
    stat = stat,
    level_names = as.character(treatment_levels),
    vce = vce,
    vce_nneighbor = vce_nneighbor,
    index = index,
    predictions = predictions,
    estimate = estimate,
    metadata = list(
      estimator = "nnmatch",
      control = control,
      treated = treated,
      sample = row_map,
      osample = .teffects_osample(prepared$n_original, row_map),
      formula = fml,
      call = call,
      nneighbor = nneighbor,
      vce_nneighbor = vce_nneighbor,
      metric = metric,
      bias_adjustment = bias_adjustment,
      bias_models = bias_models,
      generated_matches = matching$generated_matches,
      dtolerance = dtolerance
    )
  )
}

#' @export
print.teffects_nnmatch <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Treatment-effects estimation\n")
  cat("Estimator: nearest-neighbor matching\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
