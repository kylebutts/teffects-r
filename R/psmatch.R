## Main ----

#' Propensity-score matching
#'
#' Imputes each missing potential outcome from the outcomes of a small number
#' of nearest neighbors in the opposite treatment group, where proximity is
#' measured by estimated propensity-score distance. The number of neighbors
#' controls smoothing much like a bandwidth in nonparametric kernel regression.
#' Matching is performed with replacement.
#'
#' @inheritParams teffects_estimators
#' @param psmethod Propensity-score estimation method. Only maximum-likelihood
#'   `"logit"` and `"probit"` are currently implemented for matching.
#' @param weights Must be `NULL`; observation weights are not currently
#'   supported by matching estimators.
#' @return A `teffects_psmatch` object.
#' @references Imbens, G. W., and J. M. Wooldridge (2009). "Recent
#'   Developments in the Econometrics of Program Evaluation." *Journal of
#'   Economic Literature*, 47(1), 5--86.
#' @export
psmatch <- function(
  fml,
  data,
  psmethod = c("logit", "probit", "ipt", "cbps"),
  stat = c("ate", "att", "atc"),
  nneighbor = 1L,
  weights = NULL,
  vce = c("robust", "iid"),
  vce_nneighbor = 2L,
  caliper = NULL,
  pstolerance = 1e-5,
  psaction = c("error", "clamp", "trim"),
  generate = NULL,
  ...
) {
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_set_arg(psmethod, "strict match")
  dreamerr::check_set_value(
    psmethod,
    "strict match(logit, probit)",
    .arg_name = "psmethod",
    .message = "not implemented yet: use psmethod = \"logit\" or \"probit\"."
  )
  stat <- tolower(stat)
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_set_arg(vce, "strict match")
  dreamerr::check_set_arg(psaction, "strict match")
  dreamerr::check_value(nneighbor, "integer scalar GT{0}")
  dreamerr::check_value(vce_nneighbor, "integer scalar GT{0}")
  dreamerr::check_value(
    pstolerance,
    "numeric scalar GE{0} LT{0.5}",
    .arg_name = "pstolerance",
  )
  if (!is.null(caliper)) {
    dreamerr::check_arg(caliper, "class(numeric, integer)")
  }
  dreamerr::check_arg(caliper, "NULL | numeric scalar GT{0} LT{Inf}")
  dreamerr::check_arg(generate, "NULL | character scalar")
  if (!is.null(weights)) {
    stop(
      "`weights` are not currently supported by matching estimators.",
      call. = FALSE
    )
  }
  call <- match.call()
  prepared <- .teffects_prepare(
    fml,
    data,
    outcome_design = "none",
    treatment_design = "full"
  )
  y <- prepared$outcome
  treated <- prepared$treated_levels
  treatment_levels <- c(prepared$control_level, treated)
  if (length(treatment_levels) != 2L) {
    stop("`psmatch()` requires exactly two treatment levels.", call. = FALSE)
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
  treatment_candidate <- prepared$treatment_design
  row_map <- which(prepared$complete)
  repeat {
    treatment_design <- .teffects_independent_columns(
      treatment_candidate
    )
    treatment_model <- .teffects_glm(
      treatment_design,
      treatment,
      psmethod
    )
    raw_propensity <- unname(treatment_model$fitted.values)
    overlap_failure <- .teffects_overlap_check(
      raw_propensity,
      pstolerance,
      error = identical(psaction, "error")
    )
    if (!identical(psaction, "trim") || !any(overlap_failure)) {
      break
    }
    keep <- !overlap_failure
    if (length(unique(treatment[keep])) != 2L) {
      stop(
        "Trimming the propensity score removed an entire treatment level.",
        call. = FALSE
      )
    }
    row_map <- row_map[keep]
    treatment_candidate <- treatment_candidate[keep, , drop = FALSE]
    y <- y[keep]
    treatment <- treatment[keep]
  }
  propensity <- raw_propensity
  if (identical(psaction, "clamp")) {
    propensity <- .teffects_clamp_probability(propensity, pstolerance)
    if (pstolerance == 0 || any(!is.finite(propensity))) {
      .teffects_overlap_check(propensity, 0)
    }
  }
  propensity_distance <- structure(
    list(values = as.numeric(propensity), caliper = caliper),
    class = "teffects_scalar_distance"
  )
  index <- .teffects_matching_index(
    propensity_distance,
    treatment
  )
  matching <- .teffects_build_match_plan(
    index = index,
    treatment = treatment,
    stat = stat,
    nneighbor = nneighbor,
    generate = generate,
    row_map = row_map
  )
  predictions <- .teffects_predict_matches(
    y,
    treatment,
    matching
  )
  estimate <- mean(predictions$effect[matching$target])
  same_plan <- if (identical(vce, "robust")) {
    .teffects_build_neighborhood_plan(
      index,
      vce_nneighbor,
      same = TRUE,
      include_self = TRUE
    )
  } else if (identical(vce, "iid")) {
    NULL
  }
  ps_adjustment <- if (identical(vce, "iid")) {
    list(adjustment = 0)
  } else if (identical(vce, "robust")) {
    .teffects_psmatch_adjustment(
      y = y,
      treatment = treatment,
      propensity = propensity,
      treatment_model = treatment_model,
      design = treatment_design,
      stat = stat,
      estimate = estimate,
      vce_nneighbor = vce_nneighbor,
      nneighbor = nneighbor,
      index = index,
      same_plan = same_plan,
      clamped = identical(psaction, "clamp") &
        overlap_failure
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
    ps_adjustment = ps_adjustment,
    same_neighbor_offset = 0L,
    predictions = predictions,
    estimate = estimate,
    same_plan = same_plan,
    metadata = list(
      estimator = "psmatch",
      treatment_model_name = psmethod,
      treatment_model = treatment_model,
      control = control,
      treated = treated,
      sample = row_map,
      osample = .teffects_osample(prepared$n_original, row_map),
      formula = fml,
      call = call,
      nneighbor = nneighbor,
      vce_nneighbor = vce_nneighbor,
      propensity_score = propensity,
      propensity_score_raw = raw_propensity,
      psaction = psaction,
      generated_matches = matching$generated_matches
    )
  )
}

#' @export
print.teffects_psmatch <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Treatment-effects estimation\n")
  cat("Estimator: propensity-score matching\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
