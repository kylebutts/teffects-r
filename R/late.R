#' Local average treatment effects
#'
#' `late()` estimates a local average treatment effect (LATE) under a binary
#' instrument. `estimator = "kappa"` and `estimator = "balancing"` are
#' normalized inverse-instrument-probability weighting estimators, where
#' `"kappa"` fits the instrument propensity score by maximum likelihood and
#' `"balancing"` chooses coefficients that exactly balance the instrument-model
#' covariates across instrument levels. `estimator = "ipwra"` is a doubly robust
#' inverse-probability weighted regression-adjustment estimator that additionally
#' estimates the local average treatment effect on the treated (`stat = "latt"`).
#'
#' @param fml A two-sided formula with the outcome on the left, one `treat()`
#'   term identifying the actual treatment, and covariates on the right, for
#'   example `y ~ treat(d, ref = 0) + x1 + x2`.
#' @param data A data frame containing the variables used by the estimator.
#' @param instrument The binary instrument, supplied as an unquoted column
#'   name, a column-name string, or a vector with one value per row of `data`.
#' @param instrument_fml Optional one-sided formula for the instrument
#'   propensity score. When `NULL`, the covariates in `fml` are reused.
#' @param instrument_ref Optional value identifying the instrument's reference
#'   level. By default, `0` or `FALSE` is used when observed, and otherwise the
#'   first observed value is used.
#' @param estimator Estimation method: `"kappa"` for the normalized kappa
#'   estimator, `"balancing"` for normalized covariate balancing, or `"ipwra"`
#'   for doubly robust inverse-probability weighted regression adjustment.
#' @param stat Local effect estimated by `estimator = "ipwra"`: `"late"` or
#'   `"latt"`. The weighting estimators always report `"late"`.
#' @param imodel Instrument propensity-score model: `"logit"` or `"probit"`.
#' @param pstolerance Instrument propensity-score overlap tolerance.
#'
#' @return A `teffects_late_kappa`, `teffects_late_balancing`, or
#'   `teffects_late_ipwra` object. The reported coefficients are the local
#'   effect, its reduced-form numerator, and its first-stage denominator.
#' @references Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025).
#'   "Abadie's Kappa and Weighting Estimators of the Local Average Treatment
#'   Effect." *Journal of Business & Economic Statistics*, 43(1), 164--177.
#' @references Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2022).
#'   "Doubly Robust Estimation of Local Average Treatment Effects Using Inverse
#'   Probability Weighted Regression Adjustment." arXiv:2208.01300.
#' @section Development status:
#' These estimators are experimental and intentionally not exported from the
#' package API yet. The interface and returned objects may change before public
#' release.
#' @keywords internal
#' @name late
late <- function(
  fml,
  data,
  instrument,
  instrument_fml = NULL,
  instrument_ref = NULL,
  estimator = c("kappa", "balancing", "ipwra"),
  stat = c("late", "latt"),
  imodel = c("logit", "probit"),
  weights = NULL,
  pstolerance = 1e-5
) {
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_arg(instrument_fml, "NULL os formula")
  dreamerr::check_arg(instrument_ref, "NULL scalar")
  dreamerr::check_set_arg(estimator, "strict match")
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_set_arg(imodel, "strict match")
  dreamerr::check_value(
    pstolerance,
    "numeric scalar GE{0} LT{0.5}",
    .arg_name = "pstolerance",
  )
  if (!identical(estimator, "ipwra") && !identical(stat, "late")) {
    stop(
      "`stat = \"latt\"` requires `estimator = \"ipwra\"`.",
      call. = FALSE
    )
  }

  instrument_expression <- substitute(instrument)
  dreamerr::check_set_arg(
    instrument,
    "MBT evalset vector len(data) | charin",
    .data = data,
    .choices = names(data),
    .env = parent.frame()
  )
  instrument <- .teffects_late_resolve_instrument(
    instrument,
    instrument_expression,
    data
  )

  if (identical(estimator, "ipwra")) {
    return(.teffects_late_ipwra(
      fml = fml,
      data = data,
      instrument = instrument$value,
      instrument_name = instrument$name,
      instrument_fml = instrument_fml,
      instrument_ref = instrument_ref,
      stat = stat,
      imodel = imodel,
      weights = weights,
      pstolerance = pstolerance,
      call = match.call()
    ))
  }
  .teffects_late_weighting(
    estimator = estimator,
    fml = fml,
    data = data,
    instrument = instrument$value,
    instrument_name = instrument$name,
    instrument_fml = instrument_fml,
    instrument_ref = instrument_ref,
    imodel = imodel,
    weights = weights,
    pstolerance = pstolerance,
    call = match.call()
  )
}
