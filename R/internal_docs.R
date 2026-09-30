#' Treatment-effect estimators
#'
#' These functions build sparse model matrices and implement the estimators
#' directly. They share a formula interface in which `treat()` marks
#' the binary or multivalued treatment variable, for example
#' `y ~ treat(D, ref = 0) + x1 + x2`.
#'
#' The argument names and defaults follow the corresponding Stata
#' `teffects_*.sthlp` files where an R analogue is meaningful.
#'
#' @param fml A two-sided formula with the outcome on the left-hand side and a
#'   single `treat()` term on the right-hand side. Use its `ref` argument to
#'   select the control level. Additional right-hand-side terms specify
#'   covariates.
#' @param data A data frame containing the variables in `fml`.
#' @param treatment_fml Optional one-sided formula for the treatment model.
#'   When `NULL`, the outcome-model covariates in `fml` are reused.
#' @param omodel Outcome-model family. Only `"linear"` is currently
#'   implemented; the other accepted names are reserved for future releases.
#' @param tmodel Treatment-model family: `"logit"` or `"probit"`. Smooth
#'   estimators fit multinomial logit when treatment is multivalued, so
#'   multivalued treatment requires `"logit"`.
#' @param psmethod Propensity-score estimation method. `"logit"` and
#'   `"probit"` use maximum likelihood; `"ipt"` uses inverse-probability
#'   tilting; and `"cbps"` uses covariate-balancing propensity scores.
#' @param stat Target statistic. Smooth estimators support `"ate"`, `"att"`,
#'   `"atc"`, and/or `"pomeans"`; matching estimators support `"ate"`, `"att"`,
#'   and `"atc"`. Values are matched case-insensitively, so `"ATE"`, `"ATT"`,
#'   and `"ATC"` are accepted as well.
#' @param cme Conditional-mean estimation method for AIPW. Only maximum
#'   likelihood (`"ml"`) is currently implemented.
#' @param weights Optional nonnegative observation weights. Supply a numeric
#'   vector with one value per row of `data`, or the name of a column in `data`.
#' @param vce Variance estimator. Mean-based smooth estimators currently use
#'   `"robust"`; `qte()` uses `"bootstrap"`; matching estimators support
#'   `"iid"` and `"robust"`.
#' @param vce_nneighbor Number of matches used to estimate the robust
#'   Abadie--Imbens variance for a matching estimator. Stata defaults to two.
#' @param pstolerance Propensity-score overlap tolerance. The Stata default is
#'   `1e-5`.
#' @param psaction Action to take when estimated propensity scores violate
#'   `pstolerance`. `"error"` reproduces Stata's behavior and stops;
#'   `"clamp"` bounds the scores before weighting or matching; and `"trim"`
#'   removes violating observations and refits the propensity model until the
#'   retained sample satisfies the tolerance.
#' @param nneighbor Number of opposite-treatment nearest neighbors used to
#'   impute a missing potential outcome. It is the matching estimator's
#'   smoothing parameter and defaults to one.
#' @param biasadj Optional variables used for the nearest-neighbor large-sample
#'   bias adjustment.
#' @param ematch Optional variables on which nearest-neighbor matching must be
#'   exact.
#' @param caliper Optional maximum distance between potential matches. `NULL`
#'   imposes no caliper.
#' @param dtolerance Maximum distance at which covariate values are considered
#'   equal for nearest-neighbor matching.
#' @param generate Optional prefix used for retained nearest-neighbor indices.
#' @param metric Nearest-neighbor covariate-distance metric: `"mahalanobis"`,
#'   `"ivariance"`, `"euclidean"`, or `"matrix"`.
#' @param metric_matrix User-supplied scaling matrix when
#'   `metric = "matrix"`.
#' @param ... Reserved for display and optimization controls that have not yet
#'   been mapped to the R interface.
#'
#' @keywords internal
#' @name teffects_estimators
NULL
