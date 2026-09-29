# Propensity-score matching

Imputes each missing potential outcome from the outcomes of a small
number of nearest neighbors in the opposite treatment group, where
proximity is measured by estimated propensity-score distance. The number
of neighbors controls smoothing much like a bandwidth in nonparametric
kernel regression. Matching is performed with replacement.

## Usage

``` r
psmatch(
  fml,
  data,
  psmethod = c("logit", "probit", "ipt", "cbps"),
  stat = c("ate", "att", "atc"),
  nneighbor = 1L,
  weights = NULL,
  vce = c("robust", "iid"),
  vce_nneighbor = 2L,
  caliper = NULL,
  pstolerance = 1e-05,
  psaction = c("error", "clamp", "trim"),
  generate = NULL,
  ...
)
```

## Arguments

- fml:

  A two-sided formula with the outcome on the left-hand side and a
  single
  [`treat()`](https://kylebutts.github.io/teffects-r/reference/treat.md)
  term on the right-hand side. Use its `ref` argument to select the
  control level. Additional right-hand-side terms specify covariates.

- data:

  A data frame containing the variables in `fml`.

- psmethod:

  Propensity-score estimation method. Only maximum-likelihood `"logit"`
  and `"probit"` are currently implemented for matching.

- stat:

  Target statistic. Smooth estimators support `"ate"`, `"att"`, `"atc"`,
  and/or `"pomeans"`; matching estimators support `"ate"`, `"att"`, and
  `"atc"`.

- nneighbor:

  Number of opposite-treatment nearest neighbors used to impute a
  missing potential outcome. It is the matching estimator's smoothing
  parameter and defaults to one.

- weights:

  Must be `NULL`; observation weights are not currently supported by
  matching estimators.

- vce:

  Variance estimator. Mean-based smooth estimators currently use
  `"robust"`;
  [`qte()`](https://kylebutts.github.io/teffects-r/reference/qte.md)
  uses `"bootstrap"`; matching estimators support `"iid"` and
  `"robust"`.

- vce_nneighbor:

  Number of matches used to estimate the robust Abadie–Imbens variance
  for a matching estimator. Stata defaults to two.

- caliper:

  Optional maximum distance between potential matches. `NULL` imposes no
  caliper.

- pstolerance:

  Propensity-score overlap tolerance. The Stata default is `1e-5`.

- psaction:

  Action to take when estimated propensity scores violate `pstolerance`.
  `"error"` reproduces Stata's behavior and stops; `"clamp"` bounds the
  scores before weighting or matching; and `"trim"` removes violating
  observations and refits the propensity model until the retained sample
  satisfies the tolerance.

- generate:

  Optional prefix used for retained nearest-neighbor indices.

- ...:

  Reserved for display and optimization controls that have not yet been
  mapped to the R interface.

## Value

A `teffects_psmatch` object.

## References

Imbens, G. W., and J. M. Wooldridge (2009). "Recent Developments in the
Econometrics of Program Evaluation." *Journal of Economic Literature*,
47(1), 5–86.
