# Nearest-neighbor matching

Imputes each missing potential outcome from the outcomes of a small
number of nearest neighbors in the opposite treatment group, where
proximity is measured by covariate distance. The number of neighbors
controls smoothing much like a bandwidth in nonparametric kernel
regression. Matching is performed with replacement; exact matching,
calipers, and regression bias adjustment can refine the comparison sets.

## Usage

``` r
nnmatch(
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

- stat:

  Target statistic. Smooth estimators support `"ate"`, `"att"`, `"atc"`,
  and/or `"pomeans"`; matching estimators support `"ate"`, `"att"`, and
  `"atc"`.

- nneighbor:

  Number of opposite-treatment nearest neighbors used to impute a
  missing potential outcome. It is the matching estimator's smoothing
  parameter and defaults to one.

- biasadj:

  Optional variables used for the nearest-neighbor large-sample bias
  adjustment.

- ematch:

  Optional variables on which nearest-neighbor matching must be exact.

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

- dtolerance:

  Maximum distance at which covariate values are considered equal for
  nearest-neighbor matching.

- generate:

  Optional prefix used for retained nearest-neighbor indices.

- metric:

  Nearest-neighbor covariate-distance metric: `"mahalanobis"`,
  `"ivariance"`, `"euclidean"`, or `"matrix"`.

- metric_matrix:

  User-supplied scaling matrix when `metric = "matrix"`.

- ...:

  Reserved for display and optimization controls that have not yet been
  mapped to the R interface.

## Value

A `teffects_nnmatch` object.

## References

Imbens, G. W., and J. M. Wooldridge (2009). "Recent Developments in the
Econometrics of Program Evaluation." *Journal of Economic Literature*,
47(1), 5–86.
