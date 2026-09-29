# Marginal quantile treatment effects

Estimates quantiles of marginal potential-outcome distributions and
their contrasts under conditional independence. The implementation
follows the inverse-probability-weighted (IPW) and
efficient-influence-function (EIF) estimators in Cattaneo (2010). The
EIF estimator uses linear series regressions for the conditional
distribution functions.

## Usage

``` r
qte(
  fml,
  data,
  probs = c(0.25, 0.5, 0.75),
  treatment_fml = NULL,
  method = c("eif", "ipw"),
  tmodel = c("logit", "probit"),
  weights = NULL,
  vce = "bootstrap",
  reps = 999L,
  pstolerance = 1e-05,
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

- probs:

  Numeric vector of quantile indexes strictly between zero and one.

- treatment_fml:

  Optional one-sided formula for the treatment model. When `NULL`, the
  outcome-model covariates in `fml` are reused.

- method:

  Estimation method. `"eif"` is the augmented, efficient-
  influence-function estimator and `"ipw"` is normalized inverse-
  probability weighting.

- tmodel:

  Treatment-model family: `"logit"` or `"probit"`. Smooth estimators fit
  multinomial logit when treatment is multivalued, so multivalued
  treatment requires `"logit"`.

- weights:

  Optional nonnegative observation weights. Supply a numeric vector with
  one value per row of `data`, or the name of a column in `data`.

- vce:

  Variance estimator. Currently only `"bootstrap"` is supported.

- reps:

  Number of successful nonparametric bootstrap replications.

- pstolerance:

  Propensity-score overlap tolerance. The Stata default is `1e-5`.

- ...:

  Reserved for display and optimization controls that have not yet been
  mapped to the R interface.

## Value

A `teffects_qte` object. The `potential_quantiles` element in
`x$teffects` contains the estimated quantile for every treatment level.

## Details

The reported QTE at index `tau` is the difference between quantiles of
two marginal potential-outcome distributions. It is not, without
additional assumptions, the quantile of unit-level treatment effects.
Inference is bootstrap-based because analytic standard errors require
density estimation at every estimated quantile.

This implementation targets population quantiles for continuous
outcomes. It does not estimate quantile treatment effects on the
treated.

## References

Cattaneo, M. D. (2010). Efficient semiparametric estimation of
multi-valued treatment effects under ignorability. *Journal of
Econometrics*, 155(2), 138–154.

Cattaneo, M. D., Drukker, D. M., and Holland, A. D. (2013). Estimation
of multivalued treatment effects under conditional independence. *The
Stata Journal*, 13(3), 407–450.
