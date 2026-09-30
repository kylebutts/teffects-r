# Inverse-probability weighting

Estimates potential-outcome means and treatment effects from normalized
inverse-probability-weighted outcome means. Propensity scores may be fit
by maximum-likelihood logit or probit, inverse-probability tilting, or
covariate-balancing propensity scores.

## Usage

``` r
ipw(
  fml,
  data,
  psmethod = c("logit", "probit", "ipt", "cbps"),
  stat = c("ate", "att", "atc", "pomeans"),
  weights = NULL,
  vce = c("robust", "bootstrap", "jackknife"),
  pstolerance = 1e-05,
  psaction = c("error", "clamp", "trim"),
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

  Propensity-score estimation method. `"logit"` and `"probit"` use
  maximum likelihood; `"ipt"` uses inverse-probability tilting; and
  `"cbps"` uses covariate-balancing propensity scores.

- stat:

  Target statistic. Smooth estimators support `"ate"`, `"att"`, `"atc"`,
  and/or `"pomeans"`; matching estimators support `"ate"`, `"att"`, and
  `"atc"`. Values are matched case-insensitively, so `"ATE"`, `"ATT"`,
  and `"ATC"` are accepted as well.

- weights:

  Optional nonnegative observation weights. Supply a numeric vector with
  one value per row of `data`, or the name of a column in `data`.

- vce:

  Variance estimator. Mean-based smooth estimators currently use
  `"robust"`;
  [`qte()`](https://kylebutts.github.io/teffects-r/reference/qte.md)
  uses `"bootstrap"`; matching estimators support `"iid"` and
  `"robust"`.

- pstolerance:

  Propensity-score overlap tolerance. The Stata default is `1e-5`.

- psaction:

  Action to take when estimated propensity scores violate `pstolerance`.
  `"error"` reproduces Stata's behavior and stops; `"clamp"` bounds the
  scores before weighting or matching; and `"trim"` removes violating
  observations and refits the propensity model until the retained sample
  satisfies the tolerance.

- ...:

  Reserved for display and optimization controls that have not yet been
  mapped to the R interface.

## Value

A `teffects_ipw` object.

## References

Graham, B. S., C. C. D. X. Pinto, and D. Egel (2012). "Inverse
Probability Tilting for Moment Condition Models with Missing Data."
*Review of Economic Studies*, 79(3), 1053–1079.

Imai, K., and M. Ratkovic (2014). "Covariate Balancing Propensity
Score." *Journal of the Royal Statistical Society: Series B*, 76(1),
243–263.

Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025). "Covariate
Balancing and the Equivalence of Weighting and Doubly Robust Estimators
of Average Treatment Effects." arXiv:2310.18563.
