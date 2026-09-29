# Inverse-probability-weighted regression adjustment

Estimates potential-outcome means and treatment effects by fitting a
propensity-score-weighted outcome regression in each treatment arm.

## Usage

``` r
ipwra(
  fml,
  data,
  treatment_fml = NULL,
  omodel = c("linear", "logit", "probit", "poisson"),
  psmethod = c("logit", "probit", "ipt", "cbps"),
  stat = c("ate", "att", "atc", "pomeans"),
  weights = NULL,
  vce = c("robust", "bootstrap", "jackknife"),
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

- treatment_fml:

  Optional one-sided formula for the treatment model. When `NULL`, the
  outcome-model covariates in `fml` are reused.

- omodel:

  Outcome-model family. Only `"linear"` is currently implemented; the
  other accepted names are reserved for future releases.

- psmethod:

  Propensity-score estimation method. `"logit"` and `"probit"` use
  maximum likelihood; `"ipt"` uses inverse-probability tilting; and
  `"cbps"` uses covariate-balancing propensity scores.

- stat:

  Target statistic. Smooth estimators support `"ate"`, `"att"`, `"atc"`,
  and/or `"pomeans"`; matching estimators support `"ate"`, `"att"`, and
  `"atc"`.

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

- ...:

  Reserved for display and optimization controls that have not yet been
  mapped to the R interface.

## Value

A `teffects_ipwra` object.
