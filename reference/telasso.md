# Treatment-effects estimation using lasso

Estimates average treatment effects and potential-outcome means using
augmented inverse-probability weighting after plugin-lasso selection of
high-dimensional controls. Outcome controls are selected separately
within each treatment arm and all selected models are refitted without a
penalty.

## Usage

``` r
telasso(
  fml,
  data,
  treatment_fml = NULL,
  omodel = "linear",
  tmodel = "logit",
  stat = c("ate", "att", "atc", "pomeans"),
  selection = "plugin",
  xfolds = 1L,
  resamples = 1L,
  seed = NULL,
  outcome_ainclude = NULL,
  treatment_ainclude = NULL,
  weights = NULL,
  vce = "robust",
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

- tmodel:

  Treatment-model family: `"logit"` or `"probit"`. Smooth estimators fit
  multinomial logit when treatment is multivalued, so multivalued
  treatment requires `"logit"`.

- stat:

  Target statistic. Smooth estimators support `"ate"`, `"att"`, `"atc"`,
  and/or `"pomeans"`; matching estimators support `"ate"`, `"att"`, and
  `"atc"`.

- selection:

  Penalty-selection method. The rigorous plugin method is currently
  implemented.

- xfolds:

  Number of folds used for cross-fitting. One, the default, estimates
  the nuisance functions on the full estimation sample.

- resamples:

  Number of independently generated cross-fitting partitions to average.
  Values above one require `xfolds > 1`.

- seed:

  Optional integer seed used to generate cross-fitting partitions.

- outcome_ainclude, treatment_ainclude:

  Optional character vectors naming expanded model-matrix columns that
  must be retained in the corresponding post-lasso model. The intercept,
  when present, is always retained.

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

A `teffects_telasso` object.
