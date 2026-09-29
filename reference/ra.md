# Regression adjustment

Estimates potential-outcome means and treatment effects by fitting a
separate linear outcome regression in each treatment arm and averaging
its predictions over the covariate distribution of the target
population.

## Usage

``` r
ra(
  fml,
  data,
  omodel = c("linear", "logit", "probit", "poisson"),
  stat = c("ate", "att", "atc", "pomeans"),
  weights = NULL,
  vce = c("robust", "bootstrap", "jackknife"),
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

- omodel:

  Outcome-model family. Only `"linear"` is currently implemented; the
  other accepted names are reserved for future releases.

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

- ...:

  Reserved for display and optimization controls that have not yet been
  mapped to the R interface.

## Value

A `teffects_ra` object.

## References

Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025). "Covariate
Balancing and the Equivalence of Weighting and Doubly Robust Estimators
of Average Treatment Effects." arXiv:2310.18563.
