# Package index

## Regression adjustment

Fit ordinary or propensity-weighted outcome regressions and standardize
their predictions.

- [`ra()`](https://kylebutts.github.io/teffects-r/reference/ra.md) :
  Regression adjustment

## Weighting and balancing

Model treatment assignment or solve balance equations.

- [`ipw()`](https://kylebutts.github.io/teffects-r/reference/ipw.md) :
  Inverse-probability weighting

## Matching

Impute missing potential outcomes from a few opposite-group neighbors.

- [`psmatch()`](https://kylebutts.github.io/teffects-r/reference/psmatch.md)
  : Propensity-score matching
- [`nnmatch()`](https://kylebutts.github.io/teffects-r/reference/nnmatch.md)
  : Nearest-neighbor matching

## Doubly robust estimation

Combine treatment and outcome models.

- [`aipw()`](https://kylebutts.github.io/teffects-r/reference/aipw.md) :
  Augmented inverse-probability weighting
- [`ipwra()`](https://kylebutts.github.io/teffects-r/reference/ipwra.md)
  : Inverse-probability-weighted regression adjustment
- [`telasso()`](https://kylebutts.github.io/teffects-r/reference/telasso.md)
  : Treatment-effects estimation using lasso

## Quantile treatment effects

Compare quantiles of marginal potential-outcome distributions.

- [`qte()`](https://kylebutts.github.io/teffects-r/reference/qte.md) :
  Marginal quantile treatment effects

## Formula interface

- [`treat()`](https://kylebutts.github.io/teffects-r/reference/treat.md)
  : Mark a treatment variable in a treatment-effects formula

## Methods

- [`confint(`*`<teffects>`*`)`](https://kylebutts.github.io/teffects-r/reference/confint.teffects.md)
  : Confidence intervals for treatment-effect estimates
