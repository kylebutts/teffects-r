# Confidence intervals for treatment-effect estimates

Computes normal-approximation confidence intervals from the estimated
coefficients and covariance matrix of a treatment-effects model.

## Usage

``` r
# S3 method for class 'teffects'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A fitted `teffects` object.

- parm:

  Optional numeric or character vector selecting coefficients.

- level:

  Confidence level as a proportion in `(0, 1)`.

- ...:

  Additional arguments passed to
  [`stats::confint.default()`](https://rdrr.io/r/stats/confint.html).

## Value

A matrix with one row per selected coefficient and columns giving the
lower and upper confidence limits.
