# Mark a treatment variable in a treatment-effects formula

`treat()` identifies the treatment assignment and its reference level in
a `teffects` formula. The estimators parse the call directly, so
`treat()` does not accept any arguments beyond the treatment variable
and its reference level.

## Usage

``` r
treat(treatment, ref = NULL)
```

## Arguments

- treatment:

  A treatment variable.

- ref:

  The reference (control) treatment value. When `NULL`, estimators use
  `0` or `FALSE` when present and otherwise use the first observed
  value.

## Value

The result of
[`fixest::i()`](https://lrberge.github.io/fixest/reference/i.html) with
the requested reference level. This is only used when `treat()` is
evaluated directly; within `teffects` formulas the call is parsed as a
treatment marker before model-matrix construction.
