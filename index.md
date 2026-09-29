# teffects

**teffects** provides treatment-effect estimators with a consistent R
formula interface modeled on Stata’s `teffects` suite. It currently has
tests covering the complete set of examples from the excellent(!)
[Stata’s teffects manual](https://www.stata.com/manuals16/te.pdf).
Really, you should open that manual and read the introduction; it’s
wonderuflly written!

To install this package, use:

``` r

devtools::install_github("kylebutts/teffects-r")
```

## From Stata to R

This package supports regression adjustment, inverse-probability
weighting, covariate balancing, doubly robust, propensity-score
matching, and nearest-neighbor matching estimators. The estimator names
follow Stata’s `teffects` commands:

| Stata                  | R                                  |
|------------------------|------------------------------------|
| `teffects ra ...`      | `ra(...)`                          |
| `teffects ipw ...`     | `ipw(..., psmethod = "logit")`     |
| `teffects aipw ...`    | `aipw(..., psmethod = "logit")`    |
| `teffects ipwra ...`   | `ipwra(..., psmethod = "logit")`   |
| `teffects psmatch ...` | `psmatch(..., psmethod = "logit")` |
| `teffects nnmatch ...` | `nnmatch(...)`                     |

In R, `...` will be the shared outcome/treatment formula (following
fixest formula symantics) and the `data` argument. You can use
`treatment_fml` when the propensity-score model needs a different
covariate specification from the outcome model. For balance-first
weighting, use `ipw(..., psmethod = "ipt")` or
`ipw(..., psmethod = "cbps")`; the same propensity-score methods are
available in
[`aipw()`](https://kylebutts.github.io/teffects-r/reference/aipw.md) and
[`ipwra()`](https://kylebutts.github.io/teffects-r/reference/ipwra.md).

Use
[`treat()`](https://kylebutts.github.io/teffects-r/reference/treat.md)
inside the formula to identify treatment and its reference level:

``` r

y ~ treat(d, ref = 0) + x1 + x2
```

The outcome appears on the left, while the remaining right-hand-side
terms are pre-treatment adjustment covariates that we will condition on.

``` r

library(teffects)
library(fixest)

set.seed(2026)
n <- 300
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.45)
p <- plogis(-0.2 + 0.8 * x1 - 0.6 * x2)
d <- rbinom(n, 1, p)
y0 <- 1 + 1.5 * x1 - 0.8 * x2 + rnorm(n)
example_data <- data.frame(y = y0 + 2 * d, d, x1, x2)

setFixest_fml(..x = ~ x1 + x2)

## a few examples
teffects::ra(y ~ treat(d, ref = 0) + ..x, data = example_data)
#> Treatment-effects estimation
#> Estimator: regression adjustment
#> Outcome model: linear
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  1.93419    0.12682 15.2521 < 2.2e-16 ***
#> POmean[0]    0.68825    0.11796  5.8345 5.397e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
teffects::ipw(y ~ treat(d, ref = 0) + ..x, data = example_data)
#> Treatment-effects estimation
#> Estimator: inverse-probability weighting
#> Treatment model: logit 
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  2.09548    0.13205 15.8692 < 2.2e-16 ***
#> POmean[0]    0.64094    0.11555  5.5469 2.907e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
teffects::ipw(
  y ~ treat(d, ref = 0) + ..x,
  data = example_data,
  psmethod = "cbps"
)
#> Treatment-effects estimation
#> Estimator: inverse-probability weighting
#> Treatment model: cbps 
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  1.95449    0.13440 14.5422 < 2.2e-16 ***
#> POmean[0]    0.70616    0.11789  5.9899 2.099e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
teffects::aipw(y ~ treat(d, ref = 0) + ..x, data = example_data)
#> Treatment-effects estimation
#> Estimator: augmented inverse-probability weighting
#> Outcome model: linear
#> Propensity-score method: logit 
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  1.94881    0.13146 14.8245 < 2.2e-16 ***
#> POmean[0]    0.69214    0.11906  5.8133 6.126e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Estimator families

Start with the
[Introduction](https://kylebutts.github.io/teffects-r/articles/introduction.md),
then use [Mathematical
details](https://kylebutts.github.io/teffects-r/articles/mathematical-details.md)
and the [Stata matching
details](https://kylebutts.github.io/teffects-r/articles/stata-nnmatch.md)
article for derivations and implementation details.

### Traditional estimators

First, there are the general `teffects` estimators. These are the
standard treatment effect estimators under the conditional independence
assumption.

| Goal | Functions | Approach |
|----|----|----|
| Model conditional outcomes | [`ra()`](https://kylebutts.github.io/teffects-r/reference/ra.md) | Fit an outcome model in each treatment arm and standardize predictions |
| Reweight observed outcomes | [`ipw()`](https://kylebutts.github.io/teffects-r/reference/ipw.md) | Estimate generalized propensity scores by likelihood, IPT, or CBPS and form normalized weighted means |
| Match comparable units | [`psmatch()`](https://kylebutts.github.io/teffects-r/reference/psmatch.md), [`nnmatch()`](https://kylebutts.github.io/teffects-r/reference/nnmatch.md) | Impute missing potential outcomes from a few nearest neighbors in the opposite treatment group |
| Combine nuisance models | [`aipw()`](https://kylebutts.github.io/teffects-r/reference/aipw.md), [`ipwra()`](https://kylebutts.github.io/teffects-r/reference/ipwra.md) | Combine outcome regression with treatment-probability weighting |

### Machine-learning estimators

The high-dimensional workflow uses rigorous plugin lasso to select
controls, then refits the selected nuisance models without a penalty.
Cross-fitting is available when out-of-sample nuisance predictions are
desired.

| Goal | Functions | Approach |
|----|----|----|
| Select high-dimensional controls | [`telasso()`](https://kylebutts.github.io/teffects-r/reference/telasso.md) | Use plugin lasso, post-selection refits, and optional cross-fitting in an orthogonal AIPW estimator |

### Quantile treatment-effects

Here, a quantile treatment effect is a shift in marginal quantiles—for
example, the median of $`Y_i(1)`$ minus the median of $`Y_i(0)`$. It is
not generally the median of $`\tau_i = Y_i(1) - Y_i(0)`$; that
interpretation requires strong rank-invariance assumptions.

| Goal | Functions | Approach |
|----|----|----|
| Compare marginal outcome quantiles | [`qte()`](https://kylebutts.github.io/teffects-r/reference/qte.md) | Estimate IPW or efficient-influence-function quantile treatment effects with bootstrap inference |

## Supported targets

Smooth estimators support average treatment effects, effects on treated
units (`ATT`), effects on control units (`ATC`), potential-outcome
means, and marginal quantile treatment effects where applicable.
Matching estimators support ATE, ATT, and ATC. Several smooth estimators
also support multivalued treatment.

Identification still depends on the research design. The package cannot
test whether all confounders were measured, whether variables were
recorded before treatment, or whether a causal contrast is substantively
meaningful. Use the modeling tools together with careful covariate
selection and overlap checks.
