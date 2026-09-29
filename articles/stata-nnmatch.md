# Stata matching details

This vignette documents the matching rules used by Stata’s
`teffects nnmatch`. It is written for users who want to reproduce a
Stata analysis with
[`teffects::nnmatch()`](https://kylebutts.github.io/teffects-r/reference/nnmatch.md),
or who need to understand why a command requesting one neighbor can use
more than one observation.

## The basic estimator

Stata matches each observation to observations in the opposite treatment
group. It uses the matched outcomes to impute the missing potential
outcome. For an observation $`i`$, let $`\Omega(i)`$ denote its set of
matches. The imputed outcome is

``` math
\widehat{Y}_i(t) = \begin{cases}
  Y_i, & D_i = t, \\
  \displaystyle \frac{\sum_{j \in \Omega(i)} w_jY_j} {\sum_{j \in \Omega(i)} w_j}, & D_i \neq t.
\end{cases}
```

Individual treatment effects are given by
$`\widehat{Y}_{i}(1) - \widehat{Y}_{i}(0)`$. Averages of the individual
treatment effects depend on the requested `stat`: - `ate` averages over
the analysis sample - `att`/`atc` averages over the treated/control
observations.

## Stata defaults

Stata’s matching procedure looks like this. For a given treated unit: -
First, restrict the candidate pool of control units using the exact
match criteria (`ematch`) - Second, calculate distances between the
treated unit and the pool. Drop any units with distance larger than the
`caliper`. - Third, grab the closest `nneighbor` matches. In the case of
ties, include all ties.

Following Stata, matching is done with replacement. A control
observation can therefore be used as the neighbor for many treated
observations, and a treated observation can be used for many controls
when estimating the ATE.

### Number of matches and ties

The `nneighbor` argument requests the *minimum* number of matches (not a
maximum).

Stata finds the distance of the (m)-th closest eligible observation. It
then retains every eligible observation at that distance. If several
observations are tied at the cutoff, all of them are included in the
matched control group for that unit.

### Exact matching

`ematch` requires equality on the listed variables. This allows discrete
or continuous variables. For continuous variables, `dtolerance`
specifies the tolerance for an exact match. It defaults to
`sqrt(.Machine$double.eps)`.

``` r

nnmatch(
  y ~ treat(d, ref = 0) + age,
  data = dat,
  ematch = ~ female + region,
  metric = "euclidean"
)
```

### Calipers

`caliper(c)` removes candidates whose matching distance exceeds `c`. The
caliper is applied before the nearest-neighbor cutoff. It can therefore
make an observation deficient in matches even when the opposite
treatment group is large.

## Distance metrics

Stata’s default is Mahalanobis distance. If $`x_i`$ and $`x_j`$ are
covariate vectors, the distance has the form

``` math
  d(i,j) = \{(x_i-x_j)'S^{-1}(x_i-x_j)\}^{1/2},
```

where the scaling matrix $`S`$ is calculated from the sample of matching
variables. The available choices are:

| Stata option | Scaling | R [`nnmatch()`](https://kylebutts.github.io/teffects-r/reference/nnmatch.md) option |
|----|----|----|
| `metric(mahalanobis)` | Full sample covariance matrix | `metric = "mahalanobis"` |
| `metric(ivariance)` | Diagonal sample covariance matrix | `metric = "ivariance"` |
| `metric(euclidean)` | Identity matrix | `metric = "euclidean"` |
| `metric(matrix M)` | User-supplied symmetric matrix | `metric = "matrix"`, `metric_matrix = M` |

The matching variables determine the distance.

## Bias adjustment

With more than one continuous matching variable, imperfect matches can
create a biased treatment effect estimator. After matching, a regression
adjustment estimator can be used between the treated and matched control
unit to do bias adjustment. Use the `biasadj = ~ ..x` argument to do
this.

## Standard errors

There are two options for standard errors. vce = “iid” requests the iid
Abadie–Imbens standard error which typically are too small when the
number of matches is fixed (small).

Alterantively, we can use robust standard errors proposed by Abadie and
Imbens. See their work for details, but basically we look at the
variance of the error term between the focal unit and it’s neighbors of
the *same* treatment status. This means we need to do another nearest
neighbor matching procedure. We use `vce = "robust"` to use these
standard errors and `vce_nneighbor = 2L` to detemrine the number of
neighbors to match with (defaults to 2).

``` r

nnmatch(
  y ~ treat(d, ref = 0) + x1 + x2,
  data = dat,
  nneighbor = 1,
  vce = "robust",
  vce_nneighbor = 4
)
```

## References

StataCorp. (2026). *\[CAUSAL\] Causal: Treatment-effects estimation for
observational data*. See the `teffects nnmatch` sections on
nearest-neighbor matching, options, and methods and formulas.

Abadie, A., and G. W. Imbens (2006). Large sample properties of matching
estimators for average treatment effects. *Econometrica* 74(1), 235–267.

Abadie, A., and G. W. Imbens (2011). Bias-corrected matching estimators
for average treatment effects. *Journal of Business & Economic
Statistics* 29(1), 1–11.
