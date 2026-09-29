# Introduction

The `teffects` package implements several estimators for causal effects
in observational data. They all start from the same problem and use the
same underlying causal assumption: conditional independence. They differ
in how they use the observed data to form estimates.

This vignette explains that logic first. It then shows the common
formula interface and gives a short guide to the estimator families. The
mathematical derivations and inference are discussed in [Mathematical
details](https://kylebutts.github.io/teffects-r/articles/mathematical-details.md).

For treatment levels $`d \in \{0, \ldots, K \}`$, let
$`\mu_d = \mathbb{E}\{Y_i(d)\}`$ be the average value of the $`d`$-th
potential outcome for the population. Then, (overall) average treatment
effects are defined as
``` math
  ATE(d) = \mu_d - \mu_0.
```
In the case of binary treatments, this is the $`\text{ATE}`$. This
quantity averages over every unit in the population to form the overall
average.

Sometimes, we want to know the average effect of treatment on the units
that actually select into treatment. We define this as the average
treatment effect on the treated. For a binary treatment, the ATT is
``` math
  ATT = \mathbb{E}\{ Y_i(1) - Y_i(0) \ \mid \ D_i = 1\}.
```

In general, we have $`K`$ ATT parameters, one for each treatment value:
``` math
  ATT(d \ \mid \  d) = \mathbb{E}\{ Y_i(1) - Y_i(0) \ \mid \ D_i = d \}.
```

We will make two main assumptions to identify these treatment effect
parameters: conditional independence and overlap.

## The problem with a difference in means

Let $`D_i`$ denote treatment, $`X_i`$ the pre-treatment covariates, and
$`Y_i(d)`$ the potential outcome under treatment level $`d`$. We observe
only

``` math
  Y_i = \sum_{d} 1(D_i = d) Y_i(d).
```

The fundamental problem of causal inference is that we only observe a
single potential outcome for each unit. In a randomized experiment,
treatment assignment is independent of potential outcomes, so the
control group provides a good estimate of what would have happened to
treated units had they not been treated.

Most empirical work in economics uses observational data. Units select
into treatment, so treated and untreated units can differ in ways that
matter for the outcome. For example, college graduates may have higher
school smarts, more educated parents, and better high schools. Those
differences would lead to higher earnings even without college, so a
comparison of average earnings does not isolate the effect of attending
college.

For a binary treatment, the difference in observed means can be written
as

``` math
E(Y_i \ \mid \ D_i = 1) - E(Y_i \ \mid\ D_i = 0)
  = \underbrace{\mathbb{E}\{Y_i(1) - Y_i(0) \ \mid \  D_i=1\}}_{\text{ATT}}
  + \underbrace{\mathbb{E}\{Y_i(0) \ \mid \  D_i = 1\} - \mathbb{E}\{Y_i(0) \ \mid \  D_i = 0\}}_{\text{selection bias}}.
```

The second term is the problem and is why just comapring average
outcomes is not sufficient to identify treatment effects. Our conditiona
independence assumption will allow us to remove this term by comparing
units with the same relevant observed characteristics.

## Selection on observables

The key assumption is the conditional independence assumption, also
called the selection-on-observables assumption:

``` math
  ( Y_i(0), Y_i(1) )\mathrel{\perp\!\!\!\perp} D_i \ \mid \ X_i = x.
```

Among units with the same $`X_i=x`$, treatment assignment is independent
of the potential outcomes. In effect, treatment is “as good as randomly
assigned” within each covariate cell. This is like having many small
randomized experiments, one for each value of $`x`$.

There is always going to be debate about what variables are missing in
$`X_i`$. As I put it in my class:

> Author: “X_i includes a lot of important factors that drive selection
> into treatment!” Reviewer: “I think there are other omitted variables
> that drive selection into treatment!!”

The credibility of a selection-on-observables design depends on that
debate.

## From conditional comparisons to causal effects

Suppose first that $`X_i`$ is a discrete variable such as gender. Under
conditional independence assumption, we have two (as good as) random
experiments: one for men and one for women.

The solution is to estimate a treatment effect within each group:
``` math
 \tau(x) = E(Y_i \ \mid \ D_i = 1, X_i = x) - E(Y_i \ \mid \ D_i = 0, X_i =  x).
```
That is, we will estimate conditional average treatment effects for each
value of the covariate using difference in means.

Under conditional independence, this is the conditional average
treatment effect,
``` math
  \tau(x) = \mathbb{E} \{ Y_i(1) - Y_i(0) \ \mid \ X_i = x\}.
```

The overall effect is an average of these conditional effects
``` math
  ATE = \int \tau(x) \, dP(X = x).
```

This approach of performing many difference-in-means estimators works
well when $`X_i`$ is a single or a handful of variables, but quickly
becomes infeasible as $`X_i`$ becomes larger. In practice, trying to
find people with the exact same value of $`X_i`$ when you have many
control variables becomes too difficult to produce any good estimates.
The problem is called the **curse of dimensionality**.

Therefore, we will need to reach for alternative estimators. The goal of
this package is to make these available to users.

## Propensity Score

Above we discuss the difficulties with controlling for a large number of
$`X_i`$ variables. It turns out that under conditional independence, it
is sufficient to control for a single scalar variable: the propensity
score.

That is, under conditional independence, we have:
``` math
  (Y_i(0), Y_i(1)) \mathrel{\perp\!\!\!\perp} D_i \mid \pi(X_i).
```

This reduces a high-dimensional comparison ($`X_i`$) to one treatment
probability ($`\pi(X_i)`$).

## Overlap

In addition to the conditional independence assumption, we will need the
**overlap assumption**.

The overlap assumption involves the **propensity score** which is the
probability of being in the treatment group given your covariate values:
``` math
  \pi(x) = P(D_i = 1 \ \mid\  X_i = x)
```

The overlap condition requires
``` math
  0 < P(D_i = 1 \ \mid\  X_i = x) < 1
```
for all values of $`x`$. Without overlap, we will not be able to find
treated and untreated observations with the same value of $`x`$ to
compare.

In practice is it common to trim observations with poor overlap
(propensity score close to 0.), but it’s important to remember that this
changes the target population (ATE among units satisfying overlap).

## Alternative treatment effect estimators

The estimators below all are based on the conditional independence
assumption and the overlap assumption. They just differ in how they
compute treatment effects under these assumptions.

Here is a high-level summary of the estimators, but we will go into
details on each one below:

| Function | Adjustment strategy | Main modeling task |
|----|----|----|
| [`ra()`](https://kylebutts.github.io/teffects-r/reference/ra.md) | Regression adjustment | Fit an outcome model in each treatment arm |
| [`ipwra()`](https://kylebutts.github.io/teffects-r/reference/ipwra.md) | IPW regression adjustment | Fit propensity-weighted outcome models |
| [`ipw()`](https://kylebutts.github.io/teffects-r/reference/ipw.md) | Inverse-probability weighting | Estimate treatment probabilities or solve balance equations |
| [`aipw()`](https://kylebutts.github.io/teffects-r/reference/aipw.md) | Augmented IPW | Combine outcome and treatment models |
| [`psmatch()`](https://kylebutts.github.io/teffects-r/reference/psmatch.md) | Propensity-score matching | Match on estimated treatment probabilities |
| [`nnmatch()`](https://kylebutts.github.io/teffects-r/reference/nnmatch.md) | Nearest-neighbor matching | Match directly on covariates |

For this vignette, we will focus on estimating the $`\text{ATT}`$, but
the intuition will stand for $`\text{ATE}`$ parameters. The only
practical difference is that the $`\text{ATE}`$ requires us to guess
both the treated units’ untreated potential outcome *and* the control
units’ treated potential outcome. To begin, we will generate example
data to use in our discussion:

``` r

library(teffects)

set.seed(2026)
n <- 300
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.45)
probability <- plogis(-0.2 + 0.8 * x1 - 0.6 * x2)
d <- rbinom(n, 1, probability)
y0 <- 1 + 1.5 * x1 - 0.8 * x2 + rnorm(n)
example_data <- data.frame(y = y0 + 2 * d, d, x1, x2)
```

All estimators mark treatment with
[`treat()`](https://kylebutts.github.io/teffects-r/reference/treat.md):

``` r

y ~ treat(d, ref = 0) + x1 + x2
```

The left side is the outcome. `ref` selects the comparison level. The
remaining terms describe pre-treatment adjustment covariates. Use
ordinary R formula syntax for factors, transformations, and
interactions. Use `treatment_fml` when the treatment model needs a
different covariate specification from the outcome model. Smooth
estimators use `stat = "ate"` by default. Use `stat = "pomeans"` for all
potential-outcome means, `stat = "att"` for the treated population, or
`stat = "atc"` for the control population when the estimator supports
that target.

### Regression adjustment

With many covariates or continuous $`X_i`$, exact stratification becomes
impractical. The curse of dimensionality means that most covariate cells
have few or no observations. Regression adjustment models the
conditional outcome functions instead.

The main issue for estimating the $`\text{ATT}`$ is estimating the
untreated potential outcome, $`Y_i(0)`$.

To do so, we will try to estimate this conditional expectation function
``` math
  m_0(x) = E(Y_i(0) \mid X_i = x).
```
Since we only observe $`Y_i(0)`$ for untreated units, we will run a
regression of $`Y_i`$ on $`X_i`$ using the untreated units,
i.e. $`D_i = 0`$.

Then, we will estimate the treated unit’s untreated potential outcome
with $`\hat{Y}_i(0) = X_i \hat{\beta}_0`$ where the $`_0`$ is used to
emphasize this is the regression ran on the untreated group. In
practice, we will want our models to be more flexible, so we can include
transformations of the underling $`X_i`$ variables.

For a binary treatment, the ATT estimate can be formed as
``` math
 \widehat{\text{ATT}}^{RA}
  = \frac{1}{n_1} \sum_{i: D_i = 1} \{Y_i - \widehat m_0(X_i)\}.
```

``` r

ra(y ~ treat(d, ref = 0) + x1 + x2, data = example_data)
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
```

### Matching

Matching avoids a global functional-form assumption by finding
observations that are close in $`X_i`$ or in an estimated propensity
score. Let
``` math
  \mathcal J_M(i)
```
be the set of (M) “closest” observations from the opposite treatment
group (we will discuss how to define closest in a moment).

The missing potential outcome is imputed by their average:
``` math
  \hat{Y}_{i}(0) = \sum_{j \in \mathcal{J}_M(i) Y_j}
```
. The treatment effect is then formed as
``` math
 \widehat{\text{ATT}}^{Matching}
  = \frac{1}{n_1} \sum_{i: D_i = 1} \{Y_i - \sum_{j \in \mathcal{J}_M(i)}\}.
```

We can either use `nnmatch` to match directly on the full set of
covariates or use
[`psmatch()`](https://kylebutts.github.io/teffects-r/reference/psmatch.md)
matches on the estimated propensity score.

``` r

(ps_fit <- psmatch(y ~ treat(d, ref = 0) + x1 + x2, data = example_data, nneighbor = 1))
#> Treatment-effects estimation
#> Estimator: propensity-score matching
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  2.05530    0.17456  11.774 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
(nn_fit <- nnmatch(y ~ treat(d, ref = 0) + x1 + x2, data = example_data, nneighbor = 1))
#> Treatment-effects estimation
#> Estimator: nearest-neighbor matching
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]   2.0576     0.1578   13.04 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

There are a lot of options on how to do matching. The [Matching
details](https://kylebutts.github.io/teffects-r/articles/stata-nnmatch.md)
vignette documents options the nearest-neighbor rules used by
`teffects nnmatch`, including ties, exact matching, calipers, and
standard errors. Both match with replacement. More neighbors reduce
variance but can make the matches less similar. A caliper can rule out
comparisons that are too far apart. Exact matching in `ematch` is
applied before the distance calculation, and `biasadj` can correct the
remaining smooth covariate difference.

Matching makes each comparison explicit, but it does not make
incomparable units comparable. Check the distance between matches and
covariate balance.

In practice, matching is not perfect and we will have some differences
in the distribution of $`X`$ for treated and control units. It is
recommended to do a regression adjustment estimator between the treated
and the matched control group, using `biasadj` argument.

``` r

(nn_fit_bias_adj <- nnmatch(y ~ treat(d, ref = 0) + x1 + x2,
  data = example_data, nneighbor = 1,
  biasadj = ~ x1 + x2
))
#> Treatment-effects estimation
#> Estimator: nearest-neighbor matching
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  1.98668    0.15523  12.798 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Propensity-score weighting

Inverse-probability of treatment weighting then gives more weight to
observations whose treatment status is less expected given their
covariates. For the ATT, treated observations receive their usual weight
and controls receive IPTW weights:
``` math
  \frac{\pi(X_i)}{1 - \pi(X_i)}.
```
Our weighting procedure puts more weight on control units who have a
large propensity score because they look like they should be treated.

The treatment effect estimator is formed as:
``` math
 \widehat{\text{ATT}}^{IPTW}
  = \frac{1}{n_1} \sum_{i: D_i = 1} Y_i - \frac{1}{n_0} \sum_{j: D_j = 0} \frac{\pi(X_j)}{1 - \pi(X_j)} Y_j.
```

That is, we do a difference in means using a weighted control group.

There are different ways we could choose to estimate the propensity
score.
[`ipw()`](https://kylebutts.github.io/teffects-r/reference/ipw.md) uses
normalized weighted means. The options `psmethod = "logit"` and
`"probit"` estimate treatment probabilities by maximum likelihood.
`"ipt"` and `"cbps"` choose probabilities by solving equations that
target covariate balance.

``` r

(ipw_fit <- ipw(y ~ treat(d, ref = 0) + x1 + x2, data = example_data))
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
(ipt_fit <- ipw(y ~ treat(d, ref = 0) + x1 + x2, data = example_data, psmethod = "ipt"))
#> Treatment-effects estimation
#> Estimator: inverse-probability weighting
#> Treatment model: ipt 
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  1.95719    0.13382  14.626 < 2.2e-16 ***
#> POmean[0]    0.69266    0.11957   5.793 6.913e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
(ipt_fit <- ipw(y ~ treat(d, ref = 0) + x1 + x2, data = example_data, psmethod = "cbps"))
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
```

Extreme propensity scores create extreme weights (e.g. a propensity
score of 0.999 gives that control unit a weight of 999).  
Therefore, it is important to inspect the fitted scores and their
overlap before interpreting an IPTW estimate. The default behavior is to
stop when a score violates `pstolerance`. The opt-in
`psaction = "clamp"` changes the propensity scores to be in some range
(e.g. \[0.05, 0.95\]) and `psaction = "trim"` drops any units outside of
a given range.

### Doubly robust estimators

The final two estimators are called “doubly robust estimators”. Both use
an outcome model $`m_d(x)`$ and a treatment model $`\pi_d(x)`$.

First, we have the IPWRA estimator.
[`ipwra()`](https://kylebutts.github.io/teffects-r/reference/ipwra.md)
uses the treatment probabilities to weight the arm-specific outcome
regressions. After running the IPW-weighted least squares regression
using IPTW weights, regression adjustment proceeds as normal.

Second, we have the augmented IPW estimator.
[`aipw()`](https://kylebutts.github.io/teffects-r/reference/aipw.md)
combines outcome predictions with inverse-probability-weighted
residuals. The AIPW estimator for a potential-outcome mean is

``` math
  \widehat{\mu}_d^{AIPW} = \frac{1}{n}\sum_{i=1}^n \widehat{m}_d(X_i) + 
  \frac{1}{n}\sum_{i=1}^n  \frac{1(D_i=d)}{\widehat e_d(X_i)} \{Y_i-\ \widehat{m}_d(X_i)\}.
```

You can view it as a regression adjustment estimate for the average
potential outcome with a “correction” for modelling errors: the
estimator then rebalances $`Y_i-\ \widehat{m}_d(X_i)`$

Its double-robust property has a precise meaning: under the other
identifying and regularity conditions, it remains consistent if either
the outcome model or the treatment model is correctly specified. It does
not protect against an omitted confounder, failure of overlap, or
misspecification of both models.

``` r

(aipw_fit <- aipw(y ~ treat(d, ref = 0) + x1 + x2, data = example_data))
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
(ipwra_fit <- ipwra(y ~ treat(d, ref = 0) + x1 + x2, data = example_data))
#> Treatment-effects estimation
#> Estimator: inverse-probability-weighted regression adjustment
#> Propensity-score method: logit 
#> Outcome model: linear
#> Statistic: ATE 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> ATE[1 vs 0]  1.95314    0.13257 14.7331 < 2.2e-16 ***
#> POmean[0]    0.69230    0.11922  5.8069 6.366e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
