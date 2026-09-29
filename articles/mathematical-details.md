# Mathematical details

This vignette defines the smooth, mean-based estimators in `teffects`
and gives the influence function used for their robust standard errors.
The estimators are ordinary regression adjustment (RA), normalized
inverse- probability weighting (IPW), inverse-probability-weighted
regression adjustment (IPWRA), augmented IPW (AIPW), and the post-lasso
AIPW estimator in
[`telasso()`](https://kylebutts.github.io/teffects-r/reference/telasso.md).

Matching estimators are excluded because their asymptotic theory is
different. The final section covers
[`qte()`](https://kylebutts.github.io/teffects-r/reference/qte.md). It
estimates quantiles, so inversion of an estimated distribution function
replaces the smooth sample-mean calculation and the package uses
bootstrap inference. The experimental LATE estimators are not part of
the exported API and require an additional instrument, so they are
outside the $`(X_i,D_i,Y_i)`$ setup used here.

## Data, identification, and targets

Let the observations

``` math
O_i=(X_i,D_i,Y_i), \qquad i=1,\ldots,n,
```

be independent and identically distributed. Treatment takes values in
$`\mathcal K=\{0,1,\ldots,K-1\}`$, with 0 denoting the reference level,
and $`Y_i(k)`$ is the potential outcome under level $`k`$. Define

``` math
m_k(X_i)=E(Y_i\mid D_i=k,X_i), \qquad
\pi_k(X_i)=P(D_i=k\mid X_i).
```

We write $`1(\cdot)`$ for the indicator function.

For a binary treatment, the notation reduces to

``` math
\pi(X_i)=P(D_i=1\mid X_i), \qquad
\pi_1(X_i)=\pi(X_i), \qquad
\pi_0(X_i)=1-\pi(X_i).
```

The causal interpretation uses the usual conditions:

``` math
Y_i=Y_i(D_i), \qquad
Y_i(k)\mathrel{\perp\!\!\!\perp}D_i\mid X_i, \qquad
0<\pi_k(X_i)<1
```

almost surely for each relevant treatment level. Consistency and
conditional exchangeability identify

``` math
E\{Y_i(k)\mid X_i\}=m_k(X_i).
```

Use one notation for population, treated, and control targets. Let
$`t=\star`$ denote the full population and let $`t=g\in\mathcal K`$
denote the population receiving treatment $`g`$. Set

``` math
\begin{aligned}
a_{\star i}&=1, & r_\star(X_i)&=1, & q_\star&=1,\\
a_{gi}&=1(D_i=g), & r_g(X_i)&=\pi_g(X_i), & q_g&=P(D_i=g).
\end{aligned}
```

The potential-outcome mean in target population $`t`$ is

``` math
\theta_{k\mid t}
=
\frac{E\{a_{ti}Y_i(k)\}}{q_t}
=
\frac{E\{a_{ti}m_k(X_i)\}}{q_t}.
```

Thus $`\theta_{k\mid\star}=E\{Y_i(k)\}`$. For a focal level $`k\ne0`$,

``` math
\begin{aligned}
\operatorname{ATE}_{k,0}
  &=\theta_{k\mid\star}-\theta_{0\mid\star},\\
\operatorname{ATT}_{k,0}
  &=\theta_{k\mid k}-\theta_{0\mid k},\\
\operatorname{ATC}_{k,0}
  &=\theta_{k\mid 0}-\theta_{0\mid 0}.
\end{aligned}
```

The transport weight for outcome level $`k`$ and target $`t`$ is

``` math
w_{k\mid t}(O_i)
=
1(D_i=k)\frac{r_t(X_i)}{\pi_k(X_i)}.
```

For the full population this is the usual inverse-probability weight
$`1(D_i=k)/\pi_k(X_i)`$. For target group $`g`$, it is
$`1(D_i=k)\pi_g(X_i)/\pi_k(X_i)`$. Notice that $`w_{g\mid g}=1(D_i=g)`$
and, when the propensity score is correct,
$`E\{w_{k\mid t}(O_i)\}=q_t`$.

Throughout, $`\mathbb{E}_n[f]=n^{-1}\sum_i f(O_i)`$ denotes the sample
average. Hats denote estimates and all population quantities in an
influence function are evaluated at the probability limits of the fitted
working models. If a working model is misspecified, the formula
describes inference around the estimator’s probability limit. It
describes inference for the causal target whenever the estimator’s
stated consistency condition holds.

## Influence-function convention

For a scalar estimand $`\theta`$, `teffects` uses the normalization

``` math
\sqrt n(\widehat\theta-\theta)
=
\frac{1}{\sqrt n}\sum_{i=1}^n\phi_i+o_p(1),
\qquad E(\phi_i)=0.
```

For a vector of estimates, $`\phi_i`$ is a row vector. The package
stores the sample analogues as `fit$influence` and computes

``` math
\widehat{\operatorname{Var}}(\widehat\theta)
=
\frac{1}{n^2}\sum_{i=1}^n
\widehat\phi_i\widehat\phi_i^\top.
```

This is the variance of $`\widehat\theta`$, not the asymptotic variance
of $`\sqrt n\widehat\theta`$.

### Estimated nuisance parameters

We assume the outcome-model covariate vector $`X_i`$ includes an
intercept. The implemented outcome model is

``` math
m_k(X_i;\beta_k)=X_i^\top\beta_k.
```

Ordinary arm-specific least squares solves

``` math
\mathbb{E}_n\left[1(D_i=k)X_i
\{Y_i-X_i^\top\widehat\beta_k\}\right]=0.
```

Writing

``` math
A_k=E(1(D_i=k)X_iX_i^\top),
```

the coefficient influence function is

``` math
\phi_{\beta_k,i}
=
A_k^{-1}1(D_i=k)X_i\{Y_i-m_k(X_i)\}.
\tag{1}
```

Let $`\gamma`$ collect the parameters of the propensity-score equations,
including every arm-specific block when necessary. If

``` math
\mathbb{E}_n[s_{\pi}(O_i;\widehat\gamma)]=0, \qquad
J_\gamma
=
E\left\{\frac{\partial s_{\pi}(O_i;\gamma)}
{\partial\gamma^\top}\right\},
```

then

``` math
\phi_{\gamma,i}=-J_\gamma^{-1}s_{\pi}(O_i;\gamma).
\tag{2}
```

It is convenient to define the log-probability derivative

``` math
\ell_k(X_i)
=
\frac{\partial\log\pi_k(X_i;\gamma)}{\partial\gamma},
\qquad
\ell_\star(X_i)=0.
```

For a treatment-defined target, $`\ell_t=\ell_g`$. It follows that

``` math
\frac{\partial w_{k\mid t}(O_i;\gamma)}{\partial\gamma}
=
w_{k\mid t}(O_i;\gamma)
\{\ell_t(X_i)-\ell_k(X_i)\}.
\tag{3}
```

Equations (1)–(3) make the effect of nuisance estimation explicit.
Treating $`\widehat\pi_k`$ or \$\widehatm_k\$ as fixed generally omits
terms that are present in the package’s robust variance estimator.

## Regression adjustment

### Definition

Ordinary RA, implemented by
[`ra()`](https://kylebutts.github.io/teffects-r/reference/ra.md), fits

\$\$ \widehat\beta_k^{RA} = \arg\min_b\sum\_{i=1}^n
1(D_i=k)(Y_i-X_i^\top b)^2, \qquad
\widehatm_k(X_i)=X_i^\top\widehat\beta_k^{RA}. \$\$

It then standardizes predictions to the observed target population:

\$\$ \widehat\theta\_{k\mid t}^{RA} =
\frac{\mathbb{E}\_n\[a\_{ti}\widehatm_k(X_i)\]}
{\mathbb{E}\_n\[a\_{ti}\]}. \tag{4} \$\$

For $`t=\star`$, equation (4) averages over the full sample. For an ATT
or ATC, it averages over the observations in the corresponding treatment
group.

### Influence function

Define the target-population regressor mean

``` math
b_t=\frac{E(a_{ti}X_i)}{q_t}.
```

The RA influence function for $`\theta_{k\mid t}`$ is

``` math
\boxed{
\phi_{k\mid t,i}^{RA}
=
\frac{a_{ti}}{q_t}
\{m_k(X_i)-\theta_{k\mid t}\}
+b_t^\top\phi_{\beta_k,i}.
}
\tag{5}
```

The first term captures sampling variation in the covariate distribution
used for standardization. The second propagates estimation of the
arm-specific outcome regression. RA is consistent for the causal target
when the corresponding conditional outcome model is correctly specified.

## Propensity-score estimating equations

All propensity-score methods enter the influence functions below through
$`\phi_{\gamma,i}`$ in (2). We also write $`X_i`$ for the
treatment-model covariate vector. The exact score $`s_\pi`$ depends on
`psmethod`.

For binary maximum-likelihood logit or probit, write
$`\pi(X_i;\gamma)=G(X_i^\top\gamma)`$ and $`g_i=G'(X_i^\top\gamma)`$.
The score is

``` math
s_{\pi}(O_i;\gamma)
=
X_i
\frac{\{D_i-\pi(X_i;\gamma)\}g_i}
{\pi(X_i;\gamma)\{1-\pi(X_i;\gamma)\}}.
\tag{6}
```

For logit, this simplifies to $`X_i\{D_i-\pi(X_i;\gamma)\}`$. A
multivalued treatment uses the usual multinomial-logit score.

Inverse-probability tilting (IPT) uses a logistic tilt. For an ATE or
potential-outcome mean, it fits one arm-specific score from

``` math
\mathbb{E}_n\left[
X_i\left\{\frac{1(D_i=k)}{\pi_k(X_i;\gamma_k)}-1\right\}
\right]=0
\tag{7}
```

for every $`k`$. With an intercept, (7) makes the arm’s raw inverse-
probability weights sum to $`n`$. For a binary ATT, the same equation is
fit for the comparison arm $`k=0`$; for a binary ATC it is fit for
$`k=1`$. The arm-specific binary IPT fits used for ATEs need not add to
one, but each $`\widehat\pi_k(X_i)`$ enters the estimator and its
influence function in the same way as a generalized propensity score.

For binary covariate-balancing propensity scores (CBPS), the implemented
equations are

``` math
\begin{array}{ll}
\text{ATE or POM:}&
\mathbb{E}_n\left[X_i\dfrac{D_i-\pi(X_i)}
{\pi(X_i)\{1-\pi(X_i)\}}\right]=0,\\[10pt]
\text{ATT:}&
\mathbb{E}_n\left[X_i\dfrac{D_i-\pi(X_i)}{1-\pi(X_i)}\right]=0,\\[10pt]
\text{ATC:}&
\mathbb{E}_n\left[X_i\dfrac{D_i-\pi(X_i)}{\pi(X_i)}\right]=0.
\end{array}
\tag{8}
```

For a multivalued treatment, the CBPS implementation balances adjacent
arms:

``` math
\mathbb{E}_n\left[
X_i\left\{
\frac{1(D_i=k)}{\pi_k(X_i)}
-
\frac{1(D_i=k-1)}{\pi_{k-1}(X_i)}
\right\}
\right]=0,
\qquad k=1,\ldots,K-1.
\tag{9}
```

The Jacobian of the selected equations supplies $`J_\gamma`$ in (2).

## Normalized inverse-probability weighting

### Definition

For every outcome level and target population,
[`ipw()`](https://kylebutts.github.io/teffects-r/reference/ipw.md) uses
the normalized, or Hájek, estimator

``` math
\widehat\theta_{k\mid t}^{IPW}
=
\frac{\mathbb{E}_n[\widehat w_{k\mid t}(O_i)Y_i]}
{\mathbb{E}_n[\widehat w_{k\mid t}(O_i)]}.
\tag{10}
```

For an ATE, these are normalized inverse-probability-weighted arm means.
For an ATT or ATC, equation (10) uses normalized odds weights to
transport outcomes from arm $`k`$ to the target arm.

### Influence function

Let

``` math
H_{k\mid t}=E\{w_{k\mid t}(O_i)\}
```

and define the propensity loading

``` math
c_{k\mid t}^{IPW}
=
\frac{1}{H_{k\mid t}}
E\left[
w_{k\mid t}(O_i)
\{Y_i-\theta_{k\mid t}\}
\{\ell_t(X_i)-\ell_k(X_i)\}
\right].
```

The influence function is

``` math
\boxed{
\phi_{k\mid t,i}^{IPW}
=
\frac{w_{k\mid t}(O_i)}{H_{k\mid t}}
\{Y_i-\theta_{k\mid t}\}
+
\left(c_{k\mid t}^{IPW}\right)^\top\phi_{\gamma,i}.
}
\tag{11}
```

The first term is the influence function if the propensity score were
known. The second captures estimation of the propensity score and is
included by `teffects`. If the propensity model is correct,
$`H_{k\mid t}=q_t`$. If $`k=t=g`$, the weight is simply $`1(D_i=g)`$,
its derivative is zero, and (11) reduces to the influence function of
the observed mean in group $`g`$.

Normalized IPW is consistent for the causal target when the
propensity-score model is correctly specified.

## Inverse-probability-weighted regression adjustment

### Definition

[`ipwra()`](https://kylebutts.github.io/teffects-r/reference/ipwra.md)
implements IPWRA. For each outcome level and target, it solves the
weighted least-squares problem

``` math
\widehat\beta_{k\mid t}^{IPWRA}
=
\arg\min_b\sum_{i=1}^n
\widehat w_{k\mid t}(O_i)(Y_i-X_i^\top b)^2.
\tag{12}
```

The fitted regression is then standardized over the observed target
group:

``` math
\widehat\theta_{k\mid t}^{IPWRA}
=
\frac{\mathbb{E}_n[a_{ti}X_i^\top
\widehat\beta_{k\mid t}^{IPWRA}]}
{\mathbb{E}_n[a_{ti}]}.
\tag{13}
```

For a binary ATT, for example, the treated regression has unit weights
and the control regression has odds weights $`\pi(X_i)/\{1-\pi(X_i)\}`$.
With a multivalued treatment, the package fits a separate pair of
transported regressions for each focal ATT or ATC contrast.

### Influence function

For brevity, write $`m_k(X_i)=X_i^\top\beta_{k\mid t}`$ in this section
and define

``` math
\begin{aligned}
A_{k\mid t}^{w}
&=E\{w_{k\mid t}(O_i)X_iX_i^\top\},\\
C_{k\mid t}
&=E\left[
w_{k\mid t}(O_i)X_i
\{Y_i-m_k(X_i)\}
\{\ell_t(X_i)-\ell_k(X_i)\}^\top
\right].
\end{aligned}
```

The weighted-regression coefficient influence function is

``` math
\phi_{\beta_{k\mid t},i}^{w}
=
\left(A_{k\mid t}^{w}\right)^{-1}
\left[
w_{k\mid t}(O_i)X_i\{Y_i-m_k(X_i)\}
+C_{k\mid t}\phi_{\gamma,i}
\right].
\tag{14}
```

Consequently,

``` math
\boxed{
\phi_{k\mid t,i}^{IPWRA}
=
\frac{a_{ti}}{q_t}
\{m_k(X_i)-\theta_{k\mid t}\}
+b_t^\top\phi_{\beta_{k\mid t},i}^{w}.
}
\tag{15}
```

The $`C_{k\mid t}\phi_{\gamma,i}`$ term propagates estimation of the
propensity score through the weighted regression. It disappears if the
propensity score is treated as known or if the outcome regression is
correct, because then the relevant weighted residual expectation is
zero.

With an intercept and the required regularity conditions, IPWRA is
doubly robust: (13) is consistent when either the linear outcome model
or the propensity-score model is correctly specified.

## Augmented inverse-probability weighting

### Definition

[`aipw()`](https://kylebutts.github.io/teffects-r/reference/aipw.md)
first fits the unweighted arm-specific outcome regressions used by RA.
It then defines

\$\$ \widehat h\_{k\mid t,i} = a\_{ti}\widehatm_k(X_i) + \widehat
w\_{k\mid t}(O_i) \\Y_i-\widehatm_k(X_i)\\ \$\$

and estimates

``` math
\widehat\theta_{k\mid t}^{AIPW}
=
\frac{\mathbb{E}_n[\widehat h_{k\mid t,i}]}{\mathbb{E}_n[a_{ti}]}.
\tag{16}
```

For the full population, (16) is the familiar average of

\$\$ \widehatm_k(X_i) + \frac{1(D_i=k)}{\widehat\pi_k(X_i)}
\\Y_i-\widehatm_k(X_i)\\. \$\$

For a binary ATT, the counterfactual control mean is

\$\$ \widehat\theta\_{0\mid1}^{AIPW} = \frac{1}{\sum_iD_i} \sum\_{i=1}^n
\left\[ D_i\widehatm_0(X_i)
+(1-D_i)\frac{\widehat\pi(X_i)}{1-\widehat\pi(X_i)}
\\Y_i-\widehatm_0(X_i)\\ \right\]. \$\$

### Influence function

Define the two nuisance loadings

``` math
\begin{aligned}
c_{\beta,k\mid t}^{AIPW}
&=
\frac{1}{q_t}E[\{a_{ti}-w_{k\mid t}(O_i)\}X_i],\\
c_{\gamma,k\mid t}^{AIPW}
&=
\frac{1}{q_t}E\left[
w_{k\mid t}(O_i)\{Y_i-m_k(X_i)\}
\{\ell_t(X_i)-\ell_k(X_i)\}
\right].
\end{aligned}
```

The influence function implemented by the parametric AIPW estimator is

``` math
\boxed{
\begin{aligned}
\phi_{k\mid t,i}^{AIPW}
={}&
\frac{1}{q_t}
\left[
a_{ti}\{m_k(X_i)-\theta_{k\mid t}\}
+w_{k\mid t}(O_i)\{Y_i-m_k(X_i)\}
\right]\\
&+
\left(c_{\beta,k\mid t}^{AIPW}\right)^\top
\phi_{\beta_k,i}
+
\left(c_{\gamma,k\mid t}^{AIPW}\right)^\top
\phi_{\gamma,i}.
\end{aligned}
}
\tag{17}
```

If both nuisance models are correct, both loadings are zero. Equation
(17) then reduces to the efficient influence function

``` math
\boxed{
\phi_{k\mid t,i}^{eff}
=
\frac{1}{q_t}
\left[
a_{ti}\{m_k(X_i)-\theta_{k\mid t}\}
+w_{k\mid t}(O_i)\{Y_i-m_k(X_i)\}
\right].
}
\tag{18}
```

If exactly one working model is correct, (16) remains consistent but one
of the chain-rule corrections in (17) can remain nonzero. This is why
the parametric implementation uses the complete expression (17), rather
than always treating (18) as the influence function.

## Post-lasso and cross-fitted AIPW

[`telasso()`](https://kylebutts.github.io/teffects-r/reference/telasso.md)
is a binary-treatment implementation of the same AIPW functional in
(16). It uses plugin lasso to select outcome regressors separately in
each arm and to select propensity-score regressors, then refits the
selected models without a penalty.

With $`L`$ cross-fitting folds, let $`\ell(i)`$ be observation $`i`$’s
fold and let \$\widehatm_k^{(-\ell(i))}\$ and
$`\widehat\pi^{(-\ell(i))}`$ be trained without that fold. Substituting
these out-of-fold predictions into (16) defines the cross-fitted
estimate. If multiple random partitions are requested,
[`telasso()`](https://kylebutts.github.io/teffects-r/reference/telasso.md)
averages the estimates from the partitions.

Under overlap, suitable moment conditions, approximate sparsity, and
nuisance rates satisfying the usual product-rate condition, Neyman
orthogonality makes first-order nuisance-estimation terms disappear. Its
influence function is therefore (18). In binary notation, the
full-population versions are

``` math
\begin{aligned}
\phi_{1\mid\star,i}^{telasso}
&=
m_1(X_i)-\theta_{1\mid\star}
+\frac{D_i}{\pi(X_i)}\{Y_i-m_1(X_i)\},\\
\phi_{0\mid\star,i}^{telasso}
&=
m_0(X_i)-\theta_{0\mid\star}
+\frac{1-D_i}{1-\pi(X_i)}\{Y_i-m_0(X_i)\}.
\end{aligned}
\tag{19}
```

The stored influence values are the estimated orthogonal scores,
averaged over repeated partitions. No finite-dimensional post-selection
chain-rule term is added. With `xfolds = 1`, the same score is evaluated
using full-sample nuisance fits; its validity therefore requires
stronger conditions than the cross-fitted construction.

## Contrasts and reported coefficients

Every mean-based estimator first obtains influence functions for the
relevant potential-outcome means. Linear contrasts require no new
derivation. If

``` math
\widehat\tau_{k,0\mid t}
=
\widehat\theta_{k\mid t}-\widehat\theta_{0\mid t},
```

then, for any estimator $`E`$,

``` math
\boxed{
\phi_{\tau_{k,0\mid t},i}^{E}
=
\phi_{k\mid t,i}^{E}-\phi_{0\mid t,i}^{E}.
}
\tag{20}
```

The efficient influence functions for the binary treatment-effect
targets follow directly from (18). Write $`p_1=P(D_i=1)`$,
$`p_0=P(D_i=0)`$, and suppress the observation subscript on the
right-hand sides. Then

``` math
\begin{aligned}
\phi_{ATE}^{eff}
={}&
m_1(X)-m_0(X)-\operatorname{ATE}
+\frac{D}{\pi(X)}\{Y-m_1(X)\}
-\frac{1-D}{1-\pi(X)}\{Y-m_0(X)\},\\[4pt]
\phi_{ATT}^{eff}
={}&
\frac{1}{p_1}\left[
D\{Y-m_0(X)-\operatorname{ATT}\}
-(1-D)\frac{\pi(X)}{1-\pi(X)}\{Y-m_0(X)\}
\right],\\[4pt]
\phi_{ATC}^{eff}
={}&
\frac{1}{p_0}\left[
(1-D)\{m_1(X)-Y-\operatorname{ATC}\}
+D\frac{1-\pi(X)}{\pi(X)}\{Y-m_1(X)\}
\right].
\end{aligned}
\tag{21}
```

`stat = "pomeans"` reports the full-population $`\theta_{k\mid\star}`$.
`stat = "ate"` reports (20) with $`t=\star`$, `stat = "att"` uses
$`t=k`$, and `stat = "atc"` uses $`t=0`$. The additional reference
potential-outcome mean printed with a treatment-effect contrast uses the
same target population and its corresponding influence-function column.

The following table summarizes the causal consistency condition and the
influence function stored by the package.

| Estimator | Sufficient working-model condition | Stored influence function |
|----|----|----|
| RA | $`m_k(X)`$ correct | Equation (5) |
| normalized IPW | $`\pi_k(X)`$ correct | Equation (11) |
| IPWRA | $`m_k(X)`$ or $`\pi_k(X)`$ correct | Equation (15) |
| parametric AIPW | $`m_k(X)`$ or $`\pi_k(X)`$ correct | Equation (17) |
| [`telasso()`](https://kylebutts.github.io/teffects-r/reference/telasso.md) | orthogonal-score nuisance-rate conditions | Equation (18) |

All formulas assume nonsingular population Jacobians, finite second
moments, and a fixed estimation population. The default strict-overlap
behavior of the package fits this regular setting. Clamping introduces a
kink at the clamp boundary, while trimming changes the retained
population; inference after either operation should be interpreted with
those qualifications.

## References

Cattaneo, M. D. (2010). Efficient semiparametric estimation of
multivalued treatment effects under ignorability. *Journal of
Econometrics* 155(2), 138–154.

Chernozhukov, V., D. Chetverikov, M. Demirer, E. Duflo, C. Hansen, W.
Newey, and J. Robins (2018). Double/debiased machine learning for
treatment and structural parameters. *The Econometrics Journal* 21(1),
C1–C68.

Graham, B. S., C. C. D. X. Pinto, and D. Egel (2012). Inverse
probability tilting for moment condition models with missing data.
*Review of Economic Studies* 79(3), 1053–1079.

Imai, K., and M. Ratkovic (2014). Covariate balancing propensity score.
*Journal of the Royal Statistical Society: Series B* 76(1), 243–263.

Słoczyński, T., S. D. Uysal, and J. M. Wooldridge (2025). Covariate
balancing and the equivalence of weighting and doubly robust estimators
of average treatment effects. arXiv:2310.18563.
