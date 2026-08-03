Yes. For **hazard-scale CATE with risk-set-orthogonalized partial likelihood**, the most proper Semenova-style inference is not OLS on a scalar signal, but a **cross-fitted orthogonal Z-estimator**. The Semenova logic is the same: first build an orthogonal score using ML nuisances, then project the structural function onto (p(X)), and use a sandwich variance estimator. Semenova–Chernozhukov’s setup is exactly: construct an orthogonal signal/score with cross-fitting, project onto basis (p(X)), and use an empirical covariance formula of the form (Q^{-1}E[p p'(U+r_g)^2]Q^{-1}). 

Here the target is hazard-scale CATE / DINA:

[
\tau_0(x)=\log\lambda_1(t\mid x)-\log\lambda_0(t\mid x),
]

and under proportional hazards this does not depend on (t). Gao–Hastie define exactly this Cox DINA target. 

---

## A. Procedure

### Step 0: assume linear hazard-scale CATE

Choose low-dimensional basis

[
p(x)\in\mathbb R^d,
]

and assume

[
\boxed{
\tau_0(x)=p(x)^\top\beta_0.
}
]

Observed data:

[
O_i=(Y_i^c,\Delta_i,W_i,X_i),
]

where

[
Y_i^c=\min(T_i,C_i),\qquad \Delta_i=1(T_i\le C_i).
]

Define counting-process notation:

[
N_i(t)=1(Y_i^c\le t,\Delta_i=1),
\qquad
Y_i(t)=1(Y_i^c\ge t).
]

---

### Step 1: cross-fit nuisances

For each fold (k), estimate nuisances on the other folds:

[
\hat e(x)=P(W=1\mid X=x),
]

[
\hat\eta_0(x),\qquad \hat\eta_1(x),
]

and the risk-set survival probability

[
\hat R_w(t,x)=\widehat P(Y^c\ge t\mid W=w,X=x).
]

You can estimate

[
R_w(t,x)=S_w(t\mid x)G_w(t\mid x)
]

using survival models for event and censoring, or estimate (P(Y^c\ge t\mid W=w,X=x)) directly.

For each event time (t), construct the **risk-set modified propensity**

[
\boxed{
\hat a_t(x)
===========

\frac{
\hat e(x)\hat R_1(t,x)\exp{\hat\eta_1(x)}
}{
\hat e(x)\hat R_1(t,x)\exp{\hat\eta_1(x)}
+
{1-\hat e(x)}\hat R_0(t,x)\exp{\hat\eta_0(x)}
}.
}
]

Then define the risk-set baseline offset

[
\boxed{
\hat\nu_t(x)
============

\hat a_t(x)\hat\eta_1(x)
+
{1-\hat a_t(x)}\hat\eta_0(x).
}
]

This is the perturbation relative to Gao–Hastie’s Cox partial-likelihood algorithm: their original Cox method uses a time-invariant (a(x),\nu(x)), while this version uses (a_t(x),\nu_t(x)) inside each risk set. Their supplement shows the full-likelihood version has the (c_n^2+n^{-1/2}) robustness result, while their partial-likelihood result only holds generally under no treatment effect; the risk-set version is designed to restore orthogonality for nonzero effects. 

---

### Step 2: construct risk-set residualized treatment covariate

Define

[
\boxed{
q_i(t)
======

{W_i-\hat a_t(X_i)}p(X_i).
}
]

The Cox linear predictor is

[
\theta_i(t;\beta)
=================

\hat\nu_t(X_i)+q_i(t)^\top\beta.
]

Then define the risk-set average

[
\hat{\bar q}(t;\beta)
=====================

\frac{
\sum_{j=1}^n
Y_j(t)\exp{\theta_j(t;\beta)}q_j(t)
}{
\sum_{j=1}^n
Y_j(t)\exp{\theta_j(t;\beta)}
}.
]

---

### Step 3: solve the orthogonal partial-likelihood score

Estimate (\beta_0) by solving

[
\boxed{
0
=

\frac1n
\sum_{i=1}^n
\int
\left[
q_i(t)-\hat{\bar q}(t;\beta)
\right]dN_i(t).
}
]

Equivalently, maximize the perturbed Cox partial likelihood

[
\boxed{
\ell_n^\perp(\beta)
===================

\sum_{i:\Delta_i=1}
\left[
\hat\nu_{Y_i^c}(X_i)
+
{W_i-\hat a_{Y_i^c}(X_i)}p(X_i)^\top\beta
-----------------------------------------

\log
\sum_{j:Y_j^c\ge Y_i^c}
\exp
\left{
\hat\nu_{Y_i^c}(X_j)
+
{W_j-\hat a_{Y_i^c}(X_j)}p(X_j)^\top\beta
\right}
\right].
}
]

Then

[
\hat\tau(x)=p(x)^\top\hat\beta.
]

In R this is basically a time-dependent Cox model:

```r
coxph(
  Surv(start, stop, event) ~ q1 + q2 + ... + qd + offset(nu_t),
  data = expanded_crossfit_data,
  robust = TRUE
)
```

where

[
q_\ell(t)={W-\hat a_t(X)}p_\ell(X).
]

---

## B. CLT

Define the population score contribution

[
\psi_i(\beta_0,\eta_0)
======================

\int
\left[
q_i^0(t)-\bar q^0(t;\beta_0)
\right]dM_i(t),
]

where

[
q_i^0(t)={W_i-a_t^0(X_i)}p(X_i),
]

and

[
dM_i(t)
=======

## dN_i(t)

Y_i(t)\exp{\nu_t^0(X_i)+q_i^0(t)^\top\beta_0}
d\Lambda_0(t).
]

Let

[
Q
=

-\partial_\beta
E\left[
\int
{q_i^0(t)-\bar q^0(t;\beta)}dN_i(t)
\right]_{\beta=\beta_0}.
]

Equivalently,

[
\boxed{
Q
=

E
\int
V_q^0(t;\beta_0)dN_i(t),
}
]

where

[
V_q^0(t;\beta_0)
================

\frac{
E\left[
Y(t)e^{\theta^0(t)}
{q^0(t)-\bar q^0(t)}{q^0(t)-\bar q^0(t)}^\top
\right]
}{
E\left[
Y(t)e^{\theta^0(t)}
\right]
}.
]

Then under exact linearity

[
\tau_0(x)=p(x)^\top\beta_0,
]

cross-fitting, bounded basis, overlap, nondegenerate information, and nuisance rates fast enough so that the orthogonal remainder is (o_p(n^{-1/2})),

[
\boxed{
\sqrt n(\hat\beta-\beta_0)
==========================

Q^{-1}
\frac1{\sqrt n}
\sum_{i=1}^n
\psi_i(\beta_0,\eta_0)
+
o_p(1).
}
]

Therefore

[
\boxed{
\sqrt n(\hat\beta-\beta_0)
\Rightarrow
N(0,\Omega),
}
]

with

[
\boxed{
\Omega
======

Q^{-1}
E\left[
\psi_i(\beta_0,\eta_0)\psi_i(\beta_0,\eta_0)^\top
\right]
Q^{-1}.
}
]

For a point (x),

[
\boxed{
\frac{
\sqrt n{p(x)^\top\hat\beta-p(x)^\top\beta_0}
}{
\sqrt{p(x)^\top\Omega p(x)}
}
\Rightarrow N(0,1).
}
]

This is the Cox partial-likelihood analogue of Semenova’s

[
\Omega=Q^{-1}E[p(X)p(X)'(U+r_g(X))^2]Q^{-1}.
]

Here the scalar residual

[
Y_i(\hat\eta)-p(X_i)^\top\hat\beta
]

is replaced by the **orthogonal Cox score residual**

[
\hat\psi_i.
]

---

## C. Variance estimator

Define

[
\hat\theta_i(t)=\hat\nu_t(X_i)+\hat q_i(t)^\top\hat\beta,
]

[
\hat q_i(t)={W_i-\hat a_t(X_i)}p(X_i).
]

The empirical risk-set average is

[
\hat{\bar q}(t)
===============

\frac{
\sum_{j=1}^n
Y_j(t)e^{\hat\theta_j(t)}\hat q_j(t)
}{
\sum_{j=1}^n
Y_j(t)e^{\hat\theta_j(t)}
}.
]

The empirical risk-set variance is

[
\hat V_q(t)
===========

\frac{
\sum_{j=1}^n
Y_j(t)e^{\hat\theta_j(t)}
{\hat q_j(t)-\hat{\bar q}(t)}{\hat q_j(t)-\hat{\bar q}(t)}^\top
}{
\sum_{j=1}^n
Y_j(t)e^{\hat\theta_j(t)}
}.
]

Then

[
\boxed{
\hat Q
======

\frac1n
\sum_{i=1}^n
\int
\hat V_q(t)dN_i(t)
==================

\frac1n
\sum_{i:\Delta_i=1}
\hat V_q(Y_i^c).
}
]

Next define the Breslow baseline hazard increment

[
d\hat\Lambda_0(t)
=================

\frac{
dN_\cdot(t)
}{
\sum_{j=1}^n Y_j(t)e^{\hat\theta_j(t)}
}.
]

Then the estimated martingale increment is

[
d\hat M_i(t)
============

## dN_i(t)

Y_i(t)e^{\hat\theta_i(t)}d\hat\Lambda_0(t).
]

Define the individual orthogonal score residual

[
\boxed{
\hat\psi_i
==========

\int
{\hat q_i(t)-\hat{\bar q}(t)}d\hat M_i(t).
}
]

Then

[
\boxed{
\hat\Sigma
==========

\frac1n
\sum_{i=1}^n
\hat\psi_i\hat\psi_i^\top.
}
]

Finally,

[
\boxed{
\hat\Omega
==========

\hat Q^{-1}\hat\Sigma \hat Q^{-1}.
}
]

Pointwise standard error:

[
\boxed{
\widehat{\mathrm{se}}{\hat\tau(x)}
==================================

\sqrt{
\frac1n
p(x)^\top\hat\Omega p(x)
}.
}
]

Confidence interval:

[
\boxed{
p(x)^\top\hat\beta
\pm
1.96
\sqrt{
\frac1n p(x)^\top\hat\Omega p(x)
}.
}
]

---

## D. Connection to your formula

Your Semenova formula is

[
\Omega
======

Q^{-1}
E\left[
p(X)p(X)^\top
{U+r_g(X)}^2
\right]
Q^{-1}.
]

For the risk-set-orthogonalized Cox PL estimator, replace:

[
p(X){U+r_g(X)}
]

by

[
\psi_i
======

\int
{q_i(t)-\bar q(t)}dM_i(t).
]

So the direct analogue is

[
\boxed{
\Omega
======

Q^{-1}
E[\psi_i\psi_i^\top]
Q^{-1}.
}
]

and

[
\boxed{
\hat\Omega
==========

\hat Q^{-1}
\left[
\frac1n\sum_i\hat\psi_i\hat\psi_i^\top
\right]
\hat Q^{-1}.
}
]

If the linear model is exact,

[
r_g(x)=0.
]

If instead (\tau_0(x)) is only approximated by (p(x)^\top\beta_0), then (\beta_0) should be interpreted as the **best linear projection / score projection** target, and the approximation error enters the influence function just like Semenova’s (r_g(X)).
