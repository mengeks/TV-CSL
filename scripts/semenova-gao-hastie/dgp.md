



Yes — for **first verification**, make the DGP as simple as possible: no ML nuisance estimation, no nonlinear \(X\), known nuisance first, then gradually add complexity. The goal is to check that your **risk-set-orthogonalized PL estimator + sandwich SE** gives unbiased \(\hat\beta\), correct SE, and 95% coverage.

## Simple DGP 1: clean Cox PH, linear hazard CATE

Let

\[
X=(X_1,X_2),\qquad X_1,X_2\overset{iid}{\sim}N(0,1).
\]

Treatment:

\[
W\mid X\sim \mathrm{Bernoulli}\{e(X)\},
\]

with

\[
e(X)=\operatorname{expit}(0.3X_1-0.3X_2).
\]

Survival model:

\[
\lambda(t\mid W,X)
=
\lambda_0
\exp\left\{
\eta_0(X)+W\tau(X)
\right\},
\]

where

\[
\lambda_0=0.1,
\]

\[
\eta_0(X)=0.5X_2,
\]

and the true hazard-scale CATE is linear:

\[
\boxed{
\tau(X)=\beta_0+\beta_1X_1.
}
\]

Use

\[
\beta_0=-0.5,\qquad \beta_1=-0.5.
\]

So treatment is beneficial, and benefit is stronger for larger \(X_1\). Gao–Hastie define the survival DINA target as the log hazard ratio \(\log\lambda_1(y\mid x)-\log\lambda_0(y\mid x)\), and under PH it no longer depends on time, so this DGP exactly matches their Cox estimand. fileciteturn2file1

---

## Censoring DGP

Use independent exponential censoring:

\[
C\mid W,X\sim \mathrm{Exp}\{\lambda_C(W,X)\},
\]

with

\[
\lambda_C(W,X)=0.05\exp(0.2W+0.2X_2).
\]

Observed data:

\[
Y^c=\min(T,C),\qquad \Delta=1(T\le C).
\]

This gives moderate censoring. Later you can tune \(0.05\) to get 20%, 40%, 60% censoring.

---

## How to simulate \(T\)

Because baseline hazard is constant,

\[
\Lambda(t\mid W,X)
=
0.1t\exp\{\eta_0(X)+W\tau(X)\}.
\]

If \(E\sim \mathrm{Exp}(1)\), then

\[
T=
\frac{E}{
0.1\exp\{\eta_0(X)+W\tau(X)\}
}.
\]

Likewise,

\[
C=
\frac{E_C}{
0.05\exp(0.2W+0.2X_2)
}.
\]

This is easy and avoids numerical inversion.

---

## Estimation target

Use basis

\[
p(X)=(1,X_1)^\top.
\]

Then the true parameter is

\[
\boxed{
\beta=(\beta_0,\beta_1)^\top=(-0.5,-0.5)^\top.
}
\]

You are not estimating RMST-CATE here. This first DGP is purely for verifying **hazard-scale CATE inference**.

---

## First verification version: oracle nuisances

For the first test, use true nuisances:

\[
e(X)=\operatorname{expit}(0.3X_1-0.3X_2),
\]

\[
\eta_0(X)=0.5X_2,
\]

\[
\eta_1(X)=\eta_0(X)+\tau(X).
\]

For each event time \(t\), compute

\[
R_w(t,X)=P(Y^c\ge t\mid W=w,X)
=
S_w(t\mid X)G_w(t\mid X).
\]

Here

\[
S_w(t\mid X)
=
\exp\left[
-0.1t\exp\{\eta_0(X)+w\tau(X)\}
\right],
\]

and

\[
G_w(t\mid X)
=
\exp\left[
-0.05t\exp(0.2w+0.2X_2)
\right].
\]

Then define your risk-set modified propensity:

\[
\boxed{
a_t(X)
=
\frac{
e(X)R_1(t,X)\exp\{\eta_1(X)\}
}{
e(X)R_1(t,X)\exp\{\eta_1(X)\}
+
\{1-e(X)\}R_0(t,X)\exp\{\eta_0(X)\}
}.
}
\]

And

\[
\nu_t(X)=a_t(X)\eta_1(X)+\{1-a_t(X)\}\eta_0(X).
\]

This is the “slight perturbation” of Gao–Hastie PL: their Cox nuisance uses \(a(X),\nu(X)\), while yours uses risk-set-specific \(a_t(X),\nu_t(X)\). Their paper’s Cox section defines \(a(X)\) using not-censored probabilities and notes PL robustness only under no treatment effect; the supplement shows the full-likelihood version has the clean \(c_n^2+n^{-1/2}\) rate, while PL generally loses it for nonzero treatment effects. fileciteturn2file1 fileciteturn2file2

---

## Estimator to implement

At each event time \(t_i=Y_i^c\) with \(\Delta_i=1\), define

\[
q_j(t_i)=\{W_j-a_{t_i}(X_j)\}p(X_j).
\]

Fit Cox partial likelihood with linear predictor

\[
\theta_j(t_i;\beta)
=
\nu_{t_i}(X_j)+q_j(t_i)^\top\beta.
\]

The score is

\[
\boxed{
0=
\sum_{i:\Delta_i=1}
\left[
q_i(Y_i^c)
-
\frac{
\sum_{j:Y_j^c\ge Y_i^c}
\exp\{\theta_j(Y_i^c;\beta)\}q_j(Y_i^c)
}{
\sum_{j:Y_j^c\ge Y_i^c}
\exp\{\theta_j(Y_i^c;\beta)\}
}
\right].
}
\]

In R, easiest implementation is expanded counting-process format:

```r
coxph(Surv(start, stop, event) ~ q0 + q1 + offset(nu_t), data = expanded_data)
```

where

\[
q0(t)=W-a_t(X),
\qquad
q1(t)=\{W-a_t(X)\}X_1.
\]

---

## What to check

For each replication, compute:

\[
\hat\beta=(\hat\beta_0,\hat\beta_1).
\]

Then across replications check:

\[
\mathrm{Bias}(\hat\beta_j)=E(\hat\beta_j-\beta_j),
\]

\[
\mathrm{SD}(\hat\beta_j),
\]

\[
\overline{\mathrm{SE}}(\hat\beta_j),
\]

\[
\mathrm{Coverage}
=
P\left(
\beta_j\in
\hat\beta_j\pm1.96\widehat{\mathrm{SE}}(\hat\beta_j)
\right).
\]

Expected sanity result:

| Quantity | Expected |
|---|---|
| bias | near 0 |
| empirical SD | close to mean SE |
| coverage | around 0.95 |
| normal QQ plot | roughly straight |

---

## Sample sizes

Start with:

\[
n=1000,\quad 2000,\quad 5000.
\]

Replications:

\[
B=200
\]

for first run. Then increase to \(B=500\) once code works.

You should see bias shrink roughly like \(1/\sqrt n\), and coverage stabilize near 95%.

---

## Second verification: estimated nuisances

After oracle works, estimate nuisances by simple parametric models:

\[
\hat e(X): \text{logistic regression } W\sim X_1+X_2.
\]

\[
\hat\eta_0(X),\hat\eta_1(X): \text{separate Cox models in } W=0,W=1.
\]

\[
\hat R_w(t,X)=\hat S_w(t\mid X)\hat G_w(t\mid X).
\]

Estimate censoring survival using Cox models for censoring:

\[
1-\Delta
\]

as censoring event, separately by \(W\), or pooled with interactions.

Use 2-fold or 5-fold cross-fitting. Semenova–Chernozhukov’s procedure is exactly cross-fit nuisance estimation, construct an orthogonal signal/score, then project onto low-dimensional basis \(p(X)\), with variance estimated by an empirical sandwich analogue. fileciteturn2file0

---

## Third verification: compare to original Gao–Hastie PL

Run three estimators:

1. **Oracle risk-set orthogonal PL**  
   Uses true \(a_t,\nu_t\). This should have best behavior.

2. **Estimated risk-set orthogonal PL**  
   Uses cross-fitted \(\hat a_t,\hat\nu_t\). This checks your actual method.

3. **Original fixed-\(a(X)\) PL**  
   Uses Gao–Hastie style \(a(X),\nu(X)\). Under nonzero \(\beta\), you expect more nuisance sensitivity, especially if censoring depends on \(X,W\).

This directly tests your theoretical claim.

---

## Minimal DGP summary

Use this as the experiment specification:

\[
X_1,X_2\sim N(0,1),
\]

\[
W\mid X\sim \mathrm{Bernoulli}\{\operatorname{expit}(0.3X_1-0.3X_2)\},
\]

\[
T\mid W,X\sim \mathrm{Exp}\left(
0.1\exp\{0.5X_2+W(-0.5-0.5X_1)\}
\right),
\]

\[
C\mid W,X\sim \mathrm{Exp}\left(
0.05\exp\{0.2W+0.2X_2\}
\right),
\]

\[
Y^c=\min(T,C),\qquad \Delta=1(T\le C).
\]

Target:

\[
\boxed{
\tau(X)=-0.5-0.5X_1=p(X)^\top\beta,\qquad p(X)=(1,X_1),\quad \beta=(-0.5,-0.5).
}
\]

This is the cleanest DGP for checking your **estimation and inference procedure** before moving to nonlinear \(X\), ML nuisances, or RMST evaluation.