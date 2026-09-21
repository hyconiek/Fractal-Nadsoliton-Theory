# FR5 — global large-`J6` tail

## Statement

In the physical four-amplitude model, for all `J3,J4,J5,J6>=0` satisfying

`y = exp(-2 J6) <= 1/5,000,000`,

one has

`lambda2(M4) <= sigma_*`.

This is a four-amplitude curvature result only; it is not a transfer to the full seven-coordinate Hessian.

## 1. Shifted boundary reserve

On the exact even-parity boundary `y=0`, rerun the R7P-048--055 Bernstein machinery at the stricter threshold

`theta = sigma_* - 10^-6`.

The shifted `z^2` coefficient stays strictly positive because the certified R7P-044 trace margin exceeds `3*10^-6`.  Hence the same characteristic criterion applies: a box is safe when the cleared shifted characteristic numerator has `A<=0` or its shifted derivative has `B>=0`.

The dyadic cover has 508 terminal leaves and no unresolved leaves: 356 `SAFE_A_NONPOS`, 150 `SAFE_B_NONNEG`, and 2 leaves inside

`r in [3293/8192,3294/8192], s,t in [8191/8192,1]`.

That dyadic box lies strictly inside the projection of the already-certified FR13 box

`|r-r_*|<=1/6200, 1-s<=1/6800, 1-t<=1/6800`.

Therefore every boundary point outside the FR13 projection satisfies

`lambda2(C_+) <= sigma_* - 10^-6`.

## 2. Exact physical odd-mass bound

Write `E=Z_even` and `O=Z_odd` for the conditional partition sums with the common `J3,J4,J5` fields.  R7P-057 certifies `q0=E/(E+O)>=1/2` at `J6=0`, so `O<=E`.

With `y=exp(-2J6)`, the finite-`J6` odd mass is therefore

`e = y O/(E+y O) <= y/(1+y)`.

At `y<=1/5,000,000`, this is `e<=1/5,000,001`.

## 3. Covariance perturbation

Let `C_+`, `C_-` be the parity-conditional covariances and `m_+`, `m_-` their means.  The exact mixture identity is

`M4=(1-e)C_+ + e C_- + e(1-e)(m_+-m_-)(m_+-m_-)^T`.

For every unit vector `z`, dropping the negative term `-e Var_+(zX)` gives

`z^T(M4-C_+)z <= e[Var_-(zX)+(m_--m_+)^2]`

`= e E_-[(zX-m_+)^2] <= e D^2`,

where `D` is the diameter of the finite C4 feature set.  Exact spectral-interval enumeration over the six cyclic separations certifies `D^2<5`.  Hence

`M4 <= C_+ + 5 e I`.

For the declared tail,

`5e <= 5/5,000,001 < 10^-6`.

Thus every point outside the local equality box remains below `sigma_*` after finite-`J6` mixing.

## 4. Local box

Inside the quarantined boundary box, `e<=1/5,000,001 < 1/400`, so the point lies inside FR13 including its parity coordinate.  FR13 directly gives `lambda2(M4)<=sigma_*` there.

Combining the exterior shifted reserve and the FR13 local leaf proves the global large-`J6` tail.
