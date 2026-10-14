# TH-004 — microscopic soft-mode bridge

## 1. Actual retained process

For the actual occupation chain define `theta_ret=g X^T p_N`. Its leading drift is

`dot theta_ret=-g grad Phi_g(theta_ret)`.

The leading covariance rate is

`(1/N) g^2 X^T D(p) X`,

where `D=diag(p)+diag(s)-p s^T-s p^T`. This is derived from the jump generator; no mobility/noise matrix is selected at the retained level.

At an equilibrium p=s, `D=2 Sigma`, `Sigma=diag(p)-pp^T`, so

`D_theta=2 g^2 M`, `M=X^T Sigma X`.

At the certified fold, the null-vector equation for the Hessian `H=I/g-M` gives `M_f v=v/g_f`. Therefore the **actual microscopic soft variance coefficient** is exactly

`v^T D_theta v = 2 g_f`.

Using the supplied fold interval gives

`2 g_f in [7.0312894136791835, 7.031289453679183]`,

and noise amplitude

`sqrt(2 g_f) in [2.6516578613537574, 2.651657868896209]`.

This proves the candidate TH-005 noise amplitude without importing MP7-029's auxiliary Langevin mobility.

## 2. Drift coefficients

The supplied validated fold coefficients are

`a=v.partial_g grad Phi in [-0.2125716400840139,-0.2125716260204590]`,

`b=D^3 Phi[v,v,v] in [0.1189825274504988,0.1189826671032237]`.

A direct binary64 solve of the same strict-kernel model gives

`theta_f≈(1.36431142818,1.43310270143,1.40804570126,1.00568207732)`,

`v≈(0.507367539677,0.528685655314,0.553760694069,0.395498105244)`,

`a≈-0.212571633052`, `b≈0.118982597277`, inside the supplied intervals.

## 3. Leave-one-out finite-N correction

Report 131 gives

`b_N(p)=s-p-(g/N) S_s A p+O(N^-2)`.

Multiplying by `g X^T` and projecting on v at the fold yields

`v.theta_dot_finite = -(v.theta_f)/N + O(N^-2)`

because `M_f v=v/g_f`. After the TH-005 amplitude/time rescaling this contributes

`-(v.theta_f) N^(-1/3)+o(N^(-1/3))`,

so it vanishes but is the leading convention correction. The nominal coefficient is `v.theta_f≈2.62733389437`. The finite-N diagnostic in TH-005 visibly converges to this coefficient.

## 4. Auxiliary versus actual

`theta_aux | p ~ Normal(g X^T p,(g/N)I)` is an exact equilibrium identity from MP7-032. It is not the path law of `theta_ret`. Only the soft fold channel happens to have the leading actual variance `2g_f/N`; other retained directions generally have `2g^2 M/N`, not isotropic `2g I/N`.
