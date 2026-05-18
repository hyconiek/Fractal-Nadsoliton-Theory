# Strict Equation Sheet (Kernel -> Full Current L_total)

Status: working strict-lane equation sheet (no false-pass, no legacy-role transfer).

## 1) Strict kernel with operational parameters

Operational strict kernel used in the current strict lane:

\[
K_{\text{strict}}(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
\]

Current working tuple (QW-2049 lineage):

\[
\omega=0.18575,\quad \phi=0.16250,\quad \beta=1.0,\quad \eta=1.8.
\]

## 2) Effective couplings (strict symbolic chain)

From strict coefficient map (B1 symbolic lane):

\[
m^2_{\text{eff}} = c_0 + c_1 + c_2,
\quad
\lambda_{\text{eff}} = c_0c_1 + c_2,
\quad
y_{\text{eff}} = c_1-c_2,
\quad
\xi_{\text{eff}} = \frac{c_2}{1+c_0^2}.
\]

## 3) Full current non-skeleton Lagrangian decomposition


a) Scalar sector
\[
\mathcal L_{\text{scalar}}=
\frac12(\partial_\mu\phi)(\partial^\mu\phi)
-\frac12 m_{\text{eff}}^2\phi^2
-\frac{\lambda_{\text{eff}}}{4}\phi^4.
\]

b) Fermion sector
\[
\mathcal L_{\text{fermion}}=\bar\psi\,(i\gamma^\mu\partial_\mu-y_{\text{eff}}\phi)\,\psi.
\]

c) Gauge sector
\[
\mathcal L_{\text{gauge}}=-\frac14 F_{\mu\nu}F^{\mu\nu}.
\]

d) Gravity + nonminimal coupling density
\[
\mathcal L_{\text{gravity density}}=\sqrt{-g}\left[\frac{R-2\Lambda}{2\kappa^2}+\xi_{\text{eff}}\phi^2R\right].
\]

## 4) Current compact full statement

\[
\mathcal L_{\text{total}}=\mathcal L_{\text{scalar}}+\mathcal L_{\text{fermion}}+\mathcal L_{\text{gauge}}+\mathcal L_{\text{gravity density}}.
\]

Equivalent project statement in current repo chain:
\[
L_{\text{total}} = L_{\text{scalar}} + L_{\text{fermion}} + L_{\text{gauge}} + \sqrt{-g}\left[\frac{R-2\Lambda}{2\kappa^2}+\xi_{\text{eff}}\phi^2R\right].
\]

## 5) EOM generation rule (current gate discipline)

\[
E_X := \partial_\mu\left(\frac{\partial L_{\text{total}}}{\partial(\partial_\mu X)}\right)-\frac{\partial L_{\text{total}}}{\partial X}
\]
for each field block \(X\in\{\phi,\psi,A_\mu,g_{\mu\nu},\dots\}\).

---

This sheet is intentionally operational and strict-lane only. It does **not** claim theorem-grade global closure, selector closure (`QW-2191`), nor full ToE closure.
