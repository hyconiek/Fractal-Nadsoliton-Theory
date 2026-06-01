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

The active `P1866` symbolic lane extracts strict moments at \(d_{\text{ref}}=1\):

\[
c_0=K_{\text{strict}}(1),\quad
c_1=\left.\frac{dK_{\text{strict}}}{dd}\right|_{d=1},\quad
c_2=\frac12\left.\frac{d^2K_{\text{strict}}}{dd^2}\right|_{d=1}.
\]

The current effective-coupling map is:

```text
m2_eff = m2*(1 + c0)
lam_eff = lam*(1 + c1**2)
y_eff = y*(1 + c0/2)
g_eff = g*(1 + c2)
xi_eff = xi*(1 + c0)
```

`P2363/S1313` adds bridge-completed moment transport: the legacy bridge kernel may feed this map only after the explicit APD completion

\[
Q_{\mathrm{APD}}(d)=\frac1{\alpha_{\mathrm{geo}}}
\frac{\cos(\omega d+\phi)}{\cos(\omega_L d+\phi_L)}
\frac{1+\beta_{\mathrm{tors}}d}{1+\beta d^\eta}
\]

is multiplied into \(K_{\text{legacy\_ont}}\).  Then \(K_{\text{legacy\_ont}}Q_{\mathrm{APD}}=K_{\text{strict}}\), so the completed moments match \((c_0,c_1,c_2)\) exactly.  Raw legacy moments, including amplitude-normalized raw legacy moments, are not admissible substitutes for the strict \(L_{\text{total}}\) coefficient input.

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

## 6) P2362/S1312 explicit EOM/Lagrangian supplement

EOM/Lagrangian track is selector-independent: selector closure is a parallel problem, not a prerequisite for writing, varying, or auditing the current strict-lane \(L_{\text{total}}\).  The selector/QW-2191 track controls branch/source actualization claims; it is not an input needed to continue the termwise Euler-Lagrange export below.

### 6.1 Covariant sector EOM rows

Current sector rows inherited from the strict full-lagrangian export:

\[
\begin{aligned}
E_{\phi}:&\quad
\Box\phi + m_\phi^2\phi + \frac{\lambda_3}{2}\phi^2
+ \frac{\lambda_4}{3!}\phi^3
+ 2\lambda_{\phi H}\phi(H^\dagger H)
+ 2\xi_{\phi R}R\phi = 0,\\
E_H:&\quad
D_\mu D^\mu H + \mu_H^2H + 2\lambda_H(H^\dagger H)H
+ \xi_{HR}RH + \lambda_{\phi H}\phi^2H
+ J_{\text{Yukawa}} = 0,\\
E_A^\mu:&\quad
\nabla_\nu(Z_a F_a^{\nu\mu}) + J_a^\mu
+ \chi_{RG}\nabla_\nu(RF_a^{\nu\mu}) = 0,\\
E_\psi:&\quad
i\gamma^a e_a{}^\mu D_\mu\psi_f - y_f H\cdot\psi_f = 0,\\
E_g{}_{\mu\nu}:&\quad
\frac{M_{\text{Pl}}^2}{2}G_{\mu\nu}+\Lambda g_{\mu\nu}
+H_{\mu\nu}^{(R^2,\mathrm{Ric}^2,\mathrm{Riem}^2)}
=T_{\mu\nu}^{\mathrm{SM+mix}}.
\end{aligned}
\]

These are sector-level equations, not a full tensor-resolved theorem-grade closure.

### 6.2 Reduced termwise computational EOM

Reduced termwise computational EOM from the current \(P2086/P2087\) witness use the fields \(\psi(x),A(x),h(x)\):

\[
\begin{aligned}
0=E_\psi&=-k_\psi\psi''-m_\psi\psi-4\lambda_4\psi^3-g_{\mathrm{mix}}Ah-2\zeta A^2\psi,\\
0=E_A&=-k_AA''-m_AA-4g_AA^3-g_{\mathrm{mix}}h\psi-2\zeta A\psi^2,\\
0=E_h&=-k_hh''-m_hh-4g_hh^3-g_{\mathrm{mix}}A\psi.
\end{aligned}
\]

Equivalently, solving for the second derivatives:

\[
\begin{aligned}
\psi''&=-\frac{m_\psi\psi+4\lambda_4\psi^3+g_{\mathrm{mix}}Ah+2\zeta A^2\psi}{k_\psi},\\
A''&=-\frac{m_AA+4g_AA^3+g_{\mathrm{mix}}h\psi+2\zeta A\psi^2}{k_A},\\
h''&=-\frac{m_hh+4g_hh^3+g_{\mathrm{mix}}A\psi}{k_h}.
\end{aligned}
\]

The termwise incidence matrix over fields \((\psi,A,h)\) has full field rank \(3\), and the recomposition residuals for the three reduced EOM rows are zero in the current computational audit.  Reduced termwise computational EOM are therefore an admissible parallel continuation target even while selector closure remains open.

### 6.3 P2364/S1314 bridge-completed FRW residual table

`P2364/S1314` adds a bridge-completed FRW residual table for scalar/gauge/gravity residual rows.  It uses the `P2363` APD-completed moment source, not raw legacy moments, and instantiates a named spatially flat constant-\(H\) FRW minisuperspace slice:

\[
R_{\mathrm{FRW}}=12H^2.
\]

The bridge-completed scalar row is derived from the effective minisuperspace density and exported as:

\[
E_{\phi}^{\mathrm{FRW}}=
\ddot\phi+3H\dot\phi
m_{\mathrm{eff}}^2\phi+\lambda_{\mathrm{eff}}\phi^3
-g_{\mathrm{eff}}A^2\phi
-2\xi_{\mathrm{eff}}R_{\mathrm{FRW}}\phi.
\]

The gauge mode row is:

\[
E_A^{\mathrm{FRW}}=\ddot A+3H\dot A+g_{\mathrm{eff}}\phi^2A.
\]

The gravity residual rows use the inherited FRW Einstein-perfect-fluid component class:

\[
R_{g,00}=3H^2+\Lambda-\kappa^2\rho_{\mathrm{total}},
\quad
R_{g,ii}=-3H^2+\Lambda-\kappa^2p_{\mathrm{total}}.
\]

The exported normal forms substitute back with zero residual:

\[
\ddot\phi=-3H\dot\phi-m_{\mathrm{eff}}^2\phi-\lambda_{\mathrm{eff}}\phi^3
+g_{\mathrm{eff}}A^2\phi+2\xi_{\mathrm{eff}}R_{\mathrm{FRW}}\phi,
\]

\[
\ddot A=-3H\dot A-g_{\mathrm{eff}}\phi^2A,
\]

\[
H^2=\frac{\kappa^2(\rho_{\mathrm{total}}-p_{\mathrm{total}})}6,
\quad
\Lambda=\frac{\kappa^2(\rho_{\mathrm{total}}+p_{\mathrm{total}})}2.
\]

This advances the EOM/Lagrangian residual-table track, but it is still a named-background witness.  Full tensor-resolved nonminimal metric variation, renormalized stress-energy closure, and atlas/background-independence remain open; selector/QW-2191 remains parallel.

### 6.4 P2365/S1315 nonminimal FRW metric variation lift

`P2365/S1315` promotes the FRW gravity rows by adding the explicit nonminimal FRW metric variation of the `P1866/P2363` term

```text
xi_eff*Phi^2*R
```

with

\[
F=\xi_{\mathrm{eff}}\Phi^2.
\]

The reduced tensor variation used on the named FRW family is:

\[
\Delta_{\mu\nu}^{\mathrm{nm}}
=2\kappa^2\left(FG_{\mu\nu}
+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)F\right).
\]

For homogeneous \(\Phi(t)\) and constant \(H\):

\[
\dot F=2\xi_{\mathrm{eff}}\Phi\dot\Phi,\quad
\ddot F=2\xi_{\mathrm{eff}}(\dot\Phi^2+\Phi\ddot\Phi),
\]

\[
\Delta_{00}^{\mathrm{nm}}
=2\kappa^2(3FH^2+3H\dot F),
\]

\[
\Delta_{ii}^{\mathrm{nm}}
=2\kappa^2(-3FH^2-\ddot F-2H\dot F).
\]

Thus the lifted metric residual rows are:

\[
R_{g,00}^{\mathrm{nm}}
=3H^2+\Lambda-\kappa^2\rho_{\mathrm{matter}}+\Delta_{00}^{\mathrm{nm}},
\]

\[
R_{g,ii}^{\mathrm{nm}}
=-3H^2+\Lambda-\kappa^2p_{\mathrm{matter}}+\Delta_{ii}^{\mathrm{nm}}.
\]

`P2365` solves these two rows algebraically for \(H^2\) and \(\Lambda\), substitutes the solution back with zero residual, and verifies the limit \(\xi_{\mathrm{eff}}\to0\) returns the `P2364` minimal FRW metric normal form.  This is still a named-FRW-family lift only: Bianchi-I replay, off-FRW tensor component tables, renormalized stress-energy closure, and atlas/background-independence remain open; selector/QW-2191 remains parallel.

## 7) Parallel selector candidate audit

`P2366/S1316` keeps the selector track separate from the EOM/Lagrangian track and audits the current best selector candidate.  The grep result is:

```text
chi11_selector_source = unique top logical bottleneck / only singleton selector unlock
```

The strongest concrete phase-origin selector candidate is:

```text
chiral bispectrum -> orientation
calibrated coprime Fourier phase -> source
```

On the audited \(Z_{12}\) finite model, the recovery formula is:

\[
s=k^{-1}\,12\,(\varphi^{\mathrm{ref}}_k(\mathrm{orientation})-\varphi^{\mathrm{obs}}_k)
\pmod{12},
\]

for \(k\in\{1,5,7,11\}\).  `P2366` verifies all 24 source/orientation rows and checks negative controls: translation-invariant magnitudes cannot select source, non-coprime modes alias sources, and without the chiral marker orientation remains two-valued.

This is still premise-based.  The phase origin and handedness are not derived from strict-core data; `beta_tors -> chi11` is not proven; QW-2191 remains open.

---

This sheet is intentionally operational and strict-lane only. It does **not** claim theorem-grade global closure, selector closure (`QW-2191`), nor full ToE closure.
