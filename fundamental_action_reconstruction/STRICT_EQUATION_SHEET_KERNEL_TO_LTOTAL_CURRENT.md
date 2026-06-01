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

## P2367/S1317 selector phase-origin admissibility no-go boundary

`P2367/S1317` continues the selector audit after `P2366/S1316` by separating two facts that must not be conflated:

```text
chiral bispectrum orientation marker  !=  strict-core source origin
```

The finite computation verifies that, for each fixed chiral orientation, the twelve source choices form one free transitive `C12` translation orbit.  Hence any translation-invariant strict-core score is constant on the source orbit and cannot localize `chi11_selector_source`.

The audited candidate classes are:

```text
translation-invariant powers / histograms        -> rejected as source localizers
complete translation-invariant bispectrum        -> source-blind, but useful as chiral orientation marker
non-coprime raw phase modes                      -> rejected as full source localizers because they alias
coprime phase + calibrated phase-origin reference -> best operational candidate, but premise-based
```

This preserves the P2366 positive boundary: coprime phases can recover source only after a phase-origin reference is supplied.  It also preserves the hard limit: no strict-core phase-origin theorem, no `beta_tors -> chi11` theorem, no legacy role transfer, no QW-2191 discharge, and no ToE closure are claimed.

## P2368/S1318 self-recorded endpoint anchor selector candidate

`P2368/S1318` takes the next selector step after the P2367 no-go boundary.  Instead of asking a translation-invariant scalar to choose an absolute source, it audits a conditional **self-recorded ordered d5 endpoint anchor**:

```text
ordered d5 support + self-recorded arrow action -> ledger (2,2,2,1,1)
ledger endpoint values -> source/orientation extractor
```

The finite certificate enumerates all 35 positive five-part ledgers of total 8 and finds a unique lexicographic `(ripple, arrow)` minimizer `(2,2,2,1,1)`.  Given that ordered ledger on the d5 path, the endpoint extractor recovers all 24 source/orientation rows and is D12-equivariant over 576 shift/reflection cases.

This is stronger than a bare Fourier phase-origin premise because the source/orientation are read from structured finite self-record data.  It still does not close QW-2191: strict nadsoliton dynamics has not yet derived the d5 support, the self-recorded arrow action, the endpoint-valued ledger, or the endpoint-source convention.  No `beta_tors -> chi11` theorem, no legacy role transfer, and no ToE closure are claimed.

## P2369/S1319 closed-form ledger uniqueness theorem

`P2369/S1319` strengthens the P2368 ledger step from finite enumeration to a closed-form integer proof.  For five positive integer slots with total 8, convex smoothing gives the ripple lower bound:

```text
sum e_i^2 >= 2*1^2 + 3*2^2 = 14
```

and equality holds exactly on the ten permutations of the multiset `{2,2,2,1,1}`.  The self-recorded arrow penalty

```text
A(e)=sum_i max(0, e_{i+1}-e_i)
```

is zero exactly for the nonincreasing permutation, so the lexicographic `(ripple, arrow)` action uniquely selects `(2,2,2,1,1)`.  The same packet also checks the exact identity `(12^5)*(2^8/3^5)=4^9`, i.e. the ledger product is consistent with `eta=9/5`.

This closes only the finite ledger-uniqueness subclaim.  It still does not derive the ordered d5 support, the self-recorded arrow action, or the endpoint-source convention from strict nadsoliton dynamics; it does not discharge QW-2191 and does not transfer legacy physical roles.

## P2370/S1320 distance-5 band-pass support closed-form theorem

`P2370/S1320` addresses the support-side premise left open by P2369.  It does not derive the band-pass action, but it proves the finite consequence if such an action is admitted.  Because `gcd(5,12)=1`, the distance-5 graph is a single 12-cycle.  Therefore any five-node support has at most four internal distance-5 edges, and equality holds exactly for a connected five-node path in that cycle:

```text
max h5 = 4
argmax h5 = 12 translated distance-5 paths
```

The proof is closed-form and is cross-checked against all `binom(12,5)=792` supports.  Combined with P2369, the finite support and finite ledger uniqueness subclaims are isolated.  The remaining missing premise is now sharper: derive the distance-5 band-pass action and path-order/arrow convention from bridge-completed nadsoliton dynamics, or mark them as non-strict selector/source premises.  This does not choose a unique translate/source and does not discharge QW-2191.

## P2371/S1321 full-Aut unit-band obstruction for band-pass derivation

`P2371/S1321` proves a negative selector-side statement about the remaining P2370 premise.  If a support action is full-`Aut(Z12)` invariant on unit bands, it must weight the distance-1 and distance-5 unit bands equally.  On five-node supports the invariant score

```text
h_unit = h1 + h5
```

is maximized by 24 mixed supports with `(h1,h5)=(3,3)`, not by the 12 distance-5 paths with `(h1,h5)=(0,4)`.  For nonnegative linear scores `a*h1+b*h5`, the distance-5 path orbit is uniquely selected only in the symmetry-breaking chamber:

```text
a/b < 1/3.
```

Thus the distance-5 band-pass polarity cannot be silently obtained from full-`Aut(Z12)` unit-band symmetry.  The next real target is a bridge-completed dynamical export of this polarity inequality, or an explicit non-strict selector/source premise.  QW-2191 remains open.

## P2372/S1322 bridge-kernel direct band-polarity audit

`P2372/S1322` tests the most direct possible bridge-completed source for the P2371 polarity inequality: use the completed strict kernel values themselves as pair weights for `a*h1+b*h5`:

```text
a = K_strict(1),  b = K_strict(5).
```

The result is negative.  With the P2363/P2371 parameters,

```text
K_strict(1)/K_strict(5) ~= 19.476247 >> 1/3,
```

so direct kernel pair weights are far outside the d5 chamber and select distance-1 path supports, not the distance-5 path supports.  The amplitude-normalized legacy kernel is also not a licensed rescue: its `K_legacy(5)` is negative and legacy physical-role transfer remains forbidden.

Thus the band-polarity source is still not exported by direct `K(d)` pair weights.  The honest next target is a separate bridge-completed dynamical term with effective `h1/h5` coefficient ratio `<1/3`, or else an explicit non-strict selector/source premise.

## P2373/S1323 bridge-kernel polarity correction cone

`P2373/S1323` converts the negative P2372 direct-kernel audit into quantitative necessary correction bounds.  Starting from

```text
a0 = K_strict(1),  b0 = K_strict(5),  a0/b0 ~= 19.476247,
```

the P2371 d5 chamber `a/b<1/3` can be reached only by a large extra polarity term.  Closed-form thresholds are:

```text
pure h5 boost:        lambda > 3*a0 - b0
pure h1 suppression:  mu     > a0 - b0/3
antisymmetric +h5-h1: gamma  > (3*a0 - b0)/4
```

The probe verifies that just-above-threshold corrections select the 12 d5 path supports on all 792 supports.  This is not a derivation of such a term; it only quantifies what a bridge-completed dynamical source would have to export.  Without such a source, the band-polarity premise remains non-strict and QW-2191 remains open.

## P2374/S1324 damping-compression polarity candidate

`P2374/S1324` audits the first bridge-side feature found after the P2373 correction cone: the strict nonlinear damping/compression surplus

```text
C(d) = log((1 + beta*d^eta)/(1 + beta_tors*d)).
```

On the audited `d=1,5` pair this gives `C1/C5 ~= 0.235429`, which lies inside the P2371 d5 chamber `a/b<1/3`.  The same probe also derives the exact blend condition against the failed direct strict-kernel weights:

```text
(K1 + tau*C1)/(K5 + tau*C5) < 1/3
iff tau > (3*K1 - K5)/(C5 - 3*C1).
```

For the current canonical constants the threshold is `tau > 1.6259309081595656`, and the just-above-threshold blend selects the 12 d5 path supports on all 792 supports.  This is a candidate source direction, not a strict dynamical source theorem: without a variational/transport theorem promoting `C(d)` to the selector action, QW-2191 remains open and no legacy role transfer is licensed.

## P2375/S1325 damping-compression polarity interval robustness

`P2375/S1325` strengthens P2374 from a single canonical `beta_tors=0.01` check to an interval theorem.  For

```text
C(d) = log((1 + d^(9/5))/(1 + beta_tors*d)),   beta_tors in [0, 0.1],
```

the d5 chamber inequality `C(1)/C(5)<1/3` is equivalent to the positive-margin condition

```text
F(x) = (1 + 5^(9/5))*(1+x)^3 - 8*(1+5*x) > 0.
```

Since `F(0)>0` and `F'(x)=3*(1+5^(9/5))*(1+x)^2-40` is already positive at `x=0` and increasing on the interval, the chamber inequality holds throughout `[0,0.1]`.  Endpoint/canonical/midpoint scans again select the 12 d5 supports on all 792 supports.  This removes a fine-tuning worry for the compression-polarity candidate, but it still does not derive the variational/transport source required to promote `C(d)` into a strict selector action.

## P2376/S1326 damping-compression eta/beta rectangle robustness

`P2376/S1326` extends P2375 from a one-parameter `beta_tors` interval to the rectangle

```text
eta in [9/5, 2],   beta_tors in [0, 0.1].
```

For `C(d)=log((1+d^eta)/(1+beta_tors*d))`, the d5 chamber condition `C(1)/C(5)<1/3` is again equivalent to

```text
F(eta,x) = (1 + 5^eta)*(1+x)^3 - 8*(1+5*x) > 0.
```

The proof checks that `F` is increasing in both `eta` and `x` on the audited rectangle, so its minimum is at `(eta,x)=(9/5,0)`, where the margin is already positive.  A 3x3 grid scan over all 792 supports confirms the d5 path selection throughout the sampled rectangle.  This closes an eta/beta fine-tuning worry for the candidate, but it still does not supply the missing variational/transport source theorem.

## P2377/S1327 damping-compression transport primitive and uniform coupling

`P2377/S1327` stops widening finite robustness and tests the first requested source-level interface for the robust compression candidate.  It proves that

```text
C(d)=log((1+d^eta)/(1+beta_tors*d))
```

is the exact endpoint primitive of the damping-completion log-transport one-form along

```text
u_s(d)=(1-s)*(1+beta_tors*d)+s*(1+d^eta),
A_s(d)=partial_s log(u_s(d)).
```

Thus `integral_0^1 A_s(d) ds = C(d)`.  This gives transport provenance for the P2374 feature rather than a bare pair-feature insertion.

The same theorem computes a rectangle-uniform acceptance threshold for blending this transport primitive with the failed direct strict pair weights:

```text
(K1+tau*C1)/(K5+tau*C5)<1/3
```

through the denominator `D(eta,x)=C5-3*C1`.  On the P2376 rectangle, `D` is increasing in `eta` and decreasing in `x in [0,0.1]`, so its minimum is at `(eta,x)=(9/5,0.1)`.  Therefore a single scalar threshold `tau>(3*K1-K5)/D(9/5,0.1)` is sufficient on the audited rectangle, and grid support scans again select the 12 d5 supports.

This is progress, but it is not the missing strict variational source theorem: the scalar coupling `tau` is not derived from nadsoliton dynamics, `C(d)` is not promoted into `L_total`, and QW-2191 remains open.

## P2378/S1328 unit-normalized transport coupling insufficiency

`P2378/S1328` audits the strongest false-closure risk left by P2377.  P2377 gives exact transport provenance for `C(d)`, but if that endpoint primitive is inserted with only unit-normalized total transport mass `M<=1`, it still cannot overcome the failed direct strict-kernel pair polarity.

For a scalar transport budget

```text
a = K1 + M*C1,   b = K5 + M*C5,
```

the d5 chamber requires

```text
M > (3*K1-K5)/(C5-3*C1).
```

On the P2376 rectangle, `D(eta,x)=C5-3*C1` is positive, increasing in `eta`, and decreasing in `x`; its maximum occurs at `(eta,x)=(2,0)`.  Since even this maximum is smaller than `3*K1-K5`, every admissible point has threshold greater than one.  The certified threshold range is approximately

```text
1.1757688203 < M_threshold < 1.8435099396.
```

Thus P2377 transport provenance is not enough by itself: a future theorem must derive a super-unit source normalization above the relevant threshold, or else the extra coupling strength remains an explicit non-strict selector premise.  Grid scans confirm that `M=1` selects the mixed `(h1,h5)=(3,3)` orbit, while just-above-threshold mass selects the 12 d5 supports.

## P2379/S1329 front-loaded normalized transport profile candidate

`P2379/S1329` narrows the P2378 obstruction.  The unit-uniform primitive `C(d)` is insufficient, but a normalized nonuniform profile on the same P2377 homotopy can be strong enough because it samples the early high-contrast part of the damping-completion transport.

For

```text
rho_lambda(s)=1+lambda*(1-2s),  0<=lambda<=1,
```

one has

```text
int_0^1 rho_lambda(s) ds = 1,
int_0^1 rho_lambda(s) A_s(d) ds = C(d)+lambda*B(d),
B(d)=C(d)*(1+2*(1+beta_tors*d)/(d^eta-beta_tors*d))-2.
```

The d5 threshold in this affine profile family is

```text
lambda > [3*(K1+C1)-(K5+C5)]/(B5-3*B1).
```

The P2379 audit tests `lambda=0.8`: uniform normalized transport still fails, but the front-loaded normalized density selects the 12 d5 supports on the audited grid and stays above the computed 81x81 lattice threshold.  This does not promote the profile into `L_total`; it only records a sharper future source obligation: derive either super-unit mass or a sufficiently front-loaded normalized density from strict dynamics.

## P2380/S1330 rectangle certificate for the affine front-loaded threshold

`P2380/S1330` converts the P2379 affine-profile lattice audit into a rectangle monotonicity certificate.  For

```text
T(eta,x)=[3*(K1+C1)-(K5+C5)]/(B5-3*B1),
```

with the same `C(d)` and `B(d)` as P2379, interval arithmetic over the P2376 rectangle certifies

```text
T_eta<0,  T_x>0,  B5-3*B1>0.
```

Therefore the maximum threshold occurs at `(eta,x)=(9/5,0.1)` and equals approximately `0.7916644842269442`; `lambda=0.8` is rectangle-uniformly sufficient for the audited affine density family.  The equation-sheet status remains source-open: this is not an `L_total` promotion, only a sharper acceptance condition for a future theorem deriving the normalized front-loaded density from strict dynamics.

## P2381/S1331 source-obligation form of the affine frontload result

`P2381/S1331` records the source burden implied by P2380.  Rectangle-uniform sufficiency of the affine density family requires more than normalization:

```text
lambda > 0.7916644842269429.
```

For `rho_lambda(s)=1+lambda*(1-2s)`, this forces

```text
early-half mass > 0.6979161,
transport barycenter < 0.3680559,
rho(0)/rho(1) > 8.5995.
```

The equation-sheet consequence is a sharpened non-closure statement: a source theorem must derive this front-loading burden, or the affine profile remains a quantified non-strict premise.  P2381 does not promote the profile into `L_total` and does not close QW-2191.

## P2382/S1332 bounded-density bathtub frontload certificate

`P2382/S1332` generalizes the affine-only profile question to all normalized transport densities with a pointwise cap

```text
0 <= rho(s) <= M,   int_0^1 rho(s) ds = 1.
```

For the P2377 log-transport one-form, the relevant d5-vs-h1 contrast is

```text
q(s)=A_s(5)-3*A_s(1).
```

Interval arithmetic certifies `q'(s)<0` on the P2376 rectangle.  Therefore the bathtub rearrangement principle reduces the optimal `M`-capped normalized density to the early bang-bang profile

```text
rho_M(s)=M for 0<=s<=1/M,  rho_M(s)=0 otherwise,
W_M(d)=M*log((1+beta_tors*d+(d^eta-beta_tors*d)/M)/(1+beta_tors*d)).
```

The worst audited cap threshold is at `(eta,beta_tors)=(9/5,0.1)`:

```text
M > 1.5748213574353633.
```

Thus `M=1.6` is rectangle-uniformly sufficient and selects the 12 d5 supports, while a just-below-threshold cap at the worst corner fails the d5 chamber.  This is still not a strict source theorem: it only says that any future bounded-density source may be checked against the explicit cap/frontload burden, and that the optimal candidate under such a cap is known.  No `L_total` promotion, legacy role transfer, `beta_tors -> chi11` theorem, QW-2191 discharge, selector closure, or ToE closure follows.
