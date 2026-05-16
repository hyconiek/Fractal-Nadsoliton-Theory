# THEORY STATE SUMMARY (STRICT-ONLY CHAIN, P1864)

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Scope: strict-only lane (`no legacy bridge transfer claims`)  
Date anchor: 2026-05-16

## 1) Executive state

Current repository state exports a strict-only forward/reverse chain:

```text
K_strict -> coefficients -> full non-skeleton L_total -> covariant EOM blocks
-> renormalization / BRST / Cutkosky / background-independence witness obligations
```

The chain is explicit and auditable, but theorem-grade closure remains open.

---

## 2) Kernel and coefficient layer

### 2.1 Strict kernel anchor

```text
K_strict(d) = cos(omega*d + phi) / (1 + beta*d^eta)
```

Working strict tuple used in the exported chain:

- `omega = 0.18575`
- `phi = 0.16250`
- `beta = 1.0`
- `eta = 1.8`
- strict-side source upgrade used in coefficients: `alpha_geo_strict = 4 ln(2)`.

### 2.2 Coefficient map schema (strict)

```text
gravity: c_gr_i = C_gr_i[K_strict, alpha_geo_strict, beta, eta, omega, phi]
gauge:   c_g_i  = C_g_i[K_strict, beta, eta, spectral_inputs]
fermion: y_f,m_f = C_f[K_strict, selector_constraints]
higgs:   lambda_H,mu_H2 = C_H[K_strict, stability_inputs]
```

---

## 3) Full non-skeleton Lagrangian (strict)

Repo exports explicit sector-level non-skeleton structure:

```text
L_total = L_GR + L_gauge + L_fermion + L_Higgs + L_mix + L_int + L_covariant
```

with the explicit gravity density seed:

```text
L_GR_strict_nonproxy_v1
= sqrt(-g) * [ (M_Pl^2/2)R
             + c_gr_1 R^2
             + c_gr_2 R_{mu nu}R^{mu nu}
             + c_gr_3 R_{mu nu rho sigma}R^{mu nu rho sigma}
             + c_gr_4 G_GB ]
```

and standard strict blocks exported for gauge/fermion/higgs/mix/int/covariant composition.

---

## 4) EOM layer (forward direction)

### 4.1 Metric/gravity EOM

```text
E_{mu nu} := delta S_total / delta g^{mu nu} = 0
```

expanded seed form:

```text
(M_Pl^2/2) G_{mu nu}
+ c_gr_1 H^(R2)_{mu nu}
+ c_gr_2 H^(Ric2)_{mu nu}
+ c_gr_3 H^(Riem2)_{mu nu}
+ c_gr_4 H^(GB)_{mu nu}
= T^{(matter+gauge+higgs+mix)}_{mu nu}
```

with Bianchi consistency target:

```text
nabla^mu E_{mu nu} = 0
```

### 4.2 Exported componentwise curvature variations

```text
H^(R2)_{mu nu}
H^(Ric2)_{mu nu}
H^(Riem2)_{mu nu}
H^(GB)_{mu nu}
```

are explicitly present as strict componentwise objects for residual assembly.

---

## 5) Renormalization layer (1-loop strict B1)

### 5.1 Divergence projection basis

```text
Gamma_1loop^div = (1/epsilon) * ∫sqrt(-g)
[ a_R2 R^2 + a_Ric2 R_{mu nu}R^{mu nu} + a_Riem2 R_{mu nu rho sigma}R^{mu nu rho sigma} + a_GB G_GB ]
```

### 5.2 Counterterm contract

```text
delta_c_gr_1 = -(1/epsilon) a_R2
delta_c_gr_2 = -(1/epsilon) a_Ric2
delta_c_gr_3 = -(1/epsilon) a_Riem2
delta_c_gr_4 = -(1/epsilon) a_GB
```

```text
c_gr_i^ren(mu) = c_gr_i^bare + delta_c_gr_i(mu)
```

### 5.3 B1 symbolic coefficient layer

Exported symbolic forms:

```text
a_R2   = (alpha_geo_strict/(16*pi^2))*(beta + eta/2)*(omega^2 + phi^2)
a_Ric2 = (alpha_geo_strict/(16*pi^2))*(1+beta)*(eta*omega + phi)
a_Riem2= (alpha_geo_strict/(16*pi^2))*(beta*eta)*(omega^2 + omega*phi + phi^2)
a_GB   = (alpha_geo_strict/(16*pi^2))*(eta-1)*(omega-phi)^2
```

and a sympy-evaluated numeric seed set is present in `P1853`.

---

## 6) BRST / anomaly seed lane (reverse obligations)

### 6.1 BRST anomaly master seed

```text
S_BRST(Gamma_eff_B1) = A_B1
Target: A_B1 = 0
```

### 6.2 Cohain/anomaly polynomial seed

```text
A_B1 = k1*c*R^2 + k2*c*Ricci^2 + k3*c*Riemann^2 + k4*c*G_GB + k5*c*Tr(F^2)
```

where `k1..k4` are linked to gravity-side evaluated coefficients and `k5` is on the explicit gauge-fermion triangle computation lane.

### 6.3 k5 triangle seed

```text
k5 = sum_f Tr[T_a {T_b,T_c}] * Y_f * I_triangle_f(B1)
```

with strict representation/hypercharge seed table and regularization contract exported, but full family/chiral completion remains open.

---

## 7) Cutkosky/unitarity seed lane (reverse obligations)

### 7.1 Optical theorem target

```text
2 Im M_{i->i} = sum_f ∫ dPi_f |M_{i->f}|^2
```

first strict channel focus:

```text
graviton -> gauge_gauge
```

### 7.2 Explicit kernel contract

```text
Disc M(s) = ∫ dPi_gg K_cut(s,theta,phi; c_gr_i^ren, y_f, g_J)
```

plus sample seed kernel table and projected discontinuity rows are exported.

### 7.3 Pole/residue physical projection

Exported rule:

```text
retain poles with s_pole >= 0 and residue >= 0
exclude ghost/negative-residue poles
```

This produces a projected discontinuity seed certificate with uncertainty-banded positivity report (seed-level only).

---

## 8) TG gating state

- TG1 remains preflight-gated by verdict/state lock + S1 evidence completeness.
- TG2/TG3 remain gated by missing theorem-grade BRST/Cutkosky/background-independence witnesses.
- Current strict lane prevents silent pass promotion and keeps all unresolved closures as `OPEN_OBSTRUCTION_WITH_TRACE`.

---

## 9) Exact current blockers to strict-core closure

1. **Renormalization theorem-grade closure**: exact divergence cancellation witness beyond symbolic/seed contracts.  
2. **BRST nilpotency/cohomology theorem**: full anomaly cancellation (`A_B1=0`) with completed `k5` family/chiral computation.  
3. **Cutkosky unitarity theorem**: exact discontinuity integrals and dressed propagator residue certification.  
4. **Background-independence lift**: extension from B1 seed corridor to global theorem-grade background-family closure.

---

## 10) Next honest step (strict)

`P1864` recommendation:

```text
Replace seed uncertainty and seed residues with exact phase-space-integral error propagation + computed dressed pole residues,
then export one unified TG3 candidate witness object:
(physical-pole-projected Disc M >= 0 with uncertainty control + ghost exclusion trace).
```

This is the minimum credible step before any strict-core TG3 closure claim.

---

## 11) Lay explanation

Mamy już kompletną mapę teorii strict: od kernela, przez współczynniki i pełny lagranżian, do równań ruchu i testów kwantowych.  
Największa praca, która została, to zamiana „wersji roboczych” (seed) na pełne, ścisłe obliczenia i dowody (renormalizacja, unitarność, niezależność od tła).  
Dopiero wtedy można uczciwie mówić o domknięciu ToE.
