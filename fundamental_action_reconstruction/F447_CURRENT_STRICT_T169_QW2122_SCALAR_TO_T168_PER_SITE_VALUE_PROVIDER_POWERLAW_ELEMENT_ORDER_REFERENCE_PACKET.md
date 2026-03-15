# F447 Current Strict `T169` QW‑2122 Scalar → `T168` Per‑Site Value‑Provider (Power‑Law Element‑Order Reference) Packet (No False‑PASS)

Status: `F447_EXECUTED_CURRENT_STRICT_T169_QW2122_SCALAR_TO_T168_PER_SITE_VALUE_PROVIDER_POWERLAW_ELEMENT_ORDER_REFERENCE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`T166` (diagonal/local `pair1` mode‑2 defect decision) is computable only after the repo can produce strict‑derived numeric
values for:

```text
vpsi[0..11], g4[0..11], g6[0..11]          (T168)
→ Sigma_psi{k}_psi{k+6} for k=0..5         (T167 via N477/P437)
→ F2(d) and theta_*                        (T166 via P434/N468)
```

`N478/P443` state that the strict scalar closure exported by `QW‑2122/2123/2124` does **not** canonically lift to the
per‑site families required by `T168`, so a dedicated lift/value‑provider object is required (`T169`).

This packet exports one explicit deterministic constrained‑lift provider which:

1. reuses only strict‑admissible inputs already exported in the repo (`QW‑2122`, `alpha_geo_strict_derived_v1`, `R14`, `R15`,
   typed `Z_12` language from `F329`), and
2. produces a `T168`‑consumable value object (JSON) and auto‑populates the designated harness inputs (`P437` and `P434`)
   without any hidden generator/direction slot.

It does **not** claim global `QW‑2191` discharge or ToE closure.

## Strict-admissible inputs reused

1. `QW‑2122` (`report_qw2122_psi_potential_diagonal_floor_gate.json`)
   - `rho_star_sq` (broken branch),
   - `lambda_psi_strict`,
   - diagonal floor `m0^2` (embedded via `R15`).
2. `F309/N420`
   - `alpha_geo_strict_derived_v1 := 4 ln 2` (strict source‑upgrade datum).
3. `F329`
   - typed `Z_12` carrier language.
4. `N479`
   - `ord_Z12(x)` is `Aut(Z_12)`‑invariant ⇒ no marked generator/direction slot for references of the form `f(ord_Z12(x))`.
5. `N483`
   - theorem-level existence + uniqueness support for this exact constrained lift (no hidden slots; explicit residual fixing).
6. `R14`
   - strict specialization of the canonical symmetrized kernel channel to frozen numeric `K_total` on the 12‑slot carrier.
7. `R15`
   - strict host scalar floor embedding `m0^2 I` in the canonical diagonal decomposition.

## Exported strict reference datum (no marked direction)

Define the element order on `Z_12`:

$$
\operatorname{ord}_{Z_{12}}(x)=\frac{12}{\gcd(x,12)}\quad (x\neq 0),\qquad \operatorname{ord}_{Z_{12}}(0)=1.
$$

Define power‑law weights:

$$
w_x := \operatorname{ord}_{Z_{12}}(x)^{-\alpha_{\mathrm{geo}}},
\qquad
\alpha_{\mathrm{geo}} := 4\ln 2,
$$

and normalize:

$$
r_{\mathrm{ordpow}}(x) := \frac{w_x}{\sum_{y\in Z_{12}}w_y}.
$$

Exported artifact:

- `fundamental_action_reconstruction/generated/r_ordpow_z12_v1_reference_distribution.json`

Notes:

1. `r_ordpow` is `Aut(Z_12)`‑invariant (direction‑free) by `N479`.
2. `r_ordpow` is **not** translation‑invariant on the regular action (it distinguishes the identity orbit `{0}`).
   This is an explicit symmetry‑breaking datum class; it is not smuggled.

## The constrained lift (T169 → T168)

### 1) Magnitude lift (scalar norm)

Use the strict scalar vacuum scale from `QW‑2122` (broken branch):

$$
\rho_*^2 := \texttt{rho\_star\_sq}.
$$

Lift to per‑site magnitudes:

$$
|v_{\psi i}| := \sqrt{\rho_*^2\,r_{\mathrm{ordpow}}(i)}.
$$

### 2) Deterministic sign selection (no hidden slot)

Given `K_total` from `R14` and magnitudes `|v_{\psi i}|`, select a sign vector $s\in\{\pm 1\}^{12}$ by minimizing the
kernel mixing energy:

$$
E_{\mathrm{mix}}(s) := \sum_{i\lt j} K_{\mathrm{total}}[i,j]\,(s_i|v_{\psi i}|)\,(s_j|v_{\psi j}|).
$$

Since $2^{12}=4096$, this is solved by a deterministic exhaustive search (no solver tolerance/seed slot).

The unavoidable global `Z2` ambiguity is fixed deterministically by requiring $s_0=+1$ (the identity orbit is already
distinguished by the reference datum), and then selecting the lexicographically smallest minimizer.

Finally:

$$
v_{\psi i} := s_i\,|v_{\psi i}|.
$$

### 3) Quartic/sextic lift (scalar→per‑site matching along the lifted vacuum direction)

Set:

$$
g6_{\psi i} := 0,
$$

and define a uniform quartic:

$$
g4_{\psi i} := g4_{\mathrm{uniform}}
            := \frac{\lambda_{\psi}}{\sum_{i=0}^{11} r_{\mathrm{ordpow}}(i)^2},
$$

so that the quartic coefficient along the lifted vacuum direction matches the scalar `QW‑2122` quartic coefficient
`lambda_psi_strict` exactly.

## Exported provider + harness inputs

Provider (T168‑consumable arrays + computed sigma sums):

- `fundamental_action_reconstruction/generated/f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json`

Designated inputs auto‑populated:

- `fundamental_action_reconstruction/generated/p437_input_vpsi_g4_g6_candidate.json`
- `fundamental_action_reconstruction/generated/p434_input_sigma_opposite_pair_sum_values_candidate.json`

## Current computed outputs (repo-state, reproducible; not a global closure claim)

Running:

- `P437` produces a numeric six‑sum instantiation via the `N477` form:
  - `F2_abs ≈ 12.88048321986276`.
- `P434` evaluates:
  - `cuts_O2_on_pair1_by_N466 = true`,
  - `theta_star_by_N468 ≈ 0 (mod π)` on the current exported inputs.

Artifacts:

- `fundamental_action_reconstruction/generated/p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json`
- `fundamental_action_reconstruction/generated/p434_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_evaluation_probe_summary.json`
- `fundamental_action_reconstruction/generated/p448_current_strict_f447_t169_provider_provenance_audit_probe_summary.json`

## Hard limits (no false‑PASS)

This packet does **not** claim:

1. global `QW‑2191` discharge (beyond the scoped diagonal/local `pair1` accelerator lane),
2. strict-core selector closure,
3. ToE closure.

It only exports one explicit strict value-provider/lift object sufficient to make the diagonal accelerator lane
computable and auditable on current repo exports.
