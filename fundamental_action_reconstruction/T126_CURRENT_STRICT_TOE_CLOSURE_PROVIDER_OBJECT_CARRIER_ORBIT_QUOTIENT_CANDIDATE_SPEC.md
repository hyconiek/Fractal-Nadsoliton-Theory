# T126 Current Strict ToE Closure Provider-Object Carrier Orbit-Quotient Candidate Spec

Status: `T126_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_ORBIT_QUOTIENT_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`N327` diagnoses the dominant missing ingredient class for strict ToE closure as
one:

```text
source-side, observer-free, pair-indexed, noncyclic
strict selector/provider object-carrier layer.
```

`T125` names that missing layer as a **future-only target** object:

```text
Epsilon_strict_provider_object_carrier_layer_target_v1
```

`T126` does **not** claim discharge of that target.

`T126` does something narrower and audit-safe:

```text
export one explicit orbit-quotient carrier construction as a candidate form
with an explicit gauge/symmetry group (no "observer" ambiguity),
explicit pair indexing,
and an explicit noncyclicity mechanism (strict contraction).
```

This candidate is intentionally **not** claimed to be welded to the full
nad12-sigma residual semantics (`N328-N344`) or to discharge the `N302` bridge
object-support frontier. It is a concrete mathematical form that can later be
either:

1. welded to the existing strict lane by additional theorems, or
2. rejected as non-matching to the strict carrier requirements.

## Candidate construction (operator orbit-quotient carrier)

### 1. State space (carrier space)

Fix a candidate carrier Hilbert space:

```text
X_carrier := ℓ^2(ℤ).
```

This is a *carrier choice*, not a claim that the nadsoliton state space has
already been identified with `ℓ^2(ℤ)`.

### 2. Gauge/symmetry group (explicit; replaces "observer" wording)

Fix a minimal gauge/symmetry group:

```text
G_phase := U(1)
```

acting on `X_carrier` by global phase multiplication:

```text
(g · ψ)_k := g ψ_k.
```

This implements the “observer-free” discipline as an explicit invariance under
an internal symmetry action (no `K_obs`-indexed selection).

### 3. Pair index set

Fix a minimal designated strict pair-index set:

```text
PairIndex_v1 := {pair1, pair2}.
```

### 4. Provider operators (pair-indexed, strict contraction)

Choose pair-indexed bounded linear operators on `X_carrier`:

```text
P_pair1, P_pair2 : X_carrier -> X_carrier
```

of the explicit shift+contraction form:

```text
(P_pair1 ψ)_k := a ψ_{k-1},
(P_pair2 ψ)_k := b ψ_{k+1},
```

with strict contraction parameters:

```text
0 < |a| < 1,
0 < |b| < 1.
```

### 5. Noncyclicity witness (mechanism)

Because `∥P_pairi ψ∥ = |a_i| ∥ψ∥` for `ψ ≠ 0` (with `a_1=a`, `a_2=b`), the
construction admits one explicit noncyclicity mechanism:

```text
P_pairi^n(ψ) = ψ for some n>=1 implies ψ = 0.
```

This is a mathematical noncyclicity witness (no theta inputs, no populated
instance inputs).

### 6. Orbit-quotient relation (gauge-aware)

For each `pairi ∈ PairIndex_v1` define an equivalence relation on
`X_carrier \ {0}`:

```text
ψ ~_pairi φ
  :⇔  ∃ m,n >= 0, ∃ g ∈ G_phase such that
       P_pairi^m(ψ) = g · P_pairi^n(φ).
```

Denote the orbit-quotient class by:

```text
[ψ]_pairi ∈ (X_carrier \ {0}) / ~_pairi.
```

### 7. Source-side seed (minimal; reuses tau_src_candidate_v1 only as a tag)

Reuse only the already exported strict source packet:

```text
tau_src_candidate_v1 (F127/N235)
```

and define one minimal seed state representative:

```text
ψ_src := δ_0 ∈ ℓ^2(ℤ)
```

with the **declared** intended reading:

```text
ψ_src is a carrier representative for the source-side tag tau_src_candidate_v1,
not a discharged identification of tau_src_candidate_v1 with a physical state.
```

### 8. Provider-object carrier type (concrete)

Define the concrete pair-indexed orbit-quotient carrier type:

```text
ProviderObjectCarrier_pair^{cand,orbit} :=
(
  [ψ_src]_pair1,
  [ψ_src]_pair2
)
```

with attached metadata:

```text
(
  X_carrier,
  G_phase,
  PairIndex_v1,
  P_pair1,P_pair2,
  tau_src_candidate_v1_tag
).
```

### 9. Candidate carrier-layer map (typed)

Define the candidate carrier-layer map:

```text
C_strict_provider_object_carrier_orbit_quotient_candidate_v1 :
  tau_src_candidate_v1 -> ProviderObjectCarrier_pair^{cand,orbit}
```

by:

```text
C_strict_provider_object_carrier_orbit_quotient_candidate_v1(tau_src_candidate_v1)
  := ProviderObjectCarrier_pair^{cand,orbit}  (with tau_src tag attached).
```

## Scope / non-claims

`T126` is **not**:

1. an actual provider-object realization theorem,
2. a discharge of `Epsilon_strict_provider_object_carrier_layer_target_v1`,
3. a discharge of `N302`,
4. an export-map object,
5. an admissible `S_sel_int`,
6. a `QW-2191` discharge,
7. a ToE closure move.

`T126` is only an explicit candidate *form* for a future provider-object carrier
layer, expressed with:

1. an explicit gauge/symmetry group,
2. an explicit orbit-quotient carrier type,
3. an explicit noncyclicity mechanism.

