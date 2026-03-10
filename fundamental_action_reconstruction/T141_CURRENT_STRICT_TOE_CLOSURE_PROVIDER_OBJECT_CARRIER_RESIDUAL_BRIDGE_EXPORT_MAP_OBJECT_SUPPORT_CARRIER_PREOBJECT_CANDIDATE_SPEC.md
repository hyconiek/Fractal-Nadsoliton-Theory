# T141 Current Strict ToE Closure Provider-Object Carrier Residual Bridge/Export-Map Object-Support Carrier Preobject Candidate Spec

Status: `T141_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_PREOBJECT_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

After `N405` and `N407`, the strict closure lane is in a stable situation:

1. provider-object residual route reaches the bridge/export-map object-support
   **witness** frontier (`N405`),
2. the missing post-witness **object-support carrier** layer is now sharply
   named as a future-only target (`T140/N407`),
3. `N302` still blocks any honest promotion into **actual** object support
   using only the already exported material.

So the next honest move is not “closure”.
It is weaker and audit-safe:

```text
export one explicit candidate *preobject* construction for an object-support
carrier that is (a) noncyclic, (b) observer-free, (c) pair-indexed, and
(d) bridge-facing,
without claiming that it discharges N302 or exports an actual carrier.
```

`T141` specifies such a preobject candidate form on the already exported
provider-object carrier orbit-quotient lane (`T126/F279/N391`), using only
finite nad12-depth orbit summaries and without theta inputs.

## Candidate preobject object

Export one explicit preobject candidate:

```text
Pi_residual_datum_provider_object_carrier_bridge_export_map_object_support_carrier_preobject_candidate_v1
```

with intended type:

```text
Pi_residual_datum_provider_object_carrier_bridge_export_map_object_support_carrier_preobject_candidate_v1 :
  (Pi_strict_provider_object_carrier_orbit_quotient_candidate_v1, a, b)
    -> ProviderObjectResidualObjectSupportCarrierPreobject_pair^{cand}
```

where `a,b` are the strict contraction parameters already carried by the
provider-object carrier candidate (`F279`), and the codomain is a concrete,
pair-indexed, gauge-aware, finite-depth orbit summary intended to be
bridge-facing under the already exported projection candidate interface
(`T127/F280`).

## Construction (finite-depth orbit summary; noncyclic)

### 1. Carrier space and gauge/symmetry group

Reuse the carrier choice and explicit symmetry from `T126`:

```text
X_carrier := ℓ^2(ℤ)
G_phase := U(1) acting by (g·ψ)_k := g ψ_k
PairIndex_v1 := {pair1, pair2}
```

### 2. Provider operators (pair-indexed contractions)

Reuse the pair-indexed contraction operators from the provider-object carrier
candidate lane (`T126/F279`):

```text
(P_pair1 ψ)_k := a ψ_{k-1},  0<|a|<1
(P_pair2 ψ)_k := b ψ_{k+1},  0<|b|<1
```

and reuse the declared seed state:

```text
ψ_src := δ_0 ∈ ℓ^2(ℤ)
```

### 3. Finite nad12-depth weights (explicit; reuses the same noncyclic carrier)

Define `K := 11` (nad12 depth).

For each pair slot `i ∈ {1,2}` define:

```text
r_1 := |a|
r_2 := |b|
w_{i,k} := r_i^k / (sum_{j=0..K} r_i^j),   k=0..K.
```

This yields a finite noncyclic carrier field `E_pair^{cand,prov}` already used
by the bridge-facing projection candidate (`T127`).

### 4. Normalized orbit representatives (avoid double-damping)

For each pair slot `i` and depth `k`, define a normalized orbit representative:

```text
q_{i,k} := P_pairi^k(ψ_src) / ||P_pairi^k(ψ_src)||.
```

Because `0<|a|,|b|<1` and `ψ_src ≠ 0`, the norms are explicit:

```text
||P_pair1^k(ψ_src)|| = |a|^k,
||P_pair2^k(ψ_src)|| = |b|^k,
```

so `q_{i,k}` is well-defined and has `||q_{i,k}|| = 1`.

### 5. Preobject carrier vectors (finite convex combination)

Define the pair-indexed orbit-summary vectors:

```text
ψ_support_pairi^{cand} := sum_{k=0..K} w_{i,k} q_{i,k}  ∈ X_carrier.
```

This produces two explicit nonzero vectors:

```text
(ψ_support_pair1^{cand}, ψ_support_pair2^{cand})
```

and yields a concrete pair-indexed carrier preobject:

```text
ProviderObjectResidualObjectSupportCarrierPreobject_pair^{cand} :=
(
  ψ_support_pair1^{cand},
  ψ_support_pair2^{cand}
)
```

Optionally record gauge awareness by keeping the global-phase quotient
statement as metadata:

```text
g·(ψ1,ψ2) := (gψ1, gψ2),  g∈U(1).
```

## Bridge-facing intent (still below N302)

This preobject candidate is declared “bridge-facing” only in the weak sense:

1. it is built on the same pair-indexed finite-depth carrier `E_pair` already
   used to feed the residual projection candidate (`T127/F280`),
2. it can be attached to the existing projection candidate as an auxiliary
   carrier field, without claiming any actual object support.

No claim is made that this preobject is an actual object-support carrier,
and no claim is made that it satisfies `T140` acceptance tests.

## Noncyclic / observer-free discipline

`T141` is admissible only if:

1. no `theta_{1,2}` are used as inputs,
2. no populated basis-pair instance is used as input,
3. no `K_obs`-indexed selection is used as a primary source,
4. `U(1)` gauge/symmetry is treated as an internal symmetry action and not as
   an observer source,
5. no discharge of `N302` is implied.

## Hard limits

`T141` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. any export-map object export (`N300`),
4. actual theta export / pair population (`N18`),
5. actual `E_orient`,
6. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
7. ToE closure.

