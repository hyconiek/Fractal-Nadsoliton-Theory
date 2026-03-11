# T146 Current Strict Residual Datum Bridge/Export-Map Object-Support Carrier Actual Inhabitant Spec

Status: `T146_CURRENT_STRICT_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_ACTUAL_INHABITANT_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

The strict residual-datum bridge/export-map lane currently exports:

1. two witness-level arrivals at the object-support frontier (`N387`, `N405`),
2. one explicit post-witness object-support carrier target (`T140/N407`),
3. one provider-object carrier lane with source-derived contraction parameters
   (`N410`) and an actual carrier-layer inhabitant (`N412`),

while still remaining below:

- actual bridge/export-map object support above the map object (`N395`),
- any claim of selector closure or `QW-2191` discharge.

On the updated repo state, the strict-core bridge/export-map object is already
exported (`F311/N422`), so `N300/N301` are historical.

After `N412`, the next honest move on the residual branch is no longer to name
targets. It is to attempt a **post-witness carrier discharge**:

```text
discharge Omicron_residual_datum_bridge_export_map_object_support_carrier_target_v1
on a declared strict domain, at least on the provider-object witness track,
without implying actual object support above the map object or closure.
```

`T146` specifies one minimal carrier inhabitant construction that is:

1. noncyclic (no theta inputs; no populated-instance inputs),
2. observer-free (no `K_obs` as a primary selector source),
3. pair-indexed,
4. explicitly above the provider-object witness layer,
5. explicitly below actual bridge/export-map object support above the map object
   (`N395` remains future-only).

## Output object

Export one actual post-witness object-support carrier inhabitant:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_v1
```

with a declared upgrade scope:

```text
upgrades: Kappa_residual_datum_provider_object_carrier_bridge_export_map_object_support_witness_v1  (N405)
does_not_claim_upgrade_of: Kappa_residual_datum_sigma_int_bridge_export_map_object_support_witness_v1  (N387)
```

This avoids any silent use of `sigma_int_candidate` prerequisites (`T123/T124`).

## Construction (finite nad12-depth orbit-summary carrier; source-derived)

### 1. Source-derived contraction parameters

Reuse:

```text
A_strict_provider_object_contraction_parameter_source_map_v1 :
  tau_src_candidate_v1 -> (a,b)
```

from `N410`, so:

```text
(a,b) := (cos(phi), cos(phi)),
0 < |a| < 1, 0 < |b| < 1  (on the declared strict core branch).
```

### 2. Carrier space, gauge group, pair index

Reuse the already exported carrier choices from the provider-object lane:

```text
X_carrier := ℓ^2(ℤ)
G_phase   := U(1)
PairIndex_v1 := {pair1, pair2}.
```

### 3. Provider operators and seed

Define (as in `T126`):

```text
(P_pair1 ψ)_k := a ψ_{k-1},
(P_pair2 ψ)_k := b ψ_{k+1},
ψ_src := δ_0.
```

### 4. Finite nad12-depth weights

Let `K := 11` (nad12 depth) and define:

```text
r_1 := |a|
r_2 := |b|
w_{i,k} := r_i^k / (sum_{j=0..K} r_i^j),  k=0..K.
```

### 5. Normalized orbit representatives

Define:

```text
q_{i,k} := P_pairi^k(ψ_src) / ||P_pairi^k(ψ_src)||.
```

Because `0<|a|,|b|<1` and `ψ_src ≠ 0`, `q_{i,k}` are well-defined and satisfy
`||q_{i,k}||=1`.

### 6. Post-witness carrier vectors (finite convex combinations)

Define the pair-indexed carrier vectors:

```text
ψ_support_pairi := sum_{k=0..K} w_{i,k} q_{i,k}  ∈ X_carrier.
```

Package:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_v1 :=
(
  ψ_support_pair1,
  ψ_support_pair2,
  metadata: (X_carrier, G_phase, PairIndex_v1, K, (a,b)),
  upgrades_witness: Kappa_residual_datum_provider_object_carrier_bridge_export_map_object_support_witness_v1
).
```

### 7. Bridge-facing projection interface (still below export)

Attach the already exported projection layer (`N398`):

```text
Xi_residual_datum_provider_object_carrier_bridge_export_map_object_support_projection_v1
```

as the declared bridge-facing interface. This does not export any bridge/export
map object (exported by `F311/N422`) and does not export any actual object
support above the map object (`N395` remains future-only).

## Acceptance alignment (T140 / T139)

`T146` is intended to discharge:

1. `Omicron_residual_datum_bridge_export_map_object_support_carrier_target_v1`
   (`T140/N407`), and
2. `Lambda_residual_datum_provider_object_carrier_bridge_export_map_object_support_target_v1`
   (`T139/N406`),

by exporting one typed post-witness carrier object above `N405` and below
actual bridge/export-map object support above the map object (`N395`), with
explicit noncyclic + observer-free contracts.

## Hard limits

`T146` must not claim:

1. any actual bridge/export-map object support above the map object (`N395`),
2. discharge of sigma-int strict prerequisites (`T123/T124`),
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.
