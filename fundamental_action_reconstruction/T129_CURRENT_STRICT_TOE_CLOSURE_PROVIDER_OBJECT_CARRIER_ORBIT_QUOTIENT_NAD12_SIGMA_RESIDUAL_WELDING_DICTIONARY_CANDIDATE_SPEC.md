# T129 Current Strict ToE Closure Provider-Object Carrier Orbit-Quotient ↔ Nad12-Sigma Residual Welding Dictionary Candidate Spec

Status: `T129_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_ORBIT_QUOTIENT_NAD12_SIGMA_RESIDUAL_WELDING_DICTIONARY_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T128/N393` name one explicit future-only target object for the still-missing
semantic weld between:

1. the orbit-quotient provider-object carrier candidate lane (`T126/T127`),
2. the nad12-sigma residual pair-provider carrier target semantics (`T63/N328`),
3. the residual bridge/export-map object-support frontier (`N302` scope).

The next honest move after a target naming is not to pretend the weld is
already discharged.

It is weaker:

```text
export one explicit welding dictionary *candidate* form
that makes all identifications explicit as data (slot map, pair map, carrier map),
without claiming that the weld target (T128) is discharged.
```

`T129` specifies exactly such a candidate dictionary.

## Inputs

1. Orbit-quotient provider-object carrier candidate (`T126`):

   ```text
   ProviderObjectCarrier_pair^{cand,orbit}
   ```

   with explicit contraction parameters `(a,b)` and pair index set
   `PairIndex_v1 = {pair1, pair2}`.

2. Declared nad12 scaffold index carrier reused by the nad12-sigma lane (`T63`)
   through `C14/R11`:

   ```text
   Nad12Index_v1 := {0,1,...,11}
   ```

   (declared carrier only; no physical canonicalization implied).

3. Bridge-facing carrier-field type for the phasor reduction candidate (`T116`):

   ```text
   E_pair := {pair1: {paths:[(d,w)]}, pair2: {paths:[(d,w)]}}
   ```

## Candidate welding dictionary object

Export one explicit candidate dictionary:

```text
Pi_strict_provider_object_carrier_orbit_quotient_nad12_sigma_residual_welding_dictionary_candidate_v1
```

carrying as explicit data:

1. **Pair map (identity at this stage):**

   ```text
   W_pair_map_v1 : {pair1,pair2} -> {pair1,pair2}
   W_pair_map_v1(pairi) := pairi
   ```

2. **12-slot map (declared identification; no cyclic quotient):**

   ```text
   W_slot_map_v1 : Nad12Index_v1 -> Nad12Index_v1
   W_slot_map_v1(k) := k
   ```

   and one explicit finite orbit-depth truncation map (bridge-facing; no claim
   that the infinite orbit-depth index set is already canonicalized):

   ```text
   OrbitDepthIndex_v1 := {0,1,...,11}
   W_orbit_depth_to_nad12_slot_map_v1 : OrbitDepthIndex_v1 -> Nad12Index_v1
   W_orbit_depth_to_nad12_slot_map_v1(k) := k
   ```

   This matches the already exported declared carrier convention that the
   shared `12`-slot scaffold is indexed by `0..11` (`C14/R11/R8` scope), while
   keeping the physical canonicalization explicitly open (`QW-2191` remains in
   force).

3. **Carrier-field builder (bridge-facing; noncyclic):**

   Given `r_1 := |a|`, `r_2 := |b|`, define for each `i ∈ {1,2}` and
   `k ∈ Nad12Index_v1`:

   ```text
   d_{i,k} := k
   w_{i,k} := r_i^k / (sum_{j=0..11} r_i^j)
   ```

   and output:

   ```text
   E_i(pair) := { (d_{i,k}, w_{i,k}) }_{k=0..11}
   E_pair := { pair1: E_1(pair), pair2: E_2(pair) }.
   ```

This dictionary is intended as a *typed interface object* that makes the
welding attempt checkable:

- it explicitly identifies the shared 12-slot carrier,
- it explicitly defines the bridge-facing `E_pair` carrier field,
- it stays noncyclic (no theta input; no populated-instance input),
- it stays observer-free (no `K_obs`-indexed selection as primary source).

## Scope / non-claims

`T129` is **not**:

1. a discharge of the welding target `T128`,
2. a discharge of `T125` (provider-object carrier-layer target),
3. a discharge of `T63/N328` (nad12-sigma residual provider target),
4. a discharge of `N302` (bridge/export-map object support),
5. an export-map object,
6. an admissible `S_sel_int`,
7. a `QW-2191` discharge,
8. a ToE closure move.

It is only an explicit candidate dictionary form to be used by later welding
tests/probes.

## Hard limits

`T129` must not claim:

1. that the declared 12-slot carrier is physically canonical,
2. that the dictionary implies uniqueness/selector closure,
3. that any bridge/export-map object support exists above `N302`.
