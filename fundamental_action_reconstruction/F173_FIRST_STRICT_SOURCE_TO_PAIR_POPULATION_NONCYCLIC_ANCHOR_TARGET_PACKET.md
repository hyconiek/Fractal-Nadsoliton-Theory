# F173 First Strict Source-To-Pair-Population Noncyclic Anchor Target Packet

Status: `F173_EXECUTED_FIRST_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Following `T26`, the next honest move is not to pretend that the missing
strict-side selector ingredient already exists.

It is narrower:

```text
freeze one explicit future-only strict-side solution target
which could in principle break both:
  - the official extension-lane freeze from N283
  - the sandbox noncyclic-anchor freeze from N18
without using observer as selector source
and without reactivating K_legacy_ont
```

`F173` executes exactly that move.

## Fixed future-only target

Reuse the already actual source-side packet family:

```text
Pi_sel_src_actual_witness_v1
Upsilon_sel_basis_actual_witness_v1
Phi_qw2191_safe_actual_witness_v1
T14_src_selector_declared_scope_actual_witness_v1
```

Freeze one explicit future-only target:

```text
V_strict_src_pair_population_anchor_target_v1 :
  tau_src_candidate_v1 -> strict_src_pair_population_noncyclic_anchor_target_v1
```

with abstract component packet:

```text
Gamma_strict_src_pair_population_anchor_components_target_v1 :=
(
  strict_source_orientation_seed_target_v1,
  pair_indexed_population_anchor_target_v1,
  strict_internal_orientation_provider_target_v1
)
```

## Meaning of the target packet

`strict_src_pair_population_noncyclic_anchor_target_v1` is intended only as:

1. one future source-side anchor target,
2. one pair-indexed noncyclic entry point for later populated basis-pair
   construction,
3. one future route below actual `E_orient`,
4. one proposed new strict-core ingredient class,
5. not a current discharge.

## Target properties

Freeze the following properties as part of the target packet:

1. `anchor_is_source_side := yes`

   Reason:
   the proposal reuses source-topology selector support and keeps observer
   downstream.

2. `anchor_is_observer_free := yes`

   Reason:
   `N234` blocks observer-side promotion to selector closure.

3. `anchor_is_Kobs_independent := yes`

   Reason:
   `H15/N13` keep `K_obs` as a distinct hypothesis-extension lane rather than
   an exported strict-core source identity.

4. `anchor_is_kernel_split_safe := yes`

   Reason:
   the proposal does not reactivate `K_legacy_ont` as live constructive
   kernel and does not claim bridge erasure.

5. `anchor_is_pair_indexed := yes`

   Reason:
   the blocker in sandbox `N18` is not purely class-level; it is exposed at
   the populated pair-instance layer.

6. `anchor_is_noncircular := target_yes_but_not_yet_discharged`

   Reason:
   the entire point of the proposal is to provide one anchor from outside the
   existing loop, but the anchor itself is still only a target.

7. `future_route_only := yes`

## Why this packet is sharper than the current freeze

`F173` does not merely repeat:

1. one more extension lift under `N283`,
2. one more sandbox loop expansion under `N18`.

Instead it isolates one narrower proposed solution class:

```text
source-side pair-population noncyclic anchor
```

This is stronger than saying only:

```text
some new ingredient is needed
```

because it fixes the proposed shape of that ingredient while still remaining
below actual discharge.

## Result

`F173` exports one explicit future-only strict-side solution target:

```text
V_strict_src_pair_population_anchor_target_v1 :
  tau_src_candidate_v1 -> strict_src_pair_population_noncyclic_anchor_target_v1
```

with declared properties:

1. source-side,
2. observer-free,
3. `K_obs`-independent,
4. kernel-split-safe,
5. pair-indexed,
6. intended to break the current noncyclic deficit,
7. still future-only,
8. no false pass.

## Hard limits

`F173` does not discharge:

1. actual strict-side selector ingredient,
2. actual theta inputs,
3. actual populated basis-pair instance,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. actual ToE closure.
