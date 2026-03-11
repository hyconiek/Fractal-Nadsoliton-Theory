# F318 First Current Strict-Side Source-Seed Strict Sigma-Int Upgrade Candidate Instance Packet

Status: `F318_EXECUTED_FIRST_CURRENT_STRICT_SIDE_SOURCE_SEED_STRICT_SIGMA_INT_UPGRADE_CANDIDATE_INSTANCE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

The strict-side selector-seed lane historically fixed the first seed candidate:

```text
S_sel_int_candidate_seed_v0 := (QW-2206_local_topological_protection_layer_in_local_B_tilde_1_sector, sigma_int_candidate)
```

via `F35/F36`, where `sigma_int_candidate` was defined under hybrid FR support
(`B4/B5`).

On the updated repo state (`P390/P391`), the strict sigma-int lane now exports
an **actual** strict-core sigma-int source-upgrade datum:

```text
sigma_int_strict_derived_v1 ∈ {+1,-1}
```

via `F307/N418`, with explicit strict-side premise provenance and **no** silent
hybrid FR reuse.

The next honest strict-side move is therefore narrower than selector closure:

```text
export one strict-sigma-int-upgraded source-seed candidate instance
for the strict-side selector ingredient lane,
without changing strict core and without implying admissible S_sel_int.
```

`F318` executes exactly that move.

## Inputs reused

1. `F36`
   - seed-v0 construction instance (historical baseline),
2. `F307/N418`
   - strict sigma-int source-upgrade datum `sigma_int_strict_derived_v1`,
3. `QW-2206`
   - local topological protection layer label for the fixed local sector,
4. `A10`
   - anti-overclaim boundary.

## Packet result

`F318` exports one strict-sigma-int-upgraded seed candidate instance:

```text
S_sel_int_candidate_seed_v1 :=
(
  QW-2206_local_topological_protection_layer_in_local_B_tilde_1_sector,
  sigma_int_strict_derived_v1
).
```

## Exact meaning

This packet means only:

1. the sigma-int slot of the strict-side seed construction can now be filled by
   the exported strict-core datum `sigma_int_strict_derived_v1` (`F307/N418`),
2. the local-sector protection label remains the same fixed precursor label
   (`QW-2206`),
3. this is still only a seed candidate instance for future strict-side
   selector-ingredient work.

## Hard limits

`F318` must not claim:

1. that `sigma_int_candidate` is identified with `sigma_int_strict_derived_v1`,
2. admissible `S_sel_int`,
3. strict-core selector closure or `QW-2191` discharge,
4. strict-core theta export / target-slot population,
5. ToE closure.

