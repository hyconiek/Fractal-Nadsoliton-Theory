# F36 First Source-Seed Candidate Construction Instance Packet

Status: `F36_EXECUTED_FIRST_SOURCE_SEED_CANDIDATE_CONSTRUCTION_INSTANCE_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N132`, the next honest constructive question is no longer:

```text
which route should be used?
```

That is already fixed. The next question is now:

```text
what is the narrowest first candidate construction instance on the fixed route
QW-2206 local topological protection + sigma_int_candidate -> future S_sel_int?
```

## First candidate construction instance

The narrowest first candidate construction instance is:

```text
S_sel_int_candidate_seed_v0
```

defined only as the packaged local precursor pair

```text
S_sel_int_candidate_seed_v0 :=
(
  QW-2206_local_topological_protection_layer_in_local_B_tilde_1_sector,
  sigma_int_candidate
)
```

with:

- local sector support fixed by `QW-2206`,
- internal binary datum fixed by `B4`,
- local stability support fixed by `B5`.

## Why this instance is forced

The instance is forced by the current repo state:

1. `N132` already reduces the next constructive move to one precursor route,
2. that route has only one canonical current internal binary candidate:
   `sigma_int_candidate`,
3. the route is anchored on one canonical current strict-side local protection
   layer: `QW-2206`,
4. therefore the first honest construction attempt is one packaged seed
   instance built from exactly those two ingredients and nothing more.

## What F36 does count as

`F36` counts only as:

- the first explicit candidate construction instance for a future attempted
  source seed,
- a strict-side construction target packet,
- a frozen minimal carrier for the next future move.

## What F36 does not claim

`F36` does not claim:

- that `S_sel_int_candidate_seed_v0` is already admissible `S_sel_int`,
- that it is already a strict-core internal selector source object,
- that `E_orient` is already exported,
- that `B_sel`, `R_sel`, or `O_sel` are already reachable,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. test whether the current repo already reduces the next constructive move to
   this one explicit candidate construction instance,
2. and if so, freeze that instance as the only honest first future
   construction attempt before any admissibility upgrade is claimed.
