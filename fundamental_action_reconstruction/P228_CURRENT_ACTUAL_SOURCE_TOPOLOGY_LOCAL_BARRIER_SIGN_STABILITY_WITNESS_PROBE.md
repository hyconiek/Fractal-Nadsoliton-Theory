# P228 Current Actual Source Topology Local Barrier Sign Stability Witness Probe

Status: `P228_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P228` tests whether the current repo really exports one actual source-side
local barrier-sign stability witness, while keeping the result:

1. below full barrier-protected sign discharge,
2. below full source-topology nontriviality discharge,
3. below basis-independent selector promotion,
4. below quotient-safe `QW-2191` resolution,
5. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Witnesses under test:

```text
epsilon_src_local_barrier_radius_v1 := delta_src_barrier_sign_margin_v1 / 2

chi_src_local_barrier_sign_stability_witness_v1 :
for all epsilon in R,
if |epsilon| <= epsilon_src_local_barrier_radius_v1,
then sign(cos(phi + epsilon)) = +1
```

## What P228 checks

`P228` checks only:

1. the packet is explicitly exported,
2. the already exported scalar sign component witness remains `+1`,
3. the local radius is positive,
4. the declared local radius stays strictly inside the positive-cosine core
   domain `(-pi/2, pi/2)`,
5. the local witness stays source-side,
6. observer remains downstream only,
7. the result remains below full barrier-protected sign discharge,
8. the result remains below full source-topology nontriviality,
9. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P228` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_BELOW_FULL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P228
```

This means:

1. the current repo exports one actual local barrier-sign stability witness for
   `tau_src_candidate_v1`,
2. the current repo exports one explicit positive-radius neighborhood on the
   declared core branch where the sign remains `+1`,
3. but it still does not export full barrier-protected sign discharge,
4. and it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P228` does not establish:

1. full barrier-protected sign of `tau_src_candidate_v1`,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. current selector closure,
6. current global `QW-2191` discharge,
7. ToE closure.
