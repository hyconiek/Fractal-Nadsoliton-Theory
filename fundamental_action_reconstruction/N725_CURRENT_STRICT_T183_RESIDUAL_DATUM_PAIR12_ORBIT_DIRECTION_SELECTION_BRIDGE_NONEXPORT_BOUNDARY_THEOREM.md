# N725 Current Strict `T183` Residual-Datum Pair1/Pair2 Orbit-Direction Selection Bridge Nonexport Boundary Theorem

Status: `N725_CURRENT_STRICT_T183_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_BOUNDARY_THEOREM_NO_FALSE_PASS`

## Question

After `P728/N724` reduces the surviving positive corridor to the strict
residual-datum sublane `pair1,pair2`, does the current repo already export the
next finer bridge object:

`ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1`?

That target means:

1. use the already exported residual-datum carrier on `pair1,pair2`,
2. identify whether the remaining `pair1/pair2` ambiguity is already localized
   on current exports as one precise branch type,
3. and ask whether the current repo already selects one of those two remaining
   branches.

## Packaged input state

This theorem packages the strongest honest current answer from:

1. `P728` (source-side reduction to `pair1,pair2`),
2. `F301` (actual residual-datum source-side support carrier),
3. `F461` and the exported `alpha12` transition-angle object,
4. `F462` (pair1/pair2 projector-level glue).

## Theorem-level conclusion

On the current repo state:

1. `P728` already reduces the surviving positive corridor to the
   boundary-shielded sublane `pair1,pair2`;
2. inside that already selected sublane, `F301` now localizes the remaining
   split exactly as two opposite orbit-index branches on the same source-side
   carrier:
   - `pair1` carries the branch
     `q_{1,k} = \delta_k`,
   - `pair2` carries the branch
     `q_{2,k} = \delta_{-k}`;
3. the exported carrier vectors confirm that these two branches differ only by
   index inversion `k \mapsto -k`, not by support size or amplitude profile;
4. the exported `pair1↔pair2` lane transport and projector glue remain exact,
   so current lane data relate those two branches faithfully but still do **not**
   select one of them as the distinguished source-side branch;
5. therefore the remaining `pair1/pair2` ambiguity is now sharply localized on
   current exports as:
   - a positive-index orbit branch,
   - versus a negative-index orbit branch;
6. however, the current repo still exports **no** source-side rule selecting
   between those two orbit-direction branches.

So the current strict state narrows again:

```text
the residual pair1/pair2 ambiguity is no longer merely a chart-label split;
it is now localized as opposite orbit-index directions on the same
already-exported residual-datum carrier,
but one branch is still not selected
```

## Sharp next consequence

The next honest strict move should therefore no longer ask whether the current
source-side residual carrier prefers:

- `pair1,pair2` versus `pair3,pair5`,
- or even `pair1` versus `pair2` as bare chart labels.

Those reductions are already packaged.

The missing step is now finer:
export a residual-datum source-side rule distinguishing

- the positive-index orbit branch `\delta_k`
from
- the negative-index orbit branch `\delta_{-k}`

inside the already selected `pair1/pair2` sublane.

## Hard limits

This theorem does **not** claim:

1. `T183` discharge,
2. a unique chart-seed selection,
3. promotion of orbit-index direction into a physical orientation datum,
4. identification of `tau_src` with the current selector carrier,
5. kernel-alone/global `QW-2191` discharge,
6. or ToE closure.
