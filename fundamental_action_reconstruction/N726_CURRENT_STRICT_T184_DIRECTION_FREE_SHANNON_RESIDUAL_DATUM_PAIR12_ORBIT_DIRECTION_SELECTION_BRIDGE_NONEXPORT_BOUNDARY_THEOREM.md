# N726 Current Strict `T184` Direction-Free Shannon Residual-Datum Pair1/Pair2 Orbit-Direction Selection Bridge Nonexport Boundary Theorem

Status: `N726_CURRENT_STRICT_T184_DIRECTION_FREE_SHANNON_RESIDUAL_DATUM_PAIR12_ORBIT_DIRECTION_SELECTION_BRIDGE_NONEXPORT_BOUNDARY_THEOREM_NO_FALSE_PASS`

## Question

After `P729/N725` localizes the surviving `pair1/pair2` ambiguity as opposite
residual-datum orbit-direction branches `\delta_k` versus `\delta_{-k}`, does
the strongest already-exported strict-side Shannon lane now close that gap?

Concretely, does the current repo already export:

`DirectionFreeShannonResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1`?

That target means:

1. take the already-exported strict Shannon element-order reference lane
   (`F446/N480/N489`) on `pair1/pair2`,
2. evaluate it against the already-exported residual-datum carrier localized by
   `P729`,
3. and ask whether that current direction-free Shannon lane actually selects one
   of the two surviving orbit-direction branches.

## Packaged input state

This theorem packages the strongest honest current answer from:

1. `P729` (the surviving split is exactly `\delta_k` versus `\delta_{-k}`),
2. `F301` (actual residual-datum source-side support carrier),
3. `F446/N479` via the exported reference object `r_ord`,
4. `F451/N489` via the exported slot-free strict Shannon theta-pair source on
   `pair1/pair2`.

## Theorem-level conclusion

On the current repo state:

1. the strict Shannon element-order reference lane already exports a genuine
   slot-free `pair1/pair2` `O(2)\to Z_2` cut via the current theta-pair source;
2. however, that lane is built from the direction-free reference
   `r_ord(x)\propto e^{-\alpha_{\mathrm{geo}}\operatorname{ord}_{Z_{12}}(x)}`;
3. because `-1\in \mathrm{Aut}(Z_{12})`, both `\operatorname{ord}_{Z_{12}}` and
   therefore `r_ord` are invariant under orbit-direction inversion
   `k\mapsto -k`;
4. `P729/F301` show that the surviving pair1/pair2 residual-datum branches are
   exactly such inversion partners:
   - `pair1` carries `\delta_k`,
   - `pair2` carries `\delta_{-k}`,
   with identical carrier weights under inversion;
5. therefore the induced source-side Shannon scores on those two branches are
   equal:
   - the `\operatorname{ord}` expectation is the same on both branches,
   - the cross-entropy to `r_ord` is the same on both branches;
6. so the current direction-free Shannon lane does **not** pick one surviving
   branch over the other.

Hence the strongest honest current conclusion is:

```text
the current strict Shannon ord-reference lane already cuts pair1/pair2
from O(2) to residual Z2 at the pair-plane level,
but it still does not upgrade to a source-side selector between
δ_k and δ_{-k} on the surviving residual-datum carrier
```

Therefore the repo still does **not** export `T184`.

## Sharp next consequence

The next honest strict move should therefore no longer ask whether one more
reuse of the current direction-free Shannon element-order reference lane will
close the surviving `pair1/pair2` branch split.

That test is now packaged and negative.

If a strict source-side selector between `\delta_k` and `\delta_{-k}` is still
required, it must come from:

- a genuinely inversion-sensitive source-side rule on the residual-datum
  carrier,
- or another exported provider class not erased by `k\mapsto -k`,

and not from the present direction-free `ord_Z12/r_ord` lane alone.

## Hard limits

This theorem does **not** claim:

1. `T184` discharge,
2. promotion of orbit-index direction into a physical orientation datum,
3. that the current pair1/pair2 `O(2)\to Z_2` cut selects a unique source-side
   branch,
4. identification of `tau_src` with the current selector carrier,
5. kernel-alone/global `QW-2191` discharge,
6. or ToE closure.
