# T37 Current Residual Datum Sigma-Int Bridge Export Map Object Support Incompatibility Boundary Spec

Status: `T37_CURRENT_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_INCOMPATIBILITY_BOUNDARY_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N301`, the residual-datum / `sigma_int_candidate` third-provider route
already has all currently honest positive structure reachable on its
object-facing layer:

1. actual bridge-map target support through `N299`,
2. exact map-layer nonexport boundary through `N300`,
3. one explicit future-only bridge/export-map object target through `N301`.

So the strongest honest next question is no longer:

```text
can the same currently exported material be lifted once more positively
to actual bridge/export-map object support?
```

It is narrower:

```text
has the route now reached a current-state incompatibility boundary
between future-only bridge/export-map object targeting
and actual bridge/export-map object support?
```

`T37` does not decide the whole route in principle.

It does something narrower:

1. writes one packet-ready incompatibility-boundary spec for the missing
   object-support layer,
2. keeps the boundary scoped only to current repo state and current
   route-local exported material,
3. keeps open the possibility of a future genuinely new bridge/export-map
   object, one different object-support carrier, or one different blocker-cut.

## Formal target

```text
T37_CurrentResidualDatumSigmaIntBridgeExportMapObjectSupport_IncompatibilityBoundary

Assume:
  A1. `N297` still keeps the residual-datum / sigma-int branch as one
      future-only third-provider-class route;
  A2. `N299` already exports one actual support packet for the
      bridge/export-map target layer;
  A3. `N300` already freezes the exact bridge/export-map layer itself as
      nonexported on the current repo state;
  A4. `N301` already exports one future-only target object for the missing
      bridge/export map;
  A5. `P4/P5` still keep the exact missing strict-core object localized as
      `sigma_int_candidate -> residual orientation datum`;
  A6. `C40/C41/C42/C43/C44/C45/C46` export at most carrier grammar,
      template content, admission, and one minimal persisted template file,
      but not one actual bridge/export-map object support witness;
  A7. no actual bridge/export map, actual theta source, or actual component-2
      support witness is exported on the current repo state.

Then:
  C1. the route admits one current-state incompatibility-boundary theorem
      between future-only bridge/export-map object target
      and actual bridge/export-map object support;
  C2. another positive lift using only the same currently exported material
      is not the honest next move for this object layer;
  C3. this boundary remains weaker than impossibility in principle and may be
      reopened by one genuinely new bridge/export-map object, one genuinely
      new object-support carrier, or one different blocker-cut.
```

## Meaning of the theorem

If later discharged, `T37` would establish only:

1. the third-provider route is stronger than target-only on the map layer,
2. the route is also stronger than abstract missing-object language because
   the object is now sharply localized,
3. the current repo still stops at future-only object targeting,
4. the exact missing upward layer is actual object support rather than only
   target naming,
5. another same-material positive lift would overclaim the current repo state.

It would not establish:

1. actual bridge/export-map discharge,
2. actual bridge/export-map object support,
3. actual theta source,
4. actual component-2 support,
5. actual `theta_1`, `theta_2`,
6. actual populated basis-pair instance,
7. actual `E_orient`,
8. admissible `S_sel_int`,
9. strict-core selector closure,
10. ToE closure,
11. impossibility of all future third-provider routes.

## Acceptance skeleton

This spec is acceptable only if all of the following stay explicit:

1. the boundary is current-state only,
2. `N301` remains valid and is not weakened,
3. future-only object targeting is not silently relabeled as actual
   object-support,
4. template/carrier grammar is not silently relabeled as an actual
   bridge/export-map object,
5. no map export, theta export, or component-2 support is claimed,
6. no closure claim is introduced.

## Recommended next move

The correct next move is:

1. package one actual incompatibility-boundary packet for the missing
   bridge/export-map object-support layer of the residual-datum route,
2. test whether the current repo really exports that packet,
3. stop the same-material positive lifting there unless one genuinely new
   bridge/export-map object, one genuinely new object-support carrier, or one
   different blocker-cut is added.
