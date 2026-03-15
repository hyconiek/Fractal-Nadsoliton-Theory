# P113 Current Strict-Core Internal Selector Source Derivation Discharge Probe

Status: `P113_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_STATE`
As of: `2026-03-15`

## Goal

After `F28`, the last remaining higher-order strict-core question is:

```text
does the current repo already support an explicit strict-core internal selector
source derivation discharge?
```

## Result

Update (`2026-03-15`): the historical package-level non-discharge packaging is no longer supported as-is on the current repo state,
because downstream operator-stage reachability has changed (`P2` now reaches an `A_1(pair1)` operator stage via `F456`).

Therefore `P113` must be re-evaluated before making any package-level negative closure claim about strict-core internal selector source derivation.

```text
P113_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_STATE
```

## What was checked

`P113` checks only the package-level strict-core source question. It asks
whether the current repo simultaneously exports:

1. any derived generic strict-core internal orientation datum,
2. any current strict-core `psi0` selector-source discharge,
3. any current strict-core FR-route internal selector-source discharge,
4. any current strict-core `sigma_int` bridge discharge strong enough to act as
   a selector-source carrier,
5. and any strict-core downstream reachability beyond those blocked source
   routes.

## Why it fails

On the current repo state, at least one prerequisite used by the old packaging is now false:

- `P2` no longer blocks on downstream operator export; it reaches an operator stage in declared scope (`F456`).

Therefore the probe cannot be used to justify a package-level negative closure conclusion without re-deriving the package-level frontier.

## Real reduction after P113

The higher-order frontier is no longer safely describable as:

```text
the current repo already supports a package-level non-discharge conclusion
```

It is now:

```text
the package-level strict-core internal selector source frontier must be re-evaluated after the sigma-int operator-stage export
```

## What P113 does not claim

`P113` does not claim:

- that no future strict-core selector source can ever be added,
- strict-core selector closure,
- `QW-2191` discharge,
- theory-level acceptance of selector requirement,
- ToE closure.

## Recommended next move

The correct next move is now:

1. keep `QW-2191` explicit (no implied global discharge),
2. decide (and document) whether the newly exported sigma-int operator-stage object is being treated only as a minimal downstream witness
   or as part of an admissible strict-core internal selector source package in a declared scope,
3. rerun/re-derive the package-level frontier (and only then attempt any theorem-level packaging such as `N124`/`N126`).
