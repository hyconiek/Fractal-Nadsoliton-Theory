# T18 Sandbox Strict-Core Theta-Supply Population Loop Incompatibility Boundary Scope

Status: `T18_SANDBOX_STRICT_CORE_THETA_SUPPLY_POPULATION_LOOP_INCOMPATIBILITY_BOUNDARY_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F16/P16/N16` and `F17/P17/N17`, the sandbox no longer has a vague
strict-core source gap.

It has one explicit dependency loop:

```text
actual strict-core theta-source supply
    waits for
actual populated basis-pair instance

actual populated basis-pair instance
    waits for
actual theta inputs
```

The next direct question is:

```text
is this now an honest incompatibility boundary on current sandbox state and
current strict-core inputs?
```

## Boundary type

`T18` is **not** a permanent impossibility theorem.

It is only a current-state incompatibility boundary with the following scope:

1. current strict-core inputs,
2. current sandbox artifact set,
3. current blocker-cut exposed by `N16` and `N17`,
4. current noncyclic guardrail from `S2` / `QW-2381/QW-2382/QW-2383`.

## Intended move

`T18` will state only that:

1. this positive route has reached a fixed point on the present sandbox state,
2. further positive lifting inside the same loop would be cyclic repetition
   under the same blocker-cut,
3. therefore the honest next move is boundary freeze or blocker-cut change,
   not another same-loop positive lift.

## Hard limits

`T18` must not claim:

1. permanent impossibility in all future repo states,
2. impossibility under a genuinely new provider class,
3. impossibility under a different blocker-cut,
4. actual strict-core theta supply,
5. actual populated basis-pair instance,
6. actual provider emission,
7. actual strict-core selector closure,
8. actual ToE closure.
