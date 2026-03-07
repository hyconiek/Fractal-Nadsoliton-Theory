# N5 Current Strict-Core psi0-Route Obstruction Theorem

Status: `N5_DISCHARGED_CURRENT_STRICT_CORE_PSI0_ROUTE_OBSTRUCTION_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N4`, the next higher-value theorem is not a global impossibility claim
over every conceivable future strict-core route.

The right target is narrower:

- the current strict-core `psi0` route itself.

`N5` asks whether the current strict core can close the selector problem through
the `psi0` lane without importing an extra selector or symmetry-breaking datum.

## Theorem

### Informal statement

Within the current strict core:

1. kernel alone is already theorem-level obstructed from giving full physical
   uniqueness because of continuous residual `O(2)` freedom (`QW-2191`),
2. no internal orientation datum exists in strict core to replace the missing
   symmetry-breaking input (`B2`),
3. the current `psi0` lane exports at most a deterministic angle candidate plus
   local coordinate embedding, but no chart-independent selector object,
4. the only currently computed nontrivial split remains extension-only and
   anchor-imported.

Therefore the current strict-core `psi0` route is obstructed: it cannot close
selector generation without adding extra symmetry-breaking structure beyond the
current strict core.

### Formal statement

```text
N5_CurrentStrictCore_Psi0Route_Obstruction_Theorem

Let S_psi0_current denote the current strict-core psi0 route consisting of:
  kernel invariants,
  local pair charts,
  local psi0 embeddings,
  and current strict-core selector-track objects on pair1.

If:
  (i) QW-2191 proves that kernel alone leaves a continuous O(2) family and
      therefore requires extra symmetry breaking for full physical uniqueness,
  (ii) B2 proves that no internal orientation datum is currently derived in the
       strict core to discharge that obstruction,
  (iii) H30/H31/H33/H34/H35/H36/H37/H38 show that psi0 supplies at most a
        candidate angle and local projective embedding structure, but no strict
        selector target, no chart-independent reduction, no physical axis,
        no directed orientation, and no sign-sensitive selector state,
  (iv) H42 and P1 show that the first nontrivial pair1 split appears only as an
       anchor-imported extension-lane effect,

then S_psi0_current cannot discharge selector closure in strict core.

Hence any successful psi0-based selector closure requires at least one extra
symmetry-breaking structure not currently present in the strict core.
```

## Proof

### Step 1. Kernel alone is already obstructed

From `QW-2191`:

- degenerate kernel eigenspaces induce continuous `O(2)` rotation freedom,
- the rotated family preserves kernel-subspace invariance and Lie-closure
  audits,
- full uniqueness from kernel alone is obstructed,
- explicit symmetry breaking is required.

So the strict-core `psi0` lane cannot rely on kernel structure alone to choose a
unique physical selector.

### Step 2. No internal strict-core orientation datum is available

From `B2`:

- `strict_internal_selector_derivations_found = 0`,
- `strict internal orientation datum = not_found_in_strict_core`,
- `kernel invariant selecting one O(2) point = not_found_in_strict_core`.

So the strict core does not currently contain the missing datum that could
replace the external selector.

### Step 3. The current `psi0` lane does not elevate beyond local chart structure

From `H30`:

- `psi0` is deterministic from kernel invariants,
- but it is not exported as strict-core selector datum.

From `H31` and `H33`:

- `psi0 -> pair1` exists only as a coordinate embedding,
- `pair1` is only a deterministic local chart, not a proven selector target.

From `H34`:

- no basis-covariant / target-independent selector reduction is exported.

So the current `psi0` lane does not rise from local chart structure to
chart-independent strict-core selector reduction.

### Step 4. No strict orientation object exists on `pair1`

From `H35`, `H36`, `H37`, and `H38`:

- no strict physical axis is selected on `pair1`,
- no directed orientation is selected on that axis,
- no sign-sensitive state distinguishes `u` from `-u`,
- at most a local projective/ray-level representative exists.

So even the local `pair1` chart does not carry the strict orientation object
needed for selector closure.

### Step 5. Nontrivial splitting appears only outside strict core

From `H42`:

- bare `c`-based retardation is selector-trivial,
- nontrivial splitting requires imported `psi0`,
- `psi0 + c` is not strict core.

From `P1`:

- the first computed nontrivial block is classified as
  `ANCHOR_IMPORTED_SPLIT`,
- the lane is `hypothesis_extension_only`,
- `strict_core_promotion = false`.

So the repo does exhibit a nontrivial split, but only outside the current
strict core.

### Conclusion

The current strict-core `psi0` route inherits the `QW-2191` obstruction, lacks
an internal orientation datum, lacks a chart-independent selector reduction,
lacks a strict orientation object on `pair1`, and reaches a nontrivial split
only on an imported extension lane.

Therefore:

```text
the current strict-core psi0 route is obstructed and cannot discharge selector
closure.
```

The route-specific corollary follows:

```text
any successful psi0-based selector closure requires at least one extra
symmetry-breaking structure not currently present in the strict core.
```

## What is discharged

`N5` discharges:

- a theorem-level obstruction statement for the current strict-core `psi0`
  route,
- a theorem-level necessity statement that extra symmetry-breaking structure is
  required for psi0-based selector closure.

## What remains open

`N5` does not discharge:

- a global impossibility theorem for every possible future strict-core route,
- `QW-2191` itself beyond its existing scope,
- `T2`,
- `T12`,
- full selector closure,
- full ToE closure.

## Acceptance boundary

`N5` is acceptable only if all of the following stay explicit:

1. the theorem is route-specific to the current strict-core `psi0` lane,
2. it does not forbid future new strict-core structures,
3. it does not promote the extension lane into strict core,
4. it does not claim that a specific future selector mechanism is already known,
5. it does not claim full ToE closure.

## Recommended next move

Only two serious routes remain:

1. add one concrete strict-core symmetry-breaking structure and re-test the
   `psi0` lane end-to-end,
2. or escalate from this route-specific obstruction to a broader impossibility
   theorem only if a new argument avoids the old `T12` meta-blocker.
