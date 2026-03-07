# N4 Current-Repo psi0 Strict-Core Nonderivation Theorem

Status: `N4_DISCHARGED_CURRENT_REPO_PSI0_STRICT_CORE_NONDERIVATION_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P1`, the narrow question is no longer whether an extension lane can
numerically split `pair1`, but whether the current repository state already
contains a strict-core derivation turning `psi0` into a selector source without
import.

`N4` does not claim a future-proof global impossibility theorem.
It discharges a narrower theorem about the current exported repo state only.

## Theorem

### Informal statement

Within the current repository state:

1. `psi0` is deterministic from kernel invariants,
2. but no exported strict-core packet upgrades `psi0` from a kernel-invariant
   angle candidate to a selector source on `pair1`,
3. and every currently computed nontrivial selector split on `pair1` remains
   extension-only and anchor-imported.

Therefore the current repository state does not contain a strict-core
derivation of `psi0` as selector source, and any current computable selector
split requires an extension lane.

### Formal statement

```text
N4_CurrentRepo_Psi0_StrictCore_Nonderivation_Theorem

Let R_current be the current exported repository state on 2026-03-07.

If in R_current:
  (i) psi0 is only certified as a deterministic kernel-invariant candidate,
  (ii) psi0 -> pair1 is only certified as a coordinate embedding,
  (iii) no basis-covariant / target-independent selector reduction is exported,
  (iv) no strict physical axis, directed orientation, or sign-sensitive
       selector state is exported on pair1,
  (v) every computed nontrivial pair1 selector split is classified as
      extension-only and anchor-imported,

then R_current contains no strict-core derivation turning psi0 into a selector
source on pair1.

Hence any currently computable selector split on pair1 is extension-necessary
relative to R_current.
```

## Proof

### Step 1. `psi0` is deterministic but not exported as strict-core selector datum

From `H30`:

- `deterministic_from_kernel_invariants = true`,
- `strict_core_selector_export = false`,
- `residual_O2_discharge = false`.

So `psi0` is stronger than a free heuristic parameter, but it is not yet a
strict-core selector datum.

### Step 2. The current `psi0 -> pair1` move is only a coordinate embedding

From `H31`:

- `coordinate_level_embedding_present = true`,
- `strict_core_selector_reduction_present = false`,
- `pair1_as_selector_target_proven = false`,
- `psi0_equals_theta1_proven = false`.

So the current repo exports only

```text
u_psi0_pair1 = cos(psi0)c_1 + sin(psi0)s_1
```

as a legal chart embedding, not as a strict-core selector reduction.

### Step 3. The embedding is not elevated beyond chart dependence

From `H34`:

- strict core contains local chart embeddings for `psi0`,
- but no `basis-covariance / target-independence` theorem elevates those
  embeddings beyond chart dependence.

So the repo still lacks a strict-core argument that the `pair1` reduction is a
physical selector reduction rather than a chart choice.

### Step 4. No strict orientation object is exported on `pair1`

From `H35`, `H36`, and `H37`:

- no strict physical axis selection is exported on `pair1`,
- no strict directed orientation is exported on that axis,
- no sign-sensitive selector state or observable distinguishes `u` from `-u`.

So even after embedding `psi0` into `pair1`, the repo still lacks the strict
orientation object needed for a selector source.

### Step 5. The only computed nontrivial split is imported and extension-only

From `H42`:

- bare `c`-based retardation on `pair1` is selector-trivial,
- the nontrivial split appears only when `psi0` is imported into anisotropic
  path data,
- `psi0 + c` is not classified as strict core.

From `P1`:

- the first computed nontrivial block on `pair1` is classified as
  `ANCHOR_IMPORTED_SPLIT`,
- the lane is `hypothesis_extension_only`,
- `strict_core_promotion = false`.

Therefore the repo does compute a nontrivial operator split, but only after
importing the anchor through an extension lane.

### Conclusion

The current exported repo state contains:

- a deterministic kernel-invariant `psi0` candidate,
- a legal coordinate embedding into `pair1`,
- but no strict-core selector reduction, no chart-independent selector target,
  no strict orientation object, and no strict-core operator map producing the
  computed nontrivial split.

Therefore:

```text
R_current contains no strict-core derivation turning psi0 into a selector
source on pair1.
```

And the operational corollary follows:

```text
Any currently computable selector split on pair1 is extension-necessary
relative to the current repository state.
```

## What is discharged

`N4` discharges:

- a theorem-level current-repo statement that `psi0` is not currently exported
  as a strict-core selector source,
- a theorem-level operational corollary that current computable selector splits
  on `pair1` require extension.

## What remains open

`N4` does not discharge:

- a future-proof impossibility theorem for all possible strict-core extensions,
- `QW-2191`,
- `T2`,
- `T12`,
- full selector closure,
- full ToE closure.

## Acceptance boundary

`N4` is acceptable only if all of the following stay explicit:

1. the theorem quantifies only over the current exported repository state,
2. it does not claim that future strict-core work cannot add a new source,
3. it does not promote `ANCHOR_IMPORTED_SPLIT` to strict core,
4. it does not convert `psi0` into `theta_1`,
5. it does not claim theorem-level/full-closure ToE completion.

## Recommended next move

Only two high-value routes remain:

1. add one real strict-core source object that upgrades `psi0` from candidate
   angle to selector source and then export the strict-core map to `A_1(pair1)`,
2. or prove a stronger impossibility theorem beyond current-repo scope.

Any further audit-ladder growth without one of those two moves should be
treated as documentation, not closure progress.
