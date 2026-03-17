# T174 Current Strict Global Oriented Transition (Edge Sign‑Lift) Target Spec

Status: `T174_CURRENT_STRICT_GLOBAL_ORIENTED_TRANSITION_EDGE_SIGN_LIFT_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After `T170` discharged the global atlas/transition objects on `C_v1` (`F469/N515`) and `T173` established that a rooted directed state
representative is edgewise compatible with the global transition object only **up to sign** (`P686/N686`) with a **global sign coherence
obstruction** under fixed axis-only (`α mod π`) transition representatives (`P687/N687`), the next honest continuation is:

```text
export an explicit oriented (α mod 2π) edge sign‑lift of the exported global transition object,
in an explicit strict_convention scope, so that a directed representative section can be transported without sign flips on every edge.
```

This target is **not** a strict-core uniqueness discharge. It is a tracked *convention layer* that makes the required additional data explicit.

## Scope

`T174` is scoped to:

1. the already exported strict global transition object class on `C_v1` (`F469`),
2. an explicit **edgewise** oriented lift datum (choice of `±` per overlap edge),
3. a convention-scoped directed representative section (e.g. `F684`) whose transport coherence is audited under the lifted transitions.

`T174` is **below**:

- strict-core physical sign datum claims,
- `Aut(Z_12)`-invariant canonicity claims (`N462`),
- kernel-alone/global `QW-2191` discharge,
- operator-level transition groupoid promotion (`N512` boundary).

## Target objects (what counts as a discharge)

### A) Global oriented transition edge sign‑lift object (strict_convention)

Export one object of the form:

```text
SelectorTransition_global_C_v1_oriented_mod_2pi_edge_sign_lift_*_strict_convention_v1
```

containing, for each overlap edge `pairi_to_pairj`, an explicit sign choice
`s_ij ∈ {±1}` interpreted as:

```text
O_ij^(oriented) := s_ij * O_ij^(axis-only)
```

so that the transported directed representative satisfies:

```text
O_ij^(oriented) u_i ≈ u_j   (no sign flips)
```

### B) Coherence audit (probe-level)

Export one probe that:

1. loads the lifted transition object and the chosen directed representative section,
2. checks **all 10 edges** transport the directed vectors without sign flips,
3. writes an explicit `PASS_*` summary without any physical sign claim.

### C) Theorem-level packaging (boundary-safe discharge)

Export one theorem-level package that records:

1. the lift exists and is explicitly tracked as convention data,
2. it resolves the edgewise sign-flip obstruction **only** by adding oriented edge data,
3. no strict-core physical sign datum is implied.

## Acceptance tests (no false pass)

Any `T174` discharge must satisfy:

1. **Explicit data:** edge sign choices are exported as data, not smuggled into “canonical” operators.
2. **Scope marking:** must be labeled as `strict_convention` (or stricter), not strict-core physics.
3. **Audit:** must pass a full-edge coherence probe for the declared directed representative.
4. **Hard limits preserved:** must keep `QW2191_kernel_alone_discharge=false`, no `Aut(Z_12)`-invariant sign canonicity claim, and no operator-level groupoid claim.

## Hard limits

`T174` must not claim:

1. strict-core selector closure,
2. strict physical sign‑sensitive orientation datum,
3. `Aut(Z_12)`-invariant canonicity,
4. kernel-alone/global `QW-2191` discharge,
5. ToE closure.

