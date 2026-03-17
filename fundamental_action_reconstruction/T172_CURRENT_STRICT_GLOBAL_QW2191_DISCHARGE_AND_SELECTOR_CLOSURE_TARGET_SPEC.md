# T172 Current Strict Global QW-2191 Discharge + Selector Closure Target Spec

Status: `T172_CURRENT_STRICT_GLOBAL_QW2191_DISCHARGE_AND_SELECTOR_CLOSURE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After the following strict exports now present in the repo:

1. **Global atlas + transition/gluing objects on `C_v1`** (`T170` discharged: `F469/N515`),
2. **Global selector state objects on `C_v1`**
   - projective/ray-level (`F470/N516`),
   - directed/sign-sensitive in an explicit **premise-tracked** scope (`T171` discharged via `F474/N524`, with fixing datum `T164` via `F473/N523`),
3. **Seed‑v1 strict-core internal selector-source lane advanced to global promotions**
   - admissibility clauses for the exported strict-core source object packetized (`N540–N545`),
   - admissible orientation export and local `B_sel/R_sel/O_sel` operators exported (`N546–N549`),
   - global promotions on `C_v1` exported (`N550–N553`) **at projector/section level with residual sign gauge tracked**,

the repo is now past the “infrastructure export” phase for strict selector continuation on `C_v1`.

The remaining honest strict frontier is:

```text
global strict selector closure + the corresponding QW-2191 uniqueness discipline
```

meaning:

1. define what *closure* means at the global level (projective vs directed),
2. prove the closure object is well-defined under the exported atlas/transition data,
3. prove the closure object yields a unique physical selector outcome under the chosen strict scope,
4. keep residual sign / fixing-datum dependence explicit (no hidden marked direction, no Aut(Z_12)-invariant overclaim),
5. keep `QW-2191` explicit and do not claim a global discharge unless a theorem-level uniqueness statement is actually exported.

This is a **target spec only**. It exports no new object.

## Scope

`T172` is scoped only to the global uniqueness/selector closure frontier under `QW-2191` discipline.

`T172` does **not** decide:

1. ToE closure,
2. any legacy→strict role transfer beyond what is explicitly exported,
3. any operator-level transition groupoid promotion beyond section/projector level (boundary `N512`).

## Target objects (what would count as a discharge)

Export, at minimum, a theorem-level package that explicitly pins down:

### A) Global selector closure object (declared scope)

One explicit global closure object on `C_v1`, for one of the two scopes:

1. **Projective closure (ray-level):**

   ```text
   SelectorClosure_global_C_v1_projective_strict_v1
   ```

   meaning: a closure notion that treats residual `Z2` sign as gauge at state level and is defined purely on projectors/rays.

2. **Directed closure (vector-level):**

   ```text
   SelectorClosure_global_C_v1_directed_strict_v1
   ```

   meaning: a closure notion that treats the selector state as directed/sign-sensitive in the declared scope, with any fixing datum (e.g. `T164`) tracked as an explicit premise.

In both cases, the closure object must be explicitly typed on the declared domain `C_v1` and must refer to exported atlas/transition/state objects rather than “implicit overlap” semantics.

### B) Global uniqueness / QW-2191 resolution statement (theorem-level)

Export one explicit theorem-level statement that makes precise what is meant by “global QW-2191 discharge” in the strict program:

- either a uniqueness statement in the chosen closure scope (projective or directed),
- or a strict boundary theorem explaining why uniqueness cannot be promoted from the currently exported data without adding new structure.

The theorem must keep the distinction explicit:

```text
QW-2191 (kernel-alone obstruction) is not "false";
it is either bypassed by an exported strict internal source
or it remains an open boundary.
```

## Acceptance tests (no false pass)

Any honest `T172`-class discharge must satisfy all of:

1. **Global well-definedness:** the closure object is defined using the exported global atlas + transition objects on `C_v1`,
   with explicit overlap/transport/cocycle discipline.
2. **Level discipline:** any cocycle/path-independence claim must state its level:
   - projector-level is acceptable as sign-gauge-safe,
   - vector-section level is acceptable only as a tracked convention layer (unless sign is exported as physical in a premise-tracked scope),
   - operator-level groupoid promotion is forbidden unless explicitly proven (`N512` boundary).
3. **Sign discipline:** residual `Z2` sign must be handled explicitly:
   - either proven gauge-irrelevant for the closure observable,
   - or fixed by an exported premise (e.g. `T164`), with the premise tracked and the non-`Aut(Z_12)`-invariance declared (`N462` boundary).
4. **No external selector smuggling:** no silent import of `QW-2192` (selector axiom) into strict core; if used, it must be labeled as axiom-augmented scope.
5. **Kernel-split safety:** no silent transfer of legacy physical-role claims onto strict objects.
6. **No implied ToE closure:** the export must keep ToE closure and global uniqueness scope-limits explicit.

## Current-state note (why `T172` is the next honest label)

The repo already exports:

- `T170` infrastructure: global atlas + transition/gluing objects on `C_v1`,
- `T171` post-projective infrastructure: a directed/sign-sensitive state layer in an explicit premise-tracked scope,
- seed‑v1 promotions on `C_v1` (`N550–N553`) but explicitly still marked:
  - `strict_core_selector_closure = false`,
  - `QW2191_discharge = false`.

So the honest next strict move is to either:

1. export a theorem-level `T172`-class closure + uniqueness package, **or**
2. export a theorem-level boundary showing why the remaining promotion cannot be done from the current exports without new structure.

## Hard limits

`T172` must not claim:

1. strict-core selector closure is already achieved,
2. global `QW-2191` discharge is already achieved,
3. ToE closure.

