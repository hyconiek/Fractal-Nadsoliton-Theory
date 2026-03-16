# T170 Current Strict Global Selector Atlas + Transition Object Target Spec

Status: `T170_CURRENT_STRICT_GLOBAL_SELECTOR_ATLAS_AND_TRANSITION_OBJECT_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `B3/B8` and the sigma-int residual-datum bridge discharge (`N491`), the strict program is now in a clarified state:

1. the continuous `QW-2191` `O(2)` ambiguity can be cut down to residual `Z2` on declared lanes by exported strict ingredients
   (diagonal/local lane: `N487` + `F453/N492`; Shannon element‑order reference lane: `N480/N488/N496` + `F454`),
2. residual sign can be frozen as a tracked gauge/convention layer for the **currently exported downstream objects** where sign is provably gauge‑irrelevant
   (`N501`, packaged by `N502`),
3. nevertheless, strict core still exports **no global selector atlas** and **no global selector transition/gluing object** on the full strict domain `C_v1`
   (audits `H41` and `H40`).

`T170` names the next strict missing object class precisely, so we do not confuse:

- lane-scoped projector-level atlas ingredients on the `n=12` Fourier carrier,
- with a global atlas / global transition structure on `C_v1`,
- and do not silently promote section-level cocycle data into an operator-level transition groupoid (boundary `N512`).

This is a strict target spec only. It exports no new object.

## Scope

`T170` is scoped only to the global uniqueness/selector frontier under `QW-2191` discipline:

- `H41`: global selector atlas / overlap-domain declaration / cocycle-level gluing data on `C_v1`,
- `H40`: global selector transition/gluing object on `C_v1` sufficient to support a future global selector closure attempt.

`T170` does **not** decide:

1. strict-core selector closure / admissible `S_sel_int`,
2. global discharge of `QW-2191`,
3. ToE closure.

## Target objects

Export, at minimum, the following two strict objects (typed, with explicit declared domains):

1. **Global atlas (H41-class):**

   ```text
   SelectorAtlas_global_C_v1_strict_v1
   ```

   containing:
   - a declared family of local selector charts `{U_i -> Chart_i}` whose domains `U_i ⊂ C_v1` are explicit as sets/conditions,
   - explicit overlap-domain declarations `U_i ∩ U_j` for each relevant pair of charts,
   - explicit cocycle/gluing data declarations (what equations are claimed, and on which overlap domains).

2. **Global transition/gluing object (H40-class):**

   ```text
   SelectorTransition_global_C_v1_strict_v1
   ```

   containing:
   - explicit transition maps/operators on overlaps (e.g. `O_ij`),
   - an explicit statement of the level at which cocycle is claimed:
     - projector/ray level **or**
     - vector-section level **or**
     - full operator-level matrix identity,
   - with boundaries kept explicit (in particular: `N512` forbids false operator-level groupoid promotion when only section-level cocycle is available).

## Acceptance tests (discharge conditions)

An honest export in the `T170` class must satisfy all of:

1. **Domain explicitness:** chart domains and overlap domains are declared as subsets/conditions in `C_v1`,
   not merely “overlap of exported artifacts”.
2. **Typed transitions:** each transition object is exported as a typed object acting between declared chart carriers.
3. **Cocycle discipline:** any cocycle/path-independence claim must specify its level and be backed by a theorem/probe:
   - projector-level cocycle is acceptable as sign-gauge-safe,
   - vector-section cocycle is acceptable only as a tracked convention layer (no physical sign claim),
   - operator-level cocycle/groupoid requires an explicit operator-level proof and must not be inferred from section-level data (`N512` boundary).
4. **`QW-2191` compatibility:** the export must keep `QW-2191` explicit and must state what extra strict ingredient breaks kernel-alone non-uniqueness
   (e.g. diagonal/local profile, Shannon element‑order reference, sigma-int corridor transition data), without silently importing `QW-2192`.
5. **No hidden selector slot:** the atlas/transition data must not smuggle a marked-direction/generator choice equivalent to a selector knob,
   unless that premise is exported and labeled explicitly as non-strict.
6. **No false pass:** no implied claim of strict-core selector closure / admissible `S_sel_int`, no implied global `QW-2191` discharge, no implied ToE closure.

## Current-state notes (why this is the honest frontier)

The repo already exports lane-scoped ingredients on the `n=12` Fourier carrier:

- two-chart stub and overlap declaration (`F463`),
- projector-level atlases with cocycle audits up to five charts (`F464/F465/F466`, packaged `N508/N509/N510`),
- oriented transport lift at vector level as a tracked convention layer (`F467`, packaged `N511`),
- and a strict boundary theorem that operator-level cocycle identities do not hold on the full carrier (`P471`, packaged `N512`).

But `H40/H41` remain globally open because those exports are lane-scoped and do not provide a global cover/overlap declaration on `C_v1`.

So the next honest strict move is not to repeat more lane-scoped atlas enlargement under the same blocker-cut, but to either:

1. export an explicit global atlas/transition structure on `C_v1` meeting the above acceptance tests, **or**
2. export a strict non-bridge/boundary theorem showing why such global objects cannot be obtained from the currently exported strict inputs without new structure.

## Hard limits

`T170` must not claim:

1. a global selector atlas is already exported,
2. a global selector transition/gluing object is already exported,
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.

