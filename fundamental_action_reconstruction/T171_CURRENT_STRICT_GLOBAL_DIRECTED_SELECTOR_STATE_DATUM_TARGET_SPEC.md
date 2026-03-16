# T171 Current Strict Global Directed Selector State Datum Target Spec (Post‑Projective Frontier)

Status: `T171_CURRENT_STRICT_GLOBAL_DIRECTED_SELECTOR_STATE_DATUM_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After:

1. `T170` discharge: strict core exports a **global selector atlas** and a **global selector transition/gluing object** on the declared strict domain `C_v1`
   (`F469`, packaged by `N515`),
2. `H39` discharge: strict core exports a **global projective selector state object** on `C_v1`
   (`F470`, packaged by `N516`),
3. sign‑gauge hygiene: residual `Z2` sign is packaged as gauge‑irrelevant for the exported downstream objects where it is provably irrelevant (`N502`),
   and extended to the exported global projective selector atlas/transition/state objects (`N519`, supported by probe `P474`),

the remaining strict “selector state” frontier is now sharply isolated:

```text
if one demands a sign-sensitive / directed selector state datum (distinguishing u from -u as physically inequivalent),
that datum is not exported in strict core today (audits H36/H37).
```

`T171` names that missing object class precisely, so the repo cannot accidentally “upgrade” a representative sign convention
into a strict physical claim (no false pass).

This is a **target spec only**. It exports no new object.

## Scope

`T171` is scoped only to the post‑projective frontier:

- export a **directed/sign-sensitive** selector state datum or observable on `C_v1` (lifting residual `Z2`),
  *compatible* with the already exported global **projective** selector state object.

`T171` does **not** decide:

1. strict-core selector closure / admissible `S_sel_int`,
2. global discharge of `QW-2191`,
3. ToE closure.

## Current-state boundary (why this is a real frontier)

On the current repo state:

1. strict core exports the selector state only at the **projective/ray level** (`F470/N516`),
2. strict core exports **no** sign-sensitive physical observable distinguishing `u` from `-u` on `pair1` (audit `H37`),
3. strict obstructions show a large tempting class cannot work:
   - even / `Aut(Z_12)`-invariant reference weight families cannot distinguish sign on the current exported `pair1` sine axis via `Σ_x w(x) u_1(x)` (`N517`, `N518`),
   - a repo scan finds no strict(-derived) weight-like reflection-breaking per-site arrays outside non-canonical marked-site `K_total` rows (`P472`).

Therefore any future strict attempt to lift residual sign must introduce **new structure** (explicit and tracked) or prove sign irrelevance for the specific downstream observable.

## Target objects (what would count as a discharge)

Export, at minimum, one of the following two strict-core objects, with explicit typing on the declared domain `C_v1`:

### Option A: Directed selector state object (vector-level lift)

Export a global directed selector state object:

```text
SelectorState_global_C_v1_directed_strict_v1
```

with intended meaning:

```text
a vector-level (sign-sensitive) lift of the projective selector state,
so that u and -u are not identified at state level in the declared scope.
```

### Option B: Sign-sensitive physical observable class (sign distinguishes states)

Export an explicit strict observable (or observable class) whose value flips under `u -> -u` **and**
is declared to have physical meaning in the strict scope, e.g.:

```text
S_dir_pair1_strict_v1(u_1) ∈ R  with  S_dir(-u_1) = -S_dir(u_1),
```

with an explicit statement of:

1. its domain of definition,
2. its provenance classification (strict-derived vs explicit strict-side premise),
3. and why it does **not** smuggle an untracked generator/orientation choice.

## Acceptance tests (no false pass)

Any honest discharge in the `T171` class must satisfy all of:

1. **Compatibility with the exported projective state:** the directed lift must descend to the already exported projective/ray object:
   ```text
   P(u_m) := |u_m><u_m|  matches the chart-local projectors of SelectorState_global_C_v1_projective_strict_v1 (F470).
   ```
2. **No hidden marked-direction slot:** the discharge must not silently introduce a `Z_12` generator/orientation choice.
   - If it depends on a marked generator/orientation, it must explicitly discharge the fixing-datum target (`T164`) or be labeled non-strict.
3. **Strict vs extension-lane discipline:** any sign-fixing observable based on a marked direction (e.g. `AX28`) remains `strict_extension_only`
   unless a strict fixing datum is exported.
4. **Obstacle awareness:** the discharge must explicitly address why it bypasses the known obstructions (`N517`, `N518`) rather than re-labeling an even/Aut-invariant weight family.
5. **Noncyclic + observer-free:** no `K_obs` indexing; no reusing cyclic L5/L12 loops under the same blocker-cut without a new provider class.
6. **Hard limits kept explicit:** no implied strict-core selector closure, no implied global `QW-2191` discharge, no implied ToE closure.

## Recommended next honest move (before any discharge attempt)

Before attempting a strict discharge of `T171`, the repo must explicitly choose one of two honest continuations:

1. **Projective-only continuation:** treat the exported global projective selector state as the strict physical state object for the declared closure stack,
   keeping residual sign as gauge/convention where proven irrelevant (`N502`, `N519`), and proceed with strict-only ToE closure tasks that depend only on projectors/spans.
2. **Directed continuation:** introduce and track a genuine symmetry-breaking/orientation source (e.g. by discharging `T164` or by an explicit domain lift such as a cover),
   then export a directed lift satisfying the above acceptance tests.

`T171` exists so either continuation stays honest and does not smuggle a false PASS.

## Hard limits

`T171` must not claim:

1. that a strict directed selector state datum is already exported,
2. that `u` vs `-u` are physically distinguished by strict core today,
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.

