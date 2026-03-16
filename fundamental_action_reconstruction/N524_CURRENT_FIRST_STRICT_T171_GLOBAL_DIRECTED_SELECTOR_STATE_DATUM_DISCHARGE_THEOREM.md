# N524 Current First Strict `T171` Global Directed Selector State Datum Discharge Theorem (No False‑PASS)

Status: `N524_DISCHARGED_CURRENT_FIRST_STRICT_T171_GLOBAL_DIRECTED_SELECTOR_STATE_DATUM_DISCHARGE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

`T171` isolates the post‑projective strict frontier:

> export a **directed/sign-sensitive** selector state datum (or a sign-sensitive physical observable class) on `C_v1`,
> compatible with the already exported global **projective** selector state object, and without silently smuggling a `Z_12` generator/orientation choice.

This theorem packages the narrow claim that:

1. `T164` is discharged (premise‑based fixing datum exported; `F473/N523`),
2. and, using that fixing datum, `F474` exports:
   - a sign‑sensitive directed observable on `pair1`, and
   - a global directed vector‑level selector state object on `C_v1` descending to the exported projective state,

thereby discharging `T171` (no false pass).

## Strict‑admissible inputs reused

1. `F470/N516`
   - global projective selector state object on `C_v1`:
     `SelectorState_global_C_v1_projective_strict_v1`.
2. `F473/N523`
   - premise‑based strict provenance fixing datum (T164 discharge):
     `Kappa_Z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1`.
3. `F474`
   - exported directed sign‑distinction observable:
     `S_dir_pair1_strict_v1`,
   - exported global directed selector state object:
     `SelectorState_global_C_v1_directed_strict_v1`.
4. `N462`
   - boundary: no `Aut(Z_12)`‑invariant canonical generator/orientation choice; directed continuation must cite an explicit fixing datum.
5. `T171`
   - acceptance tests (compatibility + no hidden marked-direction slot).

## Theorem (T171 is discharged, premise‑based)

From `F473/N523`, an explicit fixing datum is exported as premise‑based strict provenance, reducing the admissible `Aut(Z_12)` ambiguity in a tracked way. This satisfies `T171` acceptance test (2) (“no hidden marked direction slot”): any directed object may cite this datum rather than smuggling a generator choice.

From `F474`, the repo exports:

1. `S_dir_pair1_strict_v1`, a sign-sensitive scalar observable on `pair1` with:
   - explicit definition `S_dir(u) = Σ_x w_dir(x) u(x)`,
   - explicit directed weight `w_dir(x)=exp(-alpha_geo*x)` tied to `alpha_geo_strict_derived_v1`,
   - explicit dependence on the exported fixing datum (`T164`), and
   - oddness under sign: `S_dir(-u) = -S_dir(u)`.
2. `SelectorState_global_C_v1_directed_strict_v1`, a vector-level directed lift of the already exported projective state on `C_v1`, constructed by applying the `S_dir` sign-fix to an exported oriented vector section on `{pair1..pair5}`.

Moreover, the directed state **descends** to the already exported projective state:

```text
P(u_m) := |u_m><u_m|  matches the chart-local projectors of SelectorState_global_C_v1_projective_strict_v1,
```

since a global sign flip does not change rank-one projectors/spans.

Therefore `T171` is discharged in the strict no-false-pass sense:

- an explicit directed/sign-sensitive observable and a directed global selector state datum are exported,
- compatibility with the exported projective state is maintained,
- the generator/orientation dependence is tracked via an explicit exported fixing datum (premise‑based; consistent with `N462`),
- no claim of selector closure, global `QW-2191` discharge, or ToE closure is made. ∎

## What `N524` does not claim

`N524` does not claim:

1. `Aut(Z_12)`‑invariant canonicity (ruled out by `N462`; the fixing datum is premise‑based),
2. strict‑core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.

