# T150 Current Strict Sigma-Int `E_pair` Eps-Parameter Strict-Provenance Target Spec

Status: `T150_CURRENT_STRICT_SIGMA_INT_E_PAIR_EPS_PARAMETER_STRICT_PROVENANCE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `T117/F270/N382`, the strict sigma-int lane exports an internal-datum
driven **candidate** generator for a finite carrier field `E_pair`, but that
generator depends on a free amplitude parameter:

```text
eps ∈ [0,1]
```

After the strict sigma-int lane map export discharge (`P388/P391`), and after
the strict eps-source audit (`F315/N426`), the next honest missing ingredient
is no longer a bridge/export-map object.

It is:

```text
a strict-provenance source for the amplitude parameter eps used by the
sigma-int-driven E_pair generator weight law
```

`T150` names that missing ingredient as one explicit future-only target object
with explicit acceptance tests, so that no parameter choice (e.g. `eps = 1/2`)
can be silently read as strict-core theta supply.

## Scope

`T150` is scoped only to the eps-parameter used in:

```math
w_{i,k} := \frac{1 + \sigma_{int}\,\varepsilon\, b_{i,k}}{12}.
```

It does not decide:

1. strict-core `theta_1`, `theta_2` export,
2. actual bridge/export-map object support above the exported map object,
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.

## Target object

If the current repo cannot yet export an actual eps value object with strict
provenance, export one future-only target object:

```text
Delta_sigma_int_E_pair_eps_parameter_strict_provenance_target_v1
```

with the intended meaning:

```text
export one strict-provenance (strict-derived or strict-source-upgraded),
observer-free, noncyclic eps value object supplying the amplitude parameter eps
for the sigma-int-driven E_pair generator, so that the generator ceases to be a
free-parameter family with respect to eps
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Delta_sigma_int_E_pair_eps_parameter_strict_provenance_target_v1` must at
minimum provide:

1. **A dedicated eps value object:** an exported object/value
   `eps_sigma_int_E_pair_amplitude_strict_provenance_v1` with contract:
   - `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 ∈ [0,1]`.
2. **Explicit provenance classification:** the eps value object must be
   explicitly marked as one of:
   - `strict_derived`, with a declared derivation chain on a strict domain, or
   - `strict_source_upgraded`, via an explicit strict-side premise with explicit
     provenance,
   and it must not rely on silent reuse of unrelated “epsilon” symbols.
3. **Noncyclic contract:** no `theta_{1,2}` inputs and no populated basis-pair
   instance may be used as inputs (respects `N18`).
4. **Observer-free contract:** no `K_obs`-indexed selection may serve as the
   source of eps.
5. **Semantic separation / symbol-identification discipline:** if any existing
   exported `epsilon_*` radius/stability objects are cited, an explicit
   bridge/identification witness must be exported; otherwise they are treated
   as unrelated symbols.
6. **Selector discipline:** the discharge must explicitly state whether it:
   - keeps `QW-2191` open, or
   - adds a further strict-side selector/symmetry-breaking ingredient.
   In either case, no implied selector closure is permitted.

## Hard limits

`T150` must not claim:

1. discharge of `N302` / `N395` (actual object support above the export-map
   object),
2. actual strict-core `theta_1`, `theta_2`,
3. actual populated basis-pair instance,
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

