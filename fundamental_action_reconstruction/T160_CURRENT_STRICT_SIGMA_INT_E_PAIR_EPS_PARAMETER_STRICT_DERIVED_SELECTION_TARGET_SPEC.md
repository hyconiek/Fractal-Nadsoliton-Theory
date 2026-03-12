# T160 Current Strict Sigma-Int `E_pair` Eps-Parameter Strict-Derived Selection Target Spec

Status: `T160_CURRENT_STRICT_SIGMA_INT_E_PAIR_EPS_PARAMETER_STRICT_DERIVED_SELECTION_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After the eps-provenance target (`T150`) and its discharge (`F317/N428`), the
strict sigma-int lane exports one dedicated eps value object:

```text
eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2
```

with explicit provenance classification:

```text
strict_source_upgraded (explicit premise; not strict-derived).
```

However, `P407/N441` prove that on the current strict sigma-int → theta
candidate pipeline the computed theta-pair depends on admissible eps choices.
Therefore eps remains a real selector slot and cannot be silently upgraded to
strict-core by citing only one premise-based eps value.

`T159` names a strict-core upgrade target for an `O(2)`-cut selector ingredient.
After `N443`, the invariance-based slot-elimination route is blocked on the
current exports. The only remaining honest route for slot removal in strict
core is:

```text
export a genuinely strict-derived (not premise-only) eps selection law/value object,
or change the construction class so the slot does not exist.
```

`T160` names the missing **strict-derived eps selection** ingredient sharply,
as a future-only target with explicit acceptance tests.

## Extension-lane continuation note (explicitly non-strict)

On the current repo state, the strict-core eps selection target remains open:
no strict-derived eps law/value object is exported (`N446`, `P408`).

If one insists on proceeding *today* with a single reproducible sigma-int → theta
representative without false pass, the repo explicitly separates that move into
`strict_extension_only` scope:

1. `AX21` freezes `eps := 1/2` as a declared symmetry-breaking premise (premise-based; not strict-derived).
2. `AX22` packages a publication-ready strict-extension summary of that lane.

This does **not** discharge `T160`. It only records the explicit extension-scope
premise needed for a reproducible representative while keeping strict-core eps
selection open.

## Scope

`T160` is scoped only to the amplitude parameter `eps` of the strict sigma-int
driven `E_pair` generator weight law:

```math
w_{i,k} := \frac{1 + \sigma_{int}^{in}\,\varepsilon\, b_{i,k}}{12},
\qquad \varepsilon = eps \in [0,1]
```

as used by the strict sigma-int → theta candidate lane (`T117` and successors).

It does **not** decide:

1. strict-core `theta_1`, `theta_2` export,
2. discharge of post-`T148` object-support targets (`N302/N395/T130`),
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.

## Target object

If the repo cannot yet export a strict-derived eps selection ingredient, export
one explicit future-only target object:

```text
Delta_sigma_int_E_pair_eps_parameter_strict_derived_selection_target_v1
```

with intended meaning:

```text
export one dedicated eps value object with STRICT_DERIVED provenance
(not premise-only), together with an explicit strict derivation/selection chain
on a declared strict domain, so the strict sigma-int → theta construction no
longer depends on a free eps selector slot in strict core
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Delta_sigma_int_E_pair_eps_parameter_strict_derived_selection_target_v1`
must at minimum provide:

1. **Dedicated eps value object:** an exported value object
   `eps_sigma_int_E_pair_amplitude_strict_derived_v1` with contract:
   - `eps_sigma_int_E_pair_amplitude_strict_derived_v1 ∈ [0,1]`.
2. **Strict-derived provenance (no premise-only):** the object must be explicitly
   classified as `strict_derived` with an explicit derivation/selection chain on a
   declared strict domain. It must not be only a `strict_source_upgraded` premise
   (e.g. `eps := 1/2`) and must not rely on symbol reuse of unrelated
   `epsilon_*` radius/stability objects without a bridge theorem.
3. **Slot-closure discipline:** the derivation/selection chain must not smuggle
   in any free selector slot. In particular:
   - it must not take `theta_{1,2}` outputs or any populated basis-pair instance as input (`N18`),
   - if it depends on `delta_d`, then `delta_d` must itself be supplied by a strict-derived
     value object (not a premise-only corridor convention), and the dependence must be explicit.
4. **Observer-free contract:** no `K_obs`-indexed selection as a primary source.
5. **Compatibility with `T159`:** the discharge must state how this strict-derived eps object
   removes the eps selector slot required by the `T159` strict-core upgrade acceptance test (2B).
6. **No false pass discipline:** the discharge must not imply:
   - strict-core theta export,
   - discharge of `N302/N395/T130`,
   - admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
   unless separately proved.

## Relation to existing exports (current-state discipline)

On the current repo state:

1. eps is exported only as strict provenance (premise-based):
   `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2` (`F317/N428`),
   which does **not** discharge `T160`,
2. the most literal `Z2` parity-balance / “zero-charge point” constraints on the `T117` weight law do not
   derive `eps = 1/2` (they either do not constrain eps or they force `eps = 0`), packaged as `N446`,
3. no strict exported theorem currently derives `eps = 1/2` from any typed strict “charge parity balance”
   law/objective (`P408` audit).

## Hard limits

`T160` must not claim:

1. that eps is already strict-derived on the current strict sigma-int lane,
2. that premise-based eps provenance (`F317/N428`) constitutes strict derivation,
3. strict-core theta export or strict-core selector closure,
4. ToE closure.
