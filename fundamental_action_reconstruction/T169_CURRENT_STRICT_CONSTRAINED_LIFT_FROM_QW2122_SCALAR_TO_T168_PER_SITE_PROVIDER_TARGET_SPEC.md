# T169 Current Strict Constrained Lift Target Spec: `QW-2122` Scalar Closure → `T168` Per‑Site Provider (No False‑PASS)

Status: `T169_CURRENT_STRICT_CONSTRAINED_LIFT_FROM_QW2122_SCALAR_TO_T168_PER_SITE_PROVIDER_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`T168` is the real strict bottleneck under the diagonal/local accelerator lane (`T166`):
to compute the six opposite‑pair sums (`T167`) via `N477/P437`, the repo needs strict‑derived numeric arrays:

```text
vpsi[0..11], g4[0..11], g6[0..11].
```

`QW-2122/2123/2124` export only **scalar** vacuum/self‑coupling data (e.g. `rho_star_sq`, `lambda_psi_strict`, broken‑branch
floor), and `N478` proves that those scalar values do **not** canonically lift to the per‑site families on current strict
exports.

After `F446/N480`, the repo now has an actual strict Shannon‑typed selector ingredient that is **direction‑free** on `Z_12`
(element‑order reference; `N479`) and yields a theorem‑level `pair1` `O(2)->Z2` cut.

`T169` names the next missing strict object:

```text
a strict constrained lift operator / value-provider derivation chain
that canonically maps strict scalar closure data + strict selector ingredients
into the canonical per‑site arrays required by T168,
with theorem-level existence + uniqueness (up to an explicitly declared finite residual),
and with no hidden selector slots.
```

This is a target spec only. It does **not** claim the lift already exists.

Update (`2026-03-15`):
the repo now exports an explicit strict constrained lift/provider instance as:

- `F447` (packet),
- theorem-level existence + uniqueness anchor: `N483`,
- generated value object:
  `fundamental_action_reconstruction/generated/f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json`,
- provenance audit: `P448` (summary marks `theorem_level_pass=true`).

## Strict inputs (minimum)

An intended discharge of `T169` must explicitly reuse only strict‑admissible inputs such as:

1. `QW-2122/2123/2124` scalar closure outputs and branch rule (broken branch),
2. strict kernel specialization `K_total` (`R14`) and scalar floor embedding (`R15`),
3. the strict selector ingredient exported in `F446` (element‑order reference; no marked direction).

## Target object

Export one strict‑derived value provider:

```text
Delta_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1
```

with content consumable by `P437`:

```text
{ vpsi: [12 reals], g4: [12 reals], g6: [12 reals] }
```

and with explicit provenance fields tying it to the canonical constant‑vacuum and local self‑coupling families appearing in
`QW-2165/2166` (no silent reinterpretation).

## Acceptance tests (discharge criteria)

An **actual** discharge of `Delta_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1` must provide:

1. **Declared domain (candidate family):** a precise admissible set of per‑site arrays `(vpsi,g4,g6)` compatible with the
   strict scalar inputs (at minimum: a norm constraint using `rho_star_sq`), and explicit handling of the `vpsi_k≠0`
   premise if the discharge intends to feed `N477/P437`.
   - Note: since `N477/P437` uses ratios `vpsi_j/vpsi_k`, a discharge should also state whether any quantitative
     nondegeneracy condition is assumed (e.g. a lower bound on `|vpsi_k|`), or else provide an alternative formulation
     that avoids silent division instabilities.
2. **Explicit selection principle:** a strict objective/variational rule (or a strict constraint system) whose solution
   selects a unique member of that domain.
3. **No hidden slots:** the selection rule must not introduce untracked knobs of the form:
   - choose an origin/site,
   - choose a generator/direction,
   - choose a representative of a symmetry orbit,
   - choose a solver initialization/tolerance that changes the result.
4. **Theorem‑level existence + uniqueness:** a proof that the selected solution exists and is unique, modulo an explicitly
   declared finite residual (e.g. unavoidable `Z2` sign ambiguity).
5. **Deterministic computability:** a reproducible procedure (gate/algorithm) for producing the numeric arrays, with
   stable output schema.
6. **Strict provenance classification:** the provider must be explicitly classified as `strict_derived` with a derivation
   chain referencing only strict inputs and previously discharged strict theorems/packets. No “extension lane” promotion.
7. **Pipeline compatibility:** the exported values must be consumable by `P437` to compute the six opposite‑pair sums,
   which then feed `P434` to decide `F2(d)` (`T166`).

## Hard limits (no false‑PASS)

`T169` must not claim:

1. discharge of `T168` unless the above acceptance tests are met,
2. that the scalar closure alone already fixes per‑site arrays (`N478` forbids this),
3. discharge of `T167/T166` or `QW-2191`,
4. ToE closure.
