# T168 Current Strict Canonical Vacuum + Self‑Coupling Value‑Provider Target Spec (for Diagonal Sigma Six‑Sums)

Status: `T168_CURRENT_STRICT_CANONICAL_VACUUM_AND_SELF_COUPLING_VALUE_PROVIDER_FOR_DIAGONAL_SIGMA_SIX_SUMS_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`T166` (diagonal/local `pair1` `O(2)`‑cut decision) is currently blocked because the canonical diagonal residual mode‑2
defect

$$
F_2(d)
$$

cannot be evaluated without strict-derived information that fixes the **canonical** diagonal residual profile
$d_k=(D_{\mathrm{local\_residual}})_{kk}$.

`N467/P426` reduce `F_2(d)` to six opposite‑pair sums
`Sigma_psi0_psi6..Sigma_psi5_psi11`, and `T167` names the missing strict-derived **numeric instantiation** target for
those six sums so `P434` becomes computable.

Independently, `N474/N475/N477` show a strict conditional reduction:

```text
under constant-vacuum stationarity (canonical EoM) and vpsi_k ≠ 0,
the canonical diagonal residual entries (and therefore the six opposite-pair sums)
admit an explicit Yukawa-free, K_total-numeric form (R14) with explicit m0^2 floor (R15).
```

`P437` exports an evaluation harness for the `N477` formula, but it is **NOT COMPUTABLE** today because the repo does
not export strict-derived numeric inputs for:

```text
vpsi[0..11], g4[0..11], g6[0..11]
```

Therefore the next honest missing object **beneath** `T167` (if the repo chooses the `P437` route rather than a direct
sigma‑six‑sum value export) is:

```text
strict-derived canonical vacuum + local self-coupling value provider
for vpsi[0..11], g4[0..11], g6[0..11],
with explicit provenance and without hidden selector slots,
so P437 can compute Sigma_psi{k}_psi{k+6} and feed P434.
```

`T168` names that missing provider class sharply, to prevent “toy/probe” instantiations from being verbally promoted
into strict evidence.

## Scope

`T168` is scoped only to the diagonal/local accelerator lane:

1. canonical constant-vacuum specialization (`QW-2165`),
2. Yukawa-free diagonal rewrite under stationarity (`N474/N475`),
3. strict kernel specialization + scalar floor (`R14/R15`),
4. computation of the six opposite‑pair sums via `N477/P437`,
5. feeding those sums into `P434` to decide `F2(d)` (`T166`).

It does **not** decide:

1. strict-core theta export (`T159`),
2. strict eps/delta slot elimination (`T160/T161/T162`),
3. strict-core selector closure or global `QW-2191` discharge,
4. ToE closure.

## Current strict scalar-side status (not yet a `T168` discharge)

The repo **does** already export a strict scalar-side vacuum-floor chain:

1. `QW-2122` exports a scalar psi-potential broken-branch package including:
   - `lambda_psi_strict`,
   - `rho_star_sq`,
   - the broken-branch diagonal floor `diag_floor`.
2. `QW-2123` exports a strict branch-selection rule (broken branch required) under `lambda_min(K_total)<0`.
3. `QW-2124` promotes scalar vacuum closure to a branch-resolved strict PASS.

These scalar objects are sufficient for the `R15` host floor embedding, but they still do **not** export the canonical
per-site value provider required by `T168`, because the lift/mapping into the canonical per-site families
`vpsi[0..11], g4[0..11], g6[0..11]` (with no hidden selector slot) is not exported today.

`P442` audits this frontier explicitly and confirms that the scalar-side chain is present while the per-site provider is
missing.

Additionally, `N478` (with the toy audit `P443`) closes one common misreading:
even after fixing the strict scalar outputs (e.g. `rho_star_sq`, `m0^2`), the scalar data do not canonically lift to
the canonical per-site arrays needed by `T168`, and do not decide the diagonal defect `F2(d)` without an additional
strict mapping/selector ingredient (or an independent symmetry proof forcing `F2=0`).

After `F446/N480` discharge a direction‑free Shannon selector ingredient (`T165`), the remaining strict missing object
beneath this lane is now best stated as an explicit **constrained lift** target:
`T169` names the required strict mapping/selector chain from `QW-2122` scalar closure into the canonical per‑site arrays
needed by `T168` (existence + uniqueness, no hidden slots).

Update (`2026-03-15`):
the repo now exports one explicit strict constrained lift/provider (`F447`, theorem-level anchored by `N483`, audited by
`P448`) which supplies a `T168`-consumable per-site value object `(vpsi,g4,g6)` and auto-populates the designated
`P437` harness input file.

## Target object

If the repo cannot yet export strict-derived numeric values for the needed upstream dependency set, export one explicit
future‑only target object:

```text
Delta_canonical_vacuum_and_local_self_couplings_strict_derived_value_provider_target_v1
```

with intended meaning:

```text
export a dedicated strict-derived value object instantiating:
  vpsi[0..11], g4[0..11], g6[0..11]
for the canonical FIN constant vacuum/self-coupling data,
so P437 becomes computable (N477) and can produce the six Sigma_psi{k}_psi{k+6} values required by T167/P434.
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Delta_canonical_vacuum_and_local_self_couplings_strict_derived_value_provider_target_v1`
must at minimum provide:

1. **Dedicated value object:** an exported value object (e.g. JSON) that provides numeric (finite real) values for:
   - `vpsi[0..11]`,
   - `g4[0..11]`,
   - `g6[0..11]`.
2. **Canonical identification:** a strict statement pinning down that these values correspond to the **canonical**
   constant vacuum and local self-coupling families appearing in the canonical EoM/Hessian exports (`QW-2165/2166`)
   and used by the canonical diagonal decomposition (`R15`).
3. **Strict-derived provenance (not premise-only):** the value object must be explicitly classified as `strict_derived`
   with an explicit derivation/selection chain on a declared strict domain; not only a toy/probe assignment.
4. **No hidden slot:** the derivation/selection chain must not introduce an untracked “choose a site”, “choose a
   generator”, “choose a vacuum branch”, or “choose a representative” selector slot.
   - If a branch selection is required, it must be exported as a separate strict ingredient with explicit provenance
     and scope, and must be referenced by the value provider.
5. **Stationarity compatibility:** if the discharge intends to use the `N474/N475/N477` stationarity reductions, it must
   explicitly state whether the provided vacuum satisfies the nonzero premises `vpsi_k ≠ 0` needed by those theorems,
   or else provide the alternative handling for zero entries (no silent division).
6. **Evaluation compatibility:** the exported values must be consumable by `P437` so that the computed six
   `Sigma_psi{k}_psi{k+6}` can be copied into the `P434` input object and `T166` can be decided without false
   promotion.
7. **Noncyclic + observer-free contract:** the derivation must not use:
   - `theta_{1,2}` outputs or any populated basis-pair instance as input (`N18`),
   - `K_obs`‑indexed selection as a primary source.
8. **No false pass discipline:** discharging this target alone must not imply strict-core theta export, strict-core
   selector closure, global `QW-2191` discharge, or ToE closure unless separately proved.

## Hard limits

`T168` must not claim:

1. that `vpsi/g4/g6` are already strict-derived today,
2. that any toy/probe instantiation is strict evidence,
3. discharge of `T167/T166` or `QW-2191`,
4. ToE closure.
