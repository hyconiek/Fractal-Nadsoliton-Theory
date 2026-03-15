# T167 Current Strict Canonical Local-Diagonal Opposite‑Pair Sums Sigma Strict‑Derived Value‑Instantiation Target Spec

Status: `T167_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_OPPOSITE_PAIR_SUMS_SIGMA_STRICT_DERIVED_VALUE_INSTANTIATION_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`T166` names the strict decision target:

```text
decide F2(d) for the canonical FIN D_local_residual profile (strict-derived),
either prove F2(d)=0 or prove F2(d)≠0, without hidden slots.
```

`N467/P426` reduce that decision to the six opposite‑pair sums on `n=12`:

```text
Sigma_psi0_psi6, Sigma_psi1_psi7, Sigma_psi2_psi8,
Sigma_psi3_psi9, Sigma_psi4_psi10, Sigma_psi5_psi11.
```

`P434` exports a reproducible evaluation harness which computes:

- `Re(F2)`, `Im(F2)`, `|F2|`,
- the induced `pair1` anisotropy signature, and
- (if nonzero) the canonical diagonalization angle `theta_*`,

*once numeric values for those six sums are available*.

On the current repo state:

1. `N472/P431` prove coefficient‑class underdetermination (no strict decision),
2. `P432` confirms no decision‑ready numeric instantiation is exported at the designated input location,
3. `N474/N475/N477` provide conditional Yukawa‑free reductions under constant‑vacuum stationarity + `vpsi_k≠0`,
   but still do **not** export strict‑derived vacuum/coupling value objects, hence still no strict numeric values for
   the six sums.

Therefore the next honest missing object class on the diagonal accelerator lane is:

```text
strict-derived numeric instantiation of the six opposite‑pair sums
for the canonical FIN D_local_residual profile,
with explicit provenance and without hidden selector slots.
```

`T167` names that missing strict value‑instantiation target sharply, so the repo cannot “accidentally” treat
toy/probe assignments (e.g. `P436`) as strict evidence.

Update (`2026-03-15`):
the repo now exports a strict sigma-six-sum instantiation consumable by `P434` via the strict constrained lift/provider
`F447` (theorem-level anchored by `N483`, audited by `P448`). The populated designated input is:

- `fundamental_action_reconstruction/generated/p434_input_sigma_opposite_pair_sum_values_candidate.json`

## Scope

`T167` is scoped only to:

1. the canonical FIN local diagonal residual sector (`R15`) on `n=12`,
2. the six opposite‑pair sums `Sigma_psi{k}_psi{k+6}` used by the `T166` decision target,
3. the diagonal/local `pair1` `O(2)` cut attempt (`N466`).

It does **not** decide:

1. strict-core theta export (`T159`),
2. strict eps/delta slot elimination (`T160/T161/T162`),
3. strict-core selector closure or global `QW-2191` discharge,
4. ToE closure.

## Target object

If the current repo cannot yet export strict‑derived numeric values for the six sums, export one explicit future‑only
target object:

```text
Delta_canonical_D_local_residual_sigma_opposite_pair_sums_strict_derived_value_instantiation_target_v1
```

with intended meaning:

```text
export a dedicated strict-derived value object instantiating:
  Sigma_psi0_psi6..Sigma_psi5_psi11
for the canonical FIN D_local_residual profile,
so P434 becomes computable and T166 can be discharged as a strict decision.
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Delta_canonical_D_local_residual_sigma_opposite_pair_sums_strict_derived_value_instantiation_target_v1`
must at minimum provide:

1. **Dedicated value object:** an exported value object (e.g. JSON) that provides numeric (finite real) values for all
   six keys:
   - `Sigma_psi0_psi6`,
   - `Sigma_psi1_psi7`,
   - `Sigma_psi2_psi8`,
   - `Sigma_psi3_psi9`,
   - `Sigma_psi4_psi10`,
   - `Sigma_psi5_psi11`.
2. **Canonical profile identification:** a strict statement pinning down that these sums correspond to:
   ```text
   Sigma_psi{k}_psi{k+6} := d_k + d_{k+6}
   where d_k = (D_local_residual)_{kk} for the canonical FIN diagonal decomposition (R15).
   ```
3. **Strict-derived provenance (not premise-only):** the value object must be explicitly classified as `strict_derived`
   with an explicit derivation/selection chain on a declared strict domain. It must not be only a
   `strict_source_upgraded` premise (or a toy/probe instantiation).
4. **No hidden slot:** the derivation/selection chain must not introduce an untracked “choose a site”, “choose a
   generator”, “choose a vacuum branch”, or “choose a representative” selector slot. If branch selection is required,
   it must be exported as a separate strict ingredient with explicit provenance and scope.
5. **Noncyclic + observer-free contract:** the derivation must not use:
   - `theta_{1,2}` outputs or any populated basis-pair instance as input (`N18`),
   - `K_obs`‑indexed selection as a primary source.
6. **Evaluation compatibility:** the six values must be consumable by `P434` so the repo can compute an explicit
   decision of `F2(d)` and (if nonzero) the induced canonical `pair1` axis data (`N468`).
7. **No false pass discipline:** discharging this target alone must not imply strict-core theta export, strict-core
   selector closure, global `QW-2191` discharge, or ToE closure unless separately proved.

## Relation to the current evaluation harness

`P434` reads its numeric inputs from:

- `fundamental_action_reconstruction/generated/p434_input_sigma_opposite_pair_sum_values_candidate.json`

The current file intentionally contains `null` values and warns against false promotion.

Discharging `T167` may either:

1. export a dedicated strict-derived sigma value object and then populate the `P434` input file from it, or
2. replace the designated input location with a strict-derived value object and update the probe accordingly,

but in either case provenance and scope hygiene must remain explicit.

If a future strict provider exports the remaining upstream dependency set (e.g. a strict-derived vacuum vector
`vpsi[0..11]` and self-coupling families `g4[0..11]`, `g6[0..11]`), then:

- `P437` provides an explicit evaluation harness implementing the `N477` `K_total`-numeric Yukawa-free formula to compute
  `Sigma_psi0_psi6..Sigma_psi5_psi11` reproducibly, before feeding those values into `P434`.

## Hard limits

`T167` must not claim:

1. that the six opposite‑pair sums are already strict-derived today,
2. that any toy/probe numeric instantiation is strict evidence,
3. discharge of `T166` / `QW-2191`,
4. ToE closure.
