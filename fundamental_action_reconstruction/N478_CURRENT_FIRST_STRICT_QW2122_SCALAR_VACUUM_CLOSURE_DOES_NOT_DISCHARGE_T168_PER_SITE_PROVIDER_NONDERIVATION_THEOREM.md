# N478 Current First Strict QW‑2122 Scalar Vacuum Closure Does Not Discharge `T168` Canonical Per‑Site Provider (Nonderivation) Theorem

Status: `N478_DISCHARGED_CURRENT_FIRST_STRICT_QW2122_SCALAR_VACUUM_CLOSURE_DOES_NOT_DISCHARGE_T168_PER_SITE_PROVIDER_NONDERIVATION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`T166` (diagonal/local `pair1` `O(2)` cut decision) reduces to deciding the six opposite‑pair sums
`Sigma_psi0_psi6..Sigma_psi5_psi11` (`N467/P426/T167`).

`N474/N475/N477` show that **under constant‑vacuum stationarity and `vpsi_k≠0`** those six sums admit an explicit
Yukawa‑free `K_total`‑numeric form, but still depend on the remaining unknown per‑site data:

```text
vpsi[0..11], g4[0..11], g6[0..11].
```

`T168` names the missing strict-derived value-provider class for those arrays.

After `QW-2122/2123/2124` export a strict scalar vacuum closure chain (giving scalar values such as `rho_star_sq`,
`lambda_psi_strict`, `diag_floor`), a recurring temptation remains:

```text
maybe the strict scalar vacuum closure already (essentially) fixes the canonical per-site arrays needed by T168
```

`N478` closes that temptation honestly:

```text
the strict scalar vacuum closure does not canonically lift to the canonical per-site families (vpsi,g4,g6) on current exports,
so it does not discharge T168 and does not decide the diagonal/local mode-2 defect F2(d).
```

This is a nonderivation statement: it does **not** claim a strict-derived per-site vacuum, does **not** decide `F2(d)`,
and does **not** discharge `QW-2191`.

## Strict-admissible evidence reused

1. `QW-2122/2123/2124`
   - strict scalar vacuum closure chain exports scalar values on a broken branch:
     `rho_star_sq`, `lambda_psi_strict`, `diag_floor (=m0^2)`.
2. `R14/R15`
   - strict kernel channel specialization (`K_total`) and strict numeric diagonal floor `m0^2`.
3. `N477`
   - under stationarity + nonzero premises, the Yukawa-free diagonal residual and the opposite-pair sums are explicit
     `K_total`-numeric expressions in `(vpsi,g4,g6)` and `m0^2`.
4. `N476`
   - method template: prove underdetermination via explicit toy stationary witnesses (no strict promotion).

## Theorem (scalar closure does not discharge `T168`)

Even after fixing the strict scalar outputs from `QW-2122` (in particular `rho_star_sq` and `m0^2`), the repo does not
export any strict-derived canonical value-provider instantiating the per-site families:

```text
vpsi[0..11], g4[0..11], g6[0..11].
```

Moreover, the scalar data cannot by itself decide the diagonal/local mode-2 defect because:

1. `rho_star_sq` fixes only a **scalar norm constraint** of the form `Σ_i vpsi_i^2 = rho_star_sq`,
2. the `N477` diagonal/sigma expressions depend on **per-site ratios** `vpsi_j / vpsi_k` (and per-site `g4_k, g6_k`),
3. therefore multiple per-site instantiations compatible with the same scalar norm can yield different sigma six-sums,
   hence different `F2(d)`.

### Explicit stationary toy witnesses (same strict scalar inputs, different `F2(d)`)

Using the strict kernel channel and floor (`R14/R15`) and the strict scalar value `rho_star_sq` from `QW-2122`,
`P443` constructs two explicit toy instantiations which share the same scalar norm
`Σ_i vpsi_i^2 = rho_star_sq` and satisfy constant-vacuum stationarity by per-site solving of `m2_psi{i}`,
but yield different `F2(d)` after applying the `N477` formula:

1. a uniform `vpsi` representative gives `F2(d)=0`,
2. a two-site perturbed `vpsi` (same norm) gives `F2(d)≠0`.

Therefore, strict scalar closure does not fix the canonical per-site arrays and does not decide the diagonal defect.

## Audit probe

`P443` provides a reproducible toy-level numerical audit on the **actual** strict kernel channel `K_total` from `R14`
and the **actual** strict floor `m0^2` from `R15`, using `rho_star_sq` and `lambda_psi_strict` from `QW-2122`.

This audit does **not** promote any values into strict core; it only witnesses underdetermination compatible with
current exports.

## What N478 does not prove

`N478` does not prove:

1. any strict-derived canonical per-site vacuum vector `vpsi[0..11]`,
2. any strict-derived canonical per-site self-coupling families `g4[0..11]`, `g6[0..11]`,
3. discharge of `T168/T167/T166` or any strict decision of `F2(d)`,
4. strict-core selector closure or `QW-2191` discharge,
5. ToE closure.

