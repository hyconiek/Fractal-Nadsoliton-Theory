# N418 Current First Actual Strict Sigma-Int FR-Sign Map + Source-Upgrade Theorem

Status: `N418_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_FR_SIGN_MAP_AND_SOURCE_UPGRADE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Package, theorem-level, the strongest honest current statement about the
strict sigma-int FR-sign prerequisite package required by `T149`:

```text
export C_v1, export pi_1(C_v1) ≅ Z_2, export chi_FR_strict_v1,
and export sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1),
without hybrid reuse and without implied downstream closure.
```

This theorem does **not** claim theorem-level physical derivation of the FR
sign; it records the explicit strict-side provenance used for the exported map.

## Theorem-level conclusion

From `F306/N417` and `F307`, the current repo exports:

1. a strict configuration-space object:

```text
C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1,
```

2. a strict topological witness:

```text
pi_1(C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1) ≅ Z_2,
```

with generator label:

```text
gamma_pi1_v1 := the nontrivial loop class,
```

3. a strict FR-sign map object:

```text
chi_FR_strict_v1 : pi_1(C_v1) -> {+1,-1},
```

exported with explicit provenance:

```text
explicit_strict_side_premise_nontrivial_character
```

(no hybrid-only FR quantization reuse), and

4. a strict sigma-int source-upgrade object:

```text
sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1) = -1.
```

Therefore the current repo now satisfies the strict prerequisite package
requested by the acceptance tests of `T149` (items 1–4), with explicit
premise-based provenance for the FR-sign map and without any silent hybrid
reuse.

## What N418 proves

`N418` proves only:

1. the repo exports `chi_FR_strict_v1` as an explicit strict-side premise map
   on a declared strict domain with `pi_1(C_v1) ≅ Z_2`,
2. the repo exports the induced strict sigma-int source-upgrade value
   `sigma_int_strict_derived_v1`,
3. these exports are noncyclic and observer-free (no `theta_{1,2}`, no
   populated basis-pair inputs, no `K_obs` indexing),
4. no hybrid FR source is silently reused as strict evidence.

## What N418 does not prove

`N418` does not prove:

1. theorem-level gauge-quotient safety (`T123/N388`),
2. selector-track identification beyond overlay-only (`T147/N414`),
3. export of any residual-datum bridge/export-map object (`N301`) or discharge
   of `N300`,
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

