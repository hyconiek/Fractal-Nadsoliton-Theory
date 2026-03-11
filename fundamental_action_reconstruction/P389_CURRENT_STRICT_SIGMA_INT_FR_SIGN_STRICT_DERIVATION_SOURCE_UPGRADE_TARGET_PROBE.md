# P389 Current Strict Sigma-Int FR-Sign Strict-Derivation/Source-Upgrade Target Probe

Status: `P389_EXECUTED_CURRENT_STRICT_SIGMA_INT_FR_SIGN_STRICT_DERIVATION_SOURCE_UPGRADE_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Probe whether the current repo already exports any **actual** strict-core
derivation/source-upgrade of:

```text
sigma_int_candidate := chi_FR(gamma_pi1) ∈ {+1,-1},
```

or whether the strongest honest result remains only:

```text
future-only target naming for the missing strict FR-sign derivation/source-upgrade ingredient
```

as specified by `T149`.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| sigma-int candidate definition present | YES | `B4` defines `sigma_int_candidate := chi_FR(gamma_pi1)` |
| local topological protection label present | YES (label) | `F36` uses `QW-2206_local_topological_protection_layer_in_local_B_tilde_1_sector` as precursor ingredient label |
| strict config-space object `C_v1` exported | YES | `F306/N417` export `C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1` |
| strict `pi_1(C_v1) ≅ Z_2` witness exported | YES | `F306/N417` export a strict `pi_1(C_v1) ≅ Z_2` witness + generator label |
| strict FR-sign map `chi_FR_strict_v1` exported | YES | `F307/N418` export `chi_FR_strict_v1` as an explicit strict-side premise (no hybrid FR reuse) |
| strict sigma-int status upgrade exported | YES | `F307/N418` export `sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1)` |
| future-only target naming still the strongest positive move | NO | the missing ingredient is no longer “only nameable”; an actual strict source-upgrade package is now exported (`N418`) |

## Exact verdict

The strongest honest current verdict is:

```text
actual strict FR-sign map export + sigma-int source upgrade: present (explicit strict-side premise; no theorem-level physical derivation claimed)
strict domain + pi_1(C_v1) ≅ Z_2 prerequisite: present
future-only strict FR-sign target naming: superseded by actual export (still admissible as a target name)
```
