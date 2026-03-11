# P398 Current Strict Sigma-Int `E_pair` Sign-Mask Strict-Provenance Source-Upgrade Target Probe

Status: `P398_EXECUTED_CURRENT_STRICT_SIGMA_INT_E_PAIR_SIGN_MASK_STRICT_PROVENANCE_SOURCE_UPGRADE_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Probe whether the current repo already exports an explicit strict-provenance
sign-mask object for the sigma-int-driven nad12 `E_pair` generator candidate
(`T117/F270`), as targeted by `T156`.

This probe must not introduce any false pass:
no `theta` export, no target-slot population, no selector closure.

## Evidence checklist (T156 acceptance tests)

| Test | Verdict | Evidence |
|---|---|---|
| strict domain `C_v1` + `pi_1(C_v1)≅Z_2` + generator `gamma_pi1_v1` exported | YES | `F306/N417` |
| strict FR-sign character `chi_FR_strict_v1 : pi_1(C_v1)->{±1}` exported with explicit provenance | YES | `F307/N418` (premise-based; no hybrid reuse) |
| explicit sign-mask definition derived from `s := chi_FR_strict_v1(gamma_pi1_v1)=-1` exported | YES | `F324/N435` export `b_sigma_int_E_pair_sign_mask_strict_provenance_v1` |
| noncyclic + observer-free contract satisfied | YES | the construction uses only `(C_v1, gamma_pi1_v1, chi_FR_strict_v1)`; no `theta_{1,2}`, no populated instances, no `K_obs` indexing |

## Verdict

The repo now exports one explicit strict-provenance sign-mask object for the
sigma-int-driven nad12 `E_pair` generator candidate:

```text
Xi_sigma_int_E_pair_sign_mask_strict_provenance_target_v1: DISCHARGED (F324/N435)
```

This does **not** upgrade the full generator into strict derivation and does **not**
export `theta_1,theta_2` nor discharge `N302/N395`.

## Next honest step

Use the now explicit sign-mask provenance only as a **narrow** ingredient
upgrade on the sigma-int → `E_pair` → theta pipeline, while keeping:

1. generator strict-derivation status explicit,
2. theta export absent (`C50/N1` discipline),
3. `QW-2191` nonclosure discipline.

