# QW-2058: NONABELIAN FLAVOR GENERATOR (FIRST-PRINCIPLES) - STRICT SPEC

## 1. Purpose
- Close the remaining first-principles blocker (CKM/PMNS) under one frozen kernel.
- Keep strict constraints:
- no fit,
- no sector retune,
- no post-hoc parameter edits.

## 2. Inputs (Frozen)
- Kernel (from QW-2049): `omega, phi, beta, eta`.
- Mass-chain derivational inputs (from QW-1961): `q_assignment`, `gamma`, `delta_info`.
- Fixed references:
- CKM absolute matrix,
- PMNS absolute matrix.

## 3. Core Hypothesis
- Current flavor operators fail because they are effectively abelian/near-abelian in phase structure.
- Required upgrade: derive a nonabelian generator in a 3x3 flavor space:
- `G_sector = sum_a c_a(kernel, q_sector) * lambda_a`
where `lambda_a` are fixed SU(3)-like basis matrices and `c_a` are deterministic functions of kernel invariants and sector charges.

## 4. Strict Derivation Rules
- Coefficients `c_a` must be explicit deterministic expressions from:
- kernel moments,
- signed kernel geometry,
- cyclic charge distances.
- Allowed constants:
- integers,
- simple rationals,
- fundamental constants (`pi`, `e`) if needed.
- Forbidden:
- numerical fitting to CKM/PMNS,
- optimization loops that minimize CKM/PMNS error.

## 5. Candidate Construction
- For each sector `S in {up, down, nu, lep}`:
- build `G_S` from deterministic invariant blocks:
- pairwise distance invariants `I_12, I_13, I_23`,
- signed phase invariants `P_12, P_13, P_23`,
- chirality/order invariant `C_S`.
- define:
- `U_S = exp(i * G_S)` (or equivalent unitary map with exact derivation).
- predict:
- `V_CKM = |U_up^dagger U_down|`,
- `V_PMNS = |U_nu^dagger U_lep|`.

## 6. Acceptance Gate (Hard)
- `CKM mean rel% <= 15.0`
- `PMNS mean rel% <= 15.0`
- plus inherited triad constraints from QW-2055:
- mass_mean_rel_pct <= 15.0,
- mass_max_rel_pct <= 35.0,
- tau_charm_ratio_rel_pct <= 20.0,
- GW sep >= 0.0020,
- GW adv >= 0.30,
- GW auc >= 0.75,
- GW control_gap <= 0.0025.

## 7. Required Artifacts
- `QW_2058_NONABELIAN_FLAVOR_GENERATOR_NO_FIT.py`
- `report_qw2058_nonabelian_flavor_generator_no_fit.json`
- `RAPORT_QW2058_NONABELIAN_FLAVOR_GENERATOR_NO_FIT.md`

## 8. Failure Interpretation
- If QW-2058 fails, then under current kernel class the missing part is deeper than operator choice and likely requires:
- revised foundational mapping from Nadsoliton characteristics to flavor algebra,
- or explicit new first-principles axiom for flavor emergence.
