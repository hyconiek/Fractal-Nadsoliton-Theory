# F324 First Actual Strict Sigma-Int `E_pair` Sign-Mask Strict-Provenance Source-Upgrade Packet

Status: `F324_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_E_PAIR_SIGN_MASK_STRICT_PROVENANCE_SOURCE_UPGRADE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After the strict sigma-int FR-sign package is exported (`F307/N418`), the
sigma-int-driven nad12 `E_pair` generator candidate (`T117/F270`) still carried
an explicit fixed sign-mask choice `b_{i,k}` that was not sourced on the strict lane.

The next honest move is narrower than theta export:

```text
export one explicit sign-mask value object with strict-side provenance
(premise-based via chi_FR_strict_v1 on pi_1(C_v1)),
so the mask is no longer a silent convention token on the strict lane.
```

`F324` executes exactly that move, staying below any claim of strict theta export
or selector closure.

## Inputs reused

1. `F306/N417`
   - strict domain object `C_v1` and witness `pi_1(C_v1)≅Z_2` with generator label `gamma_pi1_v1`.
2. `F307/N418`
   - strict FR-sign character `chi_FR_strict_v1 : pi_1(C_v1)->{+1,-1}`,
   - strict sigma-int value `sigma_int_strict_derived_v1 = chi_FR_strict_v1(gamma_pi1_v1) = -1`.
3. `T117/F270`
   - nad12 index set `k=0..11` and the need for a pair-indexed sign mask `b_{i,k}`.

## Packet result

`F324` exports one actual strict-provenance sign-mask value object:

```text
b_sigma_int_E_pair_sign_mask_strict_provenance_v1
```

persisted at:

```text
fundamental_action_reconstruction/generated/b_sigma_int_E_pair_sign_mask_strict_provenance_v1.json
```

Definition (premise-based strict provenance):

```text
s := chi_FR_strict_v1(gamma_pi1_v1) = -1
b_{1,k} := s^k
b_{2,k} := s^(k+1)
for k=0..11
```

## Meaning

This packet means only:

1. the sigma-int-driven `E_pair` generator candidate no longer needs a silent
   ad hoc sign-mask injection on the strict lane,
2. the sign mask is explicitly sourced (premise-based) from the exported strict FR-sign
   character on `pi_1(C_v1)`,
3. the result remains strictly below:
   - strict derivation/uniqueness of the full generator,
   - strict-core theta export,
   - object-support above the export-map object (`N302/N395`),
   - selector closure / `QW-2191` discharge,
   - ToE closure.

## Hard limits

`F324` does not claim:

1. strict derivation (physics) of `chi_FR_strict_v1` (it is premise-based),
2. strict derivation/uniqueness of `G_sigma_int_to_E_pair_generator`,
3. any discharge of `T2` / `N302` / `N395`,
4. actual `theta_1`, `theta_2` export,
5. admissible `S_sel_int` or strict-core selector closure,
6. `QW-2191` discharge,
7. ToE closure.

