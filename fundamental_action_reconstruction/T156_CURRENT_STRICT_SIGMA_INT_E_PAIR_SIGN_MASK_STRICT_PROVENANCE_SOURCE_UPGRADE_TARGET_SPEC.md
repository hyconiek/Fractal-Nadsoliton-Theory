# T156 Current Strict Sigma-Int `E_pair` Sign-Mask Strict-Provenance Source-Upgrade Target Spec

Status: `T156_PACKET_READY_CURRENT_STRICT_SIGMA_INT_E_PAIR_SIGN_MASK_STRICT_PROVENANCE_SOURCE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After the strict sigma-int FR-sign prerequisite package is exported (`T149/F307/N418`),
the sigma-int-driven `E_pair` generator candidate (`T117/F270`) still contains one
explicit **free convention slot**:

```text
the fixed pair-indexed sign mask b_{i,k} ∈ {+1,-1}.
```

This spec names one narrow next missing ingredient:

```text
export an explicit strict-provenance (premise-based) sign-mask object
for the sigma-int-driven nad12 `E_pair` generator candidate,
derived from the exported strict FR-sign character on pi_1(C_v1).
```

No claim of strict-core `theta_{1,2}` export, target-slot population, selector closure,
or `QW-2191` discharge is permitted.

## Target object

Freeze one future-only target name:

```text
Xi_sigma_int_E_pair_sign_mask_strict_provenance_target_v1
```

intended meaning:

```text
an explicit sign-mask object b_{i,k} for k=0..11 and i∈{1,2},
whose provenance is strict-side and explicit (premise-based),
and which does not require theta inputs, populated instances, or K_obs.
```

## Acceptance tests

The target is discharged only if the repo exports an object satisfying all:

1. **Declared strict domain**

   The sign mask is derived from the exported strict-domain data:

   - `C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1`,
   - `pi_1(C_v1) ≅ Z_2` with generator label `gamma_pi1_v1`,
   - `chi_FR_strict_v1 : pi_1(C_v1) -> {+1,-1}` (explicit strict-side premise).

2. **Explicit definition**

   Let:

   ```text
   s := chi_FR_strict_v1(gamma_pi1_v1) = -1.
   ```

   Then for `k = 0..11` define:

   ```text
   b_{1,k} := s^k,
   b_{2,k} := s^(k+1).
   ```

   (So the mask is induced by the nontrivial Z2-character; no ad hoc pattern injection.)

3. **Typing**

   ```text
   b_{i,k} ∈ {+1,-1} for all i∈{1,2}, k∈{0..11}.
   ```

4. **Noncyclic / observer-free**

   The construction does not use:

   - `theta_1, theta_2` as inputs,
   - any populated basis-pair instance,
   - any `K_obs`-indexed selection as a primary source.

5. **No hybrid reuse**

   The only FR-sign source is `chi_FR_strict_v1` as exported in `N418`
   (explicit strict-side premise), not hybrid FR quantization artifacts.

## Hard limits

`T156` must not claim:

1. strict derivation (physics) of `chi_FR_strict_v1` (it is premise-based),
2. strict derivation/uniqueness of the full `E_pair` generator,
3. strict-core `theta_1`, `theta_2` export,
4. discharge of `N302/N395` (object-support above the export-map object),
5. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
6. ToE closure.

