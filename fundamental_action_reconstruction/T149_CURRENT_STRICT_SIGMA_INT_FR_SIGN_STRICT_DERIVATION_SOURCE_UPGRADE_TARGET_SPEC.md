# T149 Current Strict Sigma-Int FR-Sign Strict-Derivation/Source-Upgrade Target Spec

Status: `T149_CURRENT_STRICT_SIGMA_INT_FR_SIGN_STRICT_DERIVATION_SOURCE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`T124/N389` already name one missing upstream strict-core prerequisite for any
honest sigma-int bridge/export-map claim:

```text
strict derivation / source-object upgrade for sigma_int_candidate
```

However, `T124` intentionally stays generic. The current repo also contains a
more specific concrete fact about the sigma-int candidate itself (`B4`):

```text
sigma_int_candidate := chi_FR(gamma_pi1) ∈ {+1,-1},
where gamma_pi1 is the nontrivial loop class in pi_1(C)=Z_2.
```

And `B4/B5` keep an additional constraint explicit:

```text
FR-sign support currently appears only under hybrid support (QW-1622),
so strict derivation/source-upgrade cannot be silently inferred.
```

Therefore the next honest “no-vague-object” move is narrower than re-asserting
`T124`:

```text
name the missing strict-core FR-sign derivation/source-upgrade ingredient
as one explicit future-only target object with explicit acceptance tests,
so future work cannot treat chi_FR(gamma_pi1) as strict-derived by wording.
```

`T149` does **not** claim a strict FR-sign derivation exists.

## Scope

`T149` is scoped only to the strict-core derivation/source-upgrade of the
existing sigma-int candidate definition:

```text
sigma_int_candidate := chi_FR(gamma_pi1) ∈ {+1,-1}.
```

It does not decide:

1. the residual-datum bridge/export-map object export (`N300/N301`),
2. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
3. strict-core theta export / pair population (`N18`),
4. ToE closure.

## Target object

If the answer is negative with respect to actual strict discharge but positive
with respect to sharp naming, export one future-only target object:

```text
Zeta_sigma_int_candidate_FR_sign_strict_derivation_source_upgrade_target_v1
```

with the intended meaning:

```text
export a strict-core, observer-free, noncyclic derivation/source-upgrade
showing that sigma_int_candidate := chi_FR(gamma_pi1) is strict-derived
on a declared strict domain (not hybrid-only),
and (optionally) export theorem-level gauge-quotient safety for that same
construction.
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Zeta_sigma_int_candidate_FR_sign_strict_derivation_source_upgrade_target_v1`
must at minimum provide:

1. **Declared strict domain:** an explicit declared strict-core domain
   containing a configuration-space object `C_v1` for the local protected
   sector (the repo may reference the already-used local support label
   `QW-2206_local_topological_protection_layer_in_local_B_tilde_1_sector`,
   but it must export a strict object standing for the relevant `C_v1` on the
   declared domain).
2. **Topological Z2 witness:** an exported strict witness that:

   ```text
   pi_1(C_v1) ≅ Z_2
   ```

   together with one explicit generator `gamma_pi1_v1` representing the
   nontrivial class.
3. **Strict FR-sign map:** an exported strict map object:

   ```text
   chi_FR_strict_v1 : pi_1(C_v1) -> {+1,-1}
   ```

   with an explicit statement of provenance:
   - either strict-derived on the declared domain, or
   - exported as a new explicit strict-side premise.

   In either case, no silent reuse of hybrid-only `QW-1622` is permitted.
4. **Strict sigma-int source upgrade:** an exported strict object/witness
   recording the status change:

   ```text
   sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1)
   ```

   and explicitly upgrading the sigma-int datum from candidate-only to
   strict-derived on the declared domain (no axiom-lane-only promotion counted
   as strict-core discharge).
5. **Gauge-quotient safety (optional but if claimed, must be proved):** if the
   discharge claims theorem-level gauge-quotient safety (to also discharge
   `T123/N388`), it must export:
   - a declared gauge action `G ⟲ C_v1`, and
   - an explicit invariance/quotient-level statement showing `sigma_int` is
     gauge-invariant under that action,
   while explicitly addressing the `A2/A3` gauge-mode degeneracy cautions (no
   silent gauge fixing counted as proof).
6. **Noncyclic contract:** no use of `theta_{1,2}` and no use of a populated
   basis-pair instance as input (respects `N18`).
7. **Observer-free contract:** no `K_obs`-indexed selection may serve as the
   uniqueness source for `sigma_int`.
8. **Selector discipline:** no implied admissible `S_sel_int`, selector
   closure, `QW-2191` discharge, or ToE closure.

## Relation to existing targets

If discharged honestly, this target is intended to directly satisfy:

1. `T124/N389` (sigma-int strict derivation/source upgrade), and
2. optionally `T123/N388` (gauge-quotient safety) **only if** item (5) above is
   explicitly discharged.

## Hard limits

`T149` must not claim:

1. that sigma-int strict derivation/source upgrade is already discharged,
2. that gauge-quotient safety is already discharged,
3. that any residual-datum bridge/export-map object is exported,
4. admissible `S_sel_int`, selector closure, `QW-2191` discharge,
5. ToE closure.

