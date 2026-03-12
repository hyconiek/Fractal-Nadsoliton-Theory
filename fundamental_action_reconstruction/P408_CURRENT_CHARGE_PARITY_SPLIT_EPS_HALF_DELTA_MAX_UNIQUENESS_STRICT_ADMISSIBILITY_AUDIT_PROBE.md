# P408 Current “Charge-Parity-Split eps=1/2 / delta_max Uniqueness” Strict-Admissibility Audit Probe

Status: `P408_EXECUTED_CURRENT_CHARGE_PARITY_SPLIT_EPS_HALF_DELTA_MAX_UNIQUENESS_STRICT_ADMISSIBILITY_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Audit one recurring proposed promotion that risks a strict false pass:

```text
(A) “derive eps = 1/2 from Charge-Parity-Split conservation”
(B) “declare delta_max is the unique pointwise O(2)-cut saturating the Void with information”
```

This probe answers only:

1. is any such derivation/uniqueness claim already exported as strict-core evidence, and
2. may the current repo remove the explicit non-claims (`QW-2191` discipline; no strict-core theta export)
   by citing such a principle today.

## Strict-admissible evidence reused

1. `T117`
   - sigma-int-driven `E_pair` generator contract includes `eps ∈ [0,1]` (free slot).
2. `F317/N428`
   - exported eps value object `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2`,
     explicitly classified as `strict_source_upgraded` (premise-based, not strict-derived).
3. `P407/N441`
   - theta-pair depends on admissible `eps` choices (eps is a real selector slot).
4. `T119`
   - positive-window corridor admits `delta_d ∈ (0, delta_max]` (free slot).
5. `F328/N440`
   - exported delta_d value object `delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max`,
     explicitly classified as `strict_source_upgraded` (premise-based, not strict-derived).
6. `P403/N437`
   - theta-pair depends on admissible `delta_d` choices (delta_d is a real selector slot).
7. `QW-2191`
   - kernel-alone mode-index/basis uniqueness is obstructed; a selector/symmetry-breaking ingredient is required.
8. `T159/P406`
   - strict-core upgrade target requires elimination/derivation of the exposed selector slots (`eps`, `delta_d`)
     before one may claim canonical strict-core O(2)-cut.

## Audit findings

### 1) eps = 1/2 is not strict-derived today (no “charge parity split” theorem exported)

The current strict sigma-int lane exports:

- `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2` (`F317/N428`),

but it is explicitly premise-based:

```text
classification = strict_source_upgraded
```

No strict-core theorem currently exported defines a strict “charge parity split” datum/law strong enough to
derive `eps = 1/2` as a unique consequence.

Moreover, `P407/N441` prove that the computed theta-pair varies across admissible `eps ∈ [0,1]` (at fixed
delta_d). Therefore `eps` is a real selector slot on the current strict sigma-int → theta lane.

So the statement:

```text
eps = 1/2 is derived (or forced) by strict charge-parity conservation
```

is **not** currently supported as a strict claim.

At most, it may be recorded as a *separated premise/axiom* in a non-strict scope, or as a strict-side
source-upgrade premise (already done), but it cannot be promoted to strict-derived without a new strict
derivation chain that closes the slot.

### 2) delta_max is not a unique strict-core “pointwise solution” today

The positive-window corridor (`T119`) admits a family:

```text
delta_d ∈ (0, delta_max]
```

The repo exports the specific choice `delta_d := delta_max` only as an explicit premise-based value object
(`F328/N440`), i.e.:

```text
classification = strict_source_upgraded
```

And `P403/N437` prove theta-pair dependence on admissible `delta_d` choices. Therefore `delta_d` is a real
selector slot and `delta_max` is **not** a unique strict consequence of the corridor itself.

So the statement:

```text
delta_max is the unique pointwise O(2)-cut solution
```

is **not** admissible as a strict-core theorem on the current repo state.

### 3) “Void information saturation” is not a defined strict objective today

The phrase “saturating the Void with information” is not currently exported as a strict, typed objective
functional whose maximizer/unique optimizer would be provably `delta_d = delta_max`.

Therefore it cannot be used as strict evidence unless:

1. a strict objective is defined and exported on a declared strict domain, and
2. a theorem is proved that it uniquely selects `delta_d` (and does not introduce new hidden selector slots).

Until then, it is rhetoric/interpretation, not a strict discharged ingredient.

## Exact verdict

On the current repo state:

```text
eps = 1/2: exported as premise-based strict-source-upgrade (F317/N428),
           NOT strict-derived; eps remains a selector slot (P407/N441).

delta_d = delta_max: exported as premise-based strict-source-upgrade (F328/N440),
                     NOT unique; delta_d remains a selector slot (P403/N437).

Therefore it is NOT admissible to “clear” strict non-claims (theta export / QW-2191 discipline)
by declaring parity-split derivation or delta_max uniqueness in strict core.
```

## Next honest step

If one wants these moves, the next honest step must be explicitly one of:

1. **Strict-core route:** export a genuinely new strict-core derivation/selection ingredient eliminating the
   exposed selector slots (`eps`, `delta_d`) and satisfying the strict-core upgrade target (`T159`), or
2. **Separated non-strict route:** record “charge-parity-split eps=1/2” and/or “corridor saturation delta_d=delta_max”
   as explicit premises in a separated scope (axiom-augmented / strict-extension), without promoting them into
   strict core and without claiming `QW-2191` discharge.

No intermediate “verbal promotion” is admissible.

