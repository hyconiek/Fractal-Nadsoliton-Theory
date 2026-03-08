# T16 Legacy-to-Strict Kernel Nonbridge Strengthening Theorem Spec

Status: `T16_CURRENT_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_THEOREM_SPEC_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N261`, the highest-priority frontier remains explicitly bifurcated:

```text
legacy -> strict bridge
or
explicit non-bridge strengthening
```

`N123` already discharges one current-repo-state package-level nonbridge
theorem.

What still does not exist is a symmetric future-route theorem spec for
strengthening that nonbridge result at the same kernel-comparison level where
`T15` now tracks the positive bridge branch.

`T16` does exactly that:

1. it does not claim a stronger nonbridge theorem is already discharged,
2. it does not claim no bridge can ever exist,
3. it only specifies what a stronger kernel-level nonbridge theorem would need
   to show if the positive bridge branch keeps failing.

## Starting point

Current strongest negative statement already exported:

```text
N123:
the full legacy kernel/package and the current strict side
are nonbridged at package level on the current repo state
```

That theorem is still narrower than a component-level kernel-comparison
nonbridge strengthening.

## The T16 strengthening requirements

A valid `T16` strengthening theorem would need explicit obstruction control on
the same three comparison layers isolated by `T15`:

### 1. Amplitude non-absorption obstruction

One theorem-level obstruction showing that no admissible current map

```text
A_abs : alpha_geo -> strict_amplitude_normalization
```

is exported at the declared comparison scope.

### 2. Damping non-renormalization obstruction

One theorem-level obstruction showing that no admissible current map

```text
R_damp : (1 + beta_tors * d) -> (1 + beta * d^eta)
```

is exported in a way that preserves the declared comparison contract.

### 3. Phase/frequency non-conformal obstruction

One theorem-level obstruction showing that the phase/frequency mismatch

```text
(pi/4, pi/6) vs (0.18575, 0.16250)
```

is not currently discharged as an admissible chart reparameterization or
effective shift relation.

## Strengthened target statement

If all three obstruction layers were later discharged, the strongest honest
future target would be:

```text
Actual Legacy-to-Strict Kernel Nonbridge Strengthening Witness
  NB_legacy_strict_strengthening_actual_witness_v1 :
    (K_legacy_ont, K_strict_gate)
      -> explicit_legacy_strict_kernel_nonbridge_strengthening_target_v1
```

meaning only:

1. the current positive bridge branch is obstructed at kernel-comparison
   level on the declared scope,
2. package-level nonbridge is strengthened by explicit component-level
   obstructions,
3. no silent bridge inheritance is admissible on that scope.

## Scope discipline

Even if a positive `T16`-style strengthening theorem were later discharged, it
would still need to keep all of the following explicit:

1. no claim that no bridge can ever exist in every future refinement,
2. no claim that the strict kernel/package is false,
3. no claim that selector closure follows automatically from nonbridge,
4. no claim that ToE is closed.

## Hard limits

`T16` is only a specification.

It does not:

1. discharge a stronger nonbridge theorem,
2. close the positive bridge branch,
3. prove permanent impossibility of any future bridge,
4. discharge selector closure,
5. discharge global `QW-2191`.
