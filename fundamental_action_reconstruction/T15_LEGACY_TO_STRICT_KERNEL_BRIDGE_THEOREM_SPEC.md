# T15 Legacy-to-Strict Kernel Bridge Theorem Specification

Status: `T15_CURRENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_THEOREM_SPEC_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N260`, the `T14` Source Topology lane is frozen as declared-scope
complete but closure-incomplete on the present export set.

Following `K1`, `K2`, `F2`, and `S2`, the highest-priority frontier is:

```text
legacy -> strict kernel bridge
or
explicit non-bridge strengthening
```

`T15` tracks only the positive bridge branch of that frontier.

It does not close the non-bridge branch, and it does not claim that a bridge
by itself would already yield strict-core selector closure or global
`QW-2191` discharge.

`T15` specifies the exact mathematical requirements for one theorem that would
relate `K_legacy_ont` and `K_strict_gate` at kernel-comparison level while
keeping role-transfer and closure overclaims explicit.

## The Two Kernels

The split is between:

**Legacy Ontological Kernel (K_legacy_ont):**
`K_legacy_ont(d) := alpha_geo * cos(pi/4 * d + pi/6) / (1 + 0.01 * d)`
Contains explicit ontological constants like `alpha_geo = 4 ln 2`.

**Strict Gate Working Kernel (K_strict_gate):**
`K_strict_gate(d) := cos(0.18575 * d + 0.16250) / (1 + 1.0 * d^1.8)`
A mathematically refrozen derivation structure with `eta = 1.8` non-linear damping.

## The T15 Bridge Requirements

A valid `T15` bridge theorem would still need rigorous control of three
distinct components:

### 1. Amplitude Absorption Map (The `alpha_geo` problem)

The strict gate kernel lacks the explicit info-geometry amplitude scalar
`alpha_geo`.

The bridge would need one operator `A_abs` clarifying whether:

```text
A_abs(alpha_geo) -> 1
```

or whether `alpha_geo` is preserved through some explicit normalization map.

This would still be only a kernel-comparison ingredient.
It would not by itself transfer legacy EM/EW role claims.

### 2. Damping Renormalization Map (The `beta_tors` vs `eta` problem)

Legacy uses linear torsion damping `(1 + beta_tors * d)`.
Strict uses non-linear micro-supported damping `(1 + beta * d^eta)`.

The bridge would need one renormalization flow `R_damp` relating the linear
legacy damping regime to the strict non-linear regime without silently
identifying them.

### 3. Phase/Frequency Conformal Shift

Legacy uses `(pi/4, pi/6)`.
Strict uses `(0.18575, 0.16250)`.

The bridge would need to explain whether this mismatch is:

1. one chart-level reparameterization,
2. one effective renormalized drift,
3. or one genuine non-bridge obstruction.

## Scope discipline

Even if a positive bridge theorem were later discharged, it would still need
to state its scope explicitly:

1. whether it proves only a structural kernel relation,
2. whether any legacy physical-role transfer remains separate,
3. whether it supports only local comparison or a stronger operator-level
   inheritance statement.

## Target Statement

If the positive bridge branch is ever discharged, the strongest honest
target would be:

```text
Actual Legacy-to-Strict Kernel Bridge Witness
  B_legacy_strict_actual_witness_v1 :
    K_legacy_ont -> K_strict_gate
```

but only in explicitly declared bridge scope.

By itself, such a witness would still not automatically establish:

1. legacy physical-role transfer,
2. strict-core selector closure,
3. global selector closure,
4. global `QW-2191` discharge.

## Hard limits

`T15` is only a specification.

It does not:

1. claim the bridge already exists,
2. silently merge `K_legacy_ont` with `K_strict_gate`,
3. close the non-bridge branch,
4. transfer legacy physical-role claims onto the strict kernel,
5. discharge `QW-2191` globally.
