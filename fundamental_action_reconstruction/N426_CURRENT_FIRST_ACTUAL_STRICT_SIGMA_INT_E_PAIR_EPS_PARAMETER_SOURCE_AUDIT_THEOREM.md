# N426 Current First Actual Strict Sigma-Int `E_pair` Eps-Parameter Source Audit Theorem

Status: `N426_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_E_PAIR_EPS_PARAMETER_SOURCE_AUDIT_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

The strict sigma-int lane now exports:

1. a strict sigma-int source upgrade (`F307/N418`),
2. an actual strict-core sigma-int → residual target-slot export-map object
   (`F311/N422`), and
3. strict-input candidate theta-pair instantiations (`F312/N423`, `F314/N425`),

but strict-core theta export remains absent and all candidate instantiations
fix an explicit parameter choice `eps = 1/2`.

This theorem packages the narrowest honest conclusion about that parameter:

```text
the eps amplitude parameter must not remain a silent free input on the strict
lane; it either remains an explicit parameter choice (pre-discharge), or is
supplied by a dedicated eps value object with strict provenance.
```

## Theorem-level conclusion

From `F315`, the strict sigma-int lane previously exported no dedicated eps
value object with strict provenance, and therefore all instantiations fixed
`eps` only as an explicit numeric choice.

On the current repo state (`F317/N428`), the strict sigma-int lane now exports
one dedicated eps value object with explicit strict provenance
(strict-source-upgraded by explicit premise):

```text
eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2.
```

Therefore the strict sigma-int → `E_pair` → theta pipeline remains:

```text
candidate-only with respect to generator/reduction status,
but no longer parameterized by a silent free eps choice on the strict lane.
```

## What N426 does not prove

`N426` does not prove:

1. impossibility of any future strict derivation of `eps`,
2. actual strict-core `theta_1`, `theta_2` export,
3. actual populated basis-pair instance,
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.
