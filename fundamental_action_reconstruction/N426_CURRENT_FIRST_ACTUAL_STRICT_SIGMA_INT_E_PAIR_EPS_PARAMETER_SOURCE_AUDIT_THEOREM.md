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
the current declared strict core exports no internal strict provenance source
for the amplitude parameter eps used in the sigma-int-driven E_pair generator.
```

## Theorem-level conclusion

From `F315`, on the current repo state:

1. the sigma-int driven `E_pair` generator uses `eps` only as a declared
   candidate parameter `eps ∈ [0,1]` (`T117/F270/N382`),
2. no dedicated strict-core eps value object with strict provenance is
   exported,
3. the strict-input instantiation artifacts fix `eps` only as an explicit
   parameter choice (e.g. `eps = 1/2`) and remain candidate-only outputs.

Therefore the strict sigma-int → `E_pair` → theta pipeline remains:

```text
parameterized by eps (candidate-only),
and cannot be cited as strict-core theta export or strict-core selector closure.
```

## What N426 does not prove

`N426` does not prove:

1. impossibility of any future strict derivation of `eps`,
2. actual strict-core `theta_1`, `theta_2` export,
3. actual populated basis-pair instance,
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

