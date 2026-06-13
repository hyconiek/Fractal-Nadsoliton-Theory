# P2685/S1635 strict Lagrangian/EOM reverse-closure obstruction matrix

Status: `P2685_STRICT_LAGRANGIAN_EOM_REVERSE_CLOSURE_OBSTRUCTION_MATRIX_NO_FALSE_PASS`

## Content-first grep
- `selector_independent_lagrangian_eom`: `245` hits
- `full_tensor_nonproxy_obligations`: `1863` hits
- `reverse_closure_or_helmholtz`: `201` hits
- `open_obstruction_witnesses`: `992` hits
- `forbidden_imports`: `13146` hits

## Reverse-closure rows
- `R1_selector_independent_reduced_terms`: available=`True`, satisfied_now=`False`; blocker=P2329/P2362 allow selector-independent progress, but only as reduced/parallel EOM work, not full tensor reverse closure.
- `R2_best_working_ltotal_anchor`: available=`True`, satisfied_now=`False`; blocker=P2316 explicitly keeps the full Task-3 theorem unclaimed.
- `R3_anisotropic_background_transport`: available=`True`, satisfied_now=`False`; blocker=P1974 residual rank=3 with 3 nonzero numeric samples.
- `R4_unified_nonproxy_gate`: available=`True`, satisfied_now=`False`; blocker=P1795/P1809 leave metric_full_tensor_closure open and TG1_BW locked by unified nonproxy residual.

## Anisotropic residual computation
Residual vector: `['-sigma1**2 - sigma1*sigma2 - sigma2**2', '3*H*sigma1 + dsigma1 - sigma1**2 - sigma1*sigma2 - sigma2**2', '3*H*sigma2 + dsigma2 - sigma1**2 - sigma1*sigma2 - sigma2**2', '-3*H*sigma1 - 3*H*sigma2 - dsigma1 - dsigma2 - sigma1**2 - sigma1*sigma2 - sigma2**2']`.
Symbolic Jacobian rank: `3`; isotropic limit zero: `True`; isotropic Jacobian rank: `2`; nonzero numeric samples: `3`.

## Verdict
P2685 executes the strict Lagrangian/EOM reverse-closure computation requested by P2684.  The forward/reduced route is real and selector-independent, but the reverse implication to nonproxy full tensor-resolved EOM/L_total closure fails on finite obligations: reduced rows do not supply a shared componentwise covariant lift, P2316 leaves theorem-grade full closure open, P1974 gives a nonzero anisotropic Bianchi-I residual, and P1795/P1809 keep the unified nonproxy gate locked.

## Next honest step
Do not promote L_total or claim ToE closure.  The next honest proof-grade step is P2686: a shared-background nonproxy component residual table for EA, EH, and ELg, with the Bianchi-I anisotropic residual as the first required row.  If that table cannot be made zero, export a no-go theorem for reverse closure from the current reduced/FRW scaffold.
