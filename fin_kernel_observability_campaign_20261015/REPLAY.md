# KO standalone replay

Requirements: Python 3, mpmath and numpy. All inputs are deterministic. No network, fitting, random draws or long finite-N trajectories.

1. From the extracted package root run `python independent_verify.py`.
2. It writes `INDEPENDENT_VERIFICATION.json` and prints a summary of eigenvalues, event probabilities and check count.
3. To validate final package file hashes run `python verify_package.py` (the builder/verify files are excluded from the output manifest only where explicitly stated in the manifest format).
4. Inspect KO-001/STRUCTURAL_HYPOTHESIS_AUDIT.md for comparator transfer premises, KO-002/EXACT_FAST_READOUT_KERNEL.md for conditional combinatorics/CLT, KO-003/BOUNDED_CONTRAST_THEOREM.md for the continuity event, and KO-004/DEFERRED.md for the nonactivated process proof.

Assertion types: rational projector and finite hypergeometric arithmetic (exact); interval transcendental kernel/probability bounds (mpmath.iv at 65 dps); small binary64 internal consistency (numerical tolerances). This does *not* machine-prove the asymptotic local Gibbs, martingale convergence, MT polynomial mixing or the previous LI risk theorem.

`SOURCE_INPUTS.sha256` binds supplied predecessor archives and the attached plan as inputs **outside** this handoff. The output ZIP contains only new KO work.
