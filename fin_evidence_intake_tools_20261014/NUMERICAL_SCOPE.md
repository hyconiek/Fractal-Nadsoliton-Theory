# Numerical scope

READY-001 is fail-closed tooling, not a new physics campaign.

## Decision arithmetic

`strict_decision.py` validates exact integer counts and uses `mpmath.iv` interval arithmetic for the 13-symbol confidence radius and the gain-grid model laws. Binary64 entries from the frozen F source model are converted through their exact integer ratios before interval evaluation. The grid-to-continuum payment uses the already-admitted conservative softmax Lipschitz factor. Correctness of the interval backend is still an assumption; this is not proof-assistant formalization.

Cases within the declared `numerical_guard` of the rejection threshold return `NUMERICAL_BOUND_OPEN` and are semantically inconclusive. Non-rejection is never reported as certified hypothesis membership.

The independent checker re-derived the sufficient iid threshold 40,278,676 for rho=0.00119114 and alpha_validation=0.025, and 373,539,331 for the F03 illustrative residual rho=0.00039114. These remain design counts, not observed records or feasibility claims.

## Process arithmetic

`strict_process_gate.py` binds the copied F model/process inputs by SHA-256 and checks the finite N2 software contract: circulancy, symmetry, constant diagonal, conservative generator, probability initial law, column-stochastic instrument, and the heat-bath reset shortcut only when consistent with the same bound model. Passing these checks does not certify that a physical device has those symmetries.

## Record parsing

`strict_records.py` uses decimal finite-value checks, UTC timestamp parsing, role-collision checks, trajectory monotonicity, and quarantine. Any quarantined row blocks inferential use of the package. Unique IDs do not imply iid sampling.
