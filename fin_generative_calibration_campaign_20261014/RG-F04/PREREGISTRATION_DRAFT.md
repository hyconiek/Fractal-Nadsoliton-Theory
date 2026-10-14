# Prospective preregistration draft — NOT an empirical registration

**Target:** PHYS-014-type N2 pair, q12, two distinctly addressable devices. Only a proposed engineered realization; no evidence of working device exists.

**Primary:** FIN A7 versus fixed trace-matched flat P7 for the **complete stationary 13-outcome difference law**, gain `g3`. `INVALID` retained. Independently calibrate physical realized gains and source compilation across `g in [2.97,3.03]` or **replace** that hypothetical box with measured bounds before freezing a real test. Require certified bound of observed TV error per model ≤0.005 including detector effects allocated to that category. All additional errors reduce rho as stated in RG-F03.

**Preparation:** physically verify equilibration after resetting labels; distinguish g0 uniform calibration from actual g3 steady law. Heat-bath worst-case decay exp(-2t) in model time only; timing in seconds needs independent g0 relaxation and rate-stability checks. Alternative kin­etics require their own burn-in proof.

**Calibration:** disjoint, time-stamped phase-resolved pair preparations; fibre transition rows, true readout confusion/invalids, disturbance and control/compilation/stability, detector drift, independent clock, replicate/reset protocol. Freeze channel IDs, software hashes, data schema, model variant and error budgets before first validation outcome.

**Validation:** iid independently reset full 13-alphabet observations *or* a separately proved dependence-aware procedure. Conditional demonstration: if the calibration justifies rho ≥0.00119114 and validation alpha 0.025, 40,278,676 iid draws suffice for epsilon=rho/3 using the conservative full-alphabet union bound. This is not an achieved or automatically practical target. Alternative power-efficient event test demands its own proven minimax event gap.

**Algorithm:** `decision_rule.py` returns FIN_ONLY, FLAT_ONLY, INCONCLUSIVE or NEITHER_MODEL_FITS from conservative model-set rejection; never force model selection. Never tune gains/readouts on validation. Reserve alpha=0.025 for calibration but define its actual procedure before deployment. A smaller validated rho requires a revised n *before* data opening.

**Secondary diagnostic:** complete joint 13×13 two-time law after independently characterized g0 quench and physical time conversion. Heat bath versus Metropolis is a model comparison, not an empirical kinetic identification. Do not use the already-open predecessor .4/1 results as a blind holdout.

**Stop:** if the phase-resolved interface, calibration, independence or actual uncertainty cannot be paid, stop with NO_NEW_PHYSICAL_INPUT. No acquisition or external experiment is authorized by this draft.
