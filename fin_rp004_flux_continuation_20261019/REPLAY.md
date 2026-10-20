# FIN RP004-A/B replay

`python3 verify_rp004.py` — symbolic potential/remainder/parity controls (bounded).

`python3 check_spectral_energy.py` — one constrained diagnostic for the global reference energy coefficient, not ray localization or nonlinear PDE validation.

Verify the actual exterior sign with an independent local-energy integration; read `RP004-A/RETARDED_REFERENCE_BOUND_OR_OPEN_LEMMA.md` and `RP004-B/UNIFORM_REMAINDER_OR_FAILURE.md`.

Original strategy script `python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only` cannot succeed on the uploaded ZIP extraction, because source paths outside the archive are absent; this is a custody limitation, not an altered baseline.
