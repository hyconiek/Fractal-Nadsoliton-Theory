# RG-D05 — reserved finite-window computational validation

**Execution:** TEST_OPENED_AFTER_LOCK_VERIFIED. **Scientific verdict:** 36/36 FINITE-MODEL_COMPUTATIONAL_TESTS_PASS_CONDITIONAL_ARITHMETIC.

Reserved t=.35/.8 and `hidden delta2` reference outputs were constructed **only in verify_frozen.py** after checking all four D04 SHA-256 hashes. Three models × four preparations × (two complete 144-outcome one-time laws and one complete 288-outcome joint law) = 36 independently compared finite laws. FIN largest numerical TV = **2.73382629617e-08**, against largest equivalent stationary-average Markov TV **0.354837695323**; overall largest compressed TV across all models 2.97958956787e-08. FIN hidden delta2: TV single at .8 = 2.7187025338e-08, and joint TV .35→.8 = 2.73382629617e-08.

All 36 frozen ROM outputs are nonnegative and normalized within floating tolerance; the 288-law first readout includes its complement. The numerical discrepancies are smaller than the **pre-open** conditional finite-spectral error envelopes plus declared rounding allowance; no post-test tuning. Both sparse-state and tensor refresh propagations agree to ≤1e−12 in the scored one-time comparisons. For all scored laws ROM numerically outperforms the autonomous averaged-Markov baseline. This is a synthetic finite-generator computation, not a physical prediction.

**Critical scope:** predictor's 27-mode hidden basis for FIN is derived afresh by the same fixed Krylov rule from the **declared initial preparation and fixed readout**. It is not a single preparation-independent 27-mode memory law. Source/gain/clock are supplied, not independently derived.
