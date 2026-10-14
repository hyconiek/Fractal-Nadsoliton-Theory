# Predeclared total finite-law allowance (nonempirical)

Pre-open lock is `RG-D04/LOCK.json`, which binds four files: the frozen predictor, the 144/288-outcome predictions, a-priori per-law Taylor/filter error envelopes, and the preregistration. No withheld parent reference was accessed by `freeze.py`; only predictions, source-operator residuals and training documents were read.

The full one-time and joint-law model errors are bounded by the pre-open *analytic* formulas in `RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json`. Worst model term for FIN = 1.05184314055e-06 TV. The distinct exact-decimal-input + independent reference + numeric-source allowances are recorded in `RG-D02/ERROR_CERTIFICATE.json`. Conservative frozen aggregate margin = **model term + 1e−8 reference/source + 1e−7 additional numerical margin**, which is at most 1.16184314055e-06 for FIN. All are less than the fixed 0.01 primary target.

**Arithmetic status:** bound is conditional on interval math correctness and standard binary64 elementary operations. The 1e−7 margin is deliberately explicit rather than fabricated as a rigorous interval error proof for SciPy matrix `expm`. For nonnegative predicted finite outcomes, arithmetic/roundoff is also audited by independent serial rereads; no clipping/renormalization occurs.

**Predicted scope:** all 144 outcomes at t=.35 and .8 for each of four preparations and models, and 288 ordered binary/full joint outcomes at t1=.35, t2=.8, retaining the first-event complement. Model coefficients are rederived by a *frozen algorithm* per preparation, so this is source-based computational response prediction, not a preparation-independent autonomous closure.
