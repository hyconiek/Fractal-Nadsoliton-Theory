# MP7-043 — proof dependency audit

Scientific state: **DONE_AUDIT; GRAPH ACYCLIC FOR CURRENT EXPORTED CLAIMS**.

The proposition-level graph is stored in `audit/MP7-043_dependency_graph.json`; the
machine check in `audit/MP7-043_dag_check.json` reports an acyclic graph.

The main circularity checks pass:

- phase alignment is derived from the model/Fourier positivity chain and does not assume
  a global minimum catalog;
- the `g=37/10` global-complement task is downstream of the minimizer reduction and is
  still OPEN, so no current theorem uses the expected three-root pattern as exhaustion;
- the finite-N results form a separate `CONDITIONAL_ON_ADDED_MODEL` branch and are not
  premises for Target P, phase alignment or the stationary bridge;
- the fixed-fixture 60-root census is not used to prove global-minimizer alignment;
- MP7-040 only classifies full stationarity at the same frozen amplitudes and is not used
  as a no-nearby-equilibrium theorem.

The remaining open nodes are left visible rather than bypassed: notably `G37_GLOBAL`
and `FOLD_REMAINDER`.
