# MP7-044 — targeted negative controls

Scientific state: **PARTIAL_AUDIT COMPLETE FOR THE AVAILABLE MP7 CHECKPOINT**.

`audit/MP7-044_negative_controls.json` records the executed controls. They confirm the
intended rejection behavior for:

- swapped covariance eigenvalue ordering;
- discarding the transformed quadratic metric;
- a rank-deficient compression basis;
- equal-total-volume geometry with a compensated gap/overhang;
- stale input hashes;
- applying the positive-coefficient phase proof directly at negative alternating field;
- treating phases of zero amplitudes as distinct states;
- dropping the nonlinear polar gradient term away from stationarity;
- substituting M4 for M7 in full finite-N fluctuation formulas;
- the four semantic interpretation gates (local energy equality/globality, phase root/
  full stationarity, auxiliary Gibbs law/unique dynamics, and orbit multiplicity/selector).

The current handoff does not contain the full R7O3 leaf/witness tree, so the deleted-leaf,
changed-threshold, altered-c and corrupted-witness implementation controls cannot be
freshly replayed here. They remain inherited from the R7O3 intake rather than being
counted as new MP7 tests.
