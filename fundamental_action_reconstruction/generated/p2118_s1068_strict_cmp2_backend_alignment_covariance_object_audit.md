# P2118 S1068: strict CMP2 backend alignment-covariance object audit

- Status: `OPEN_PARTIAL_PROGRESS_WITH_TRACE`
- Result kind: `PASS_STRICT_CMP2_BACKEND_ALIGNMENT_COVARIANCE_OBJECT_AUDIT_WITH_TRACE`
- Rows exported: `4`
- Constructor stability pass rate: `0.0`

This stage replaces heuristic alignment sigma by a backend-exported alignment-covariance object (from quadrature error tensors) and reruns constructor audit with propagated backend alignment uncertainty.
No theorem-grade global closure claim is made.
