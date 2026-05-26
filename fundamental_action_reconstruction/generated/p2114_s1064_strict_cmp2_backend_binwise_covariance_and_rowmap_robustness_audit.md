# P2114 S1064: strict CMP2 backend binwise covariance + row-map robustness audit

- Status: `OPEN_PARTIAL_PROGRESS_WITH_TRACE`
- Result kind: `PASS_STRICT_CMP2_BACKEND_BINWISE_COVARIANCE_AND_ROWMAP_ROBUSTNESS_AUDIT_WITH_TRACE`
- Backend binwise object available: `True`
- Exported robustness rows: `4`

This stage replaces spread-based transport factor with backend-exported binwise covariance traces (when available) and audits interval stability under small row-map perturbations.
No global C3/D3 transport theorem or ToE closure claim is made.
