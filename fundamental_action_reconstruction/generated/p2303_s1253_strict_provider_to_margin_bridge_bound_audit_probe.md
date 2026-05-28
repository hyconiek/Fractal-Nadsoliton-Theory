# P2303/S1253 — provider-to-margin bridge bound audit

- Status: `OPEN_BRIDGE_OBLIGATION_STRICT_DIRECT_BOUND_REFUTED_WITH_TRACE`
- Required P2302 lift: `0.0068`
- Strict direct coefficient floor: `0.005223026434758646`
- Direct floor gap: `-0.0015769735652413535`
- Signed coefficient total: `-0.13524712375150189`
- Verdict: strict direct bridge is not proven; current direct contracts refute the requested bound.
- Theorem fingerprint: `6ac8e1bf616df11691193a0819f79f5680b780e02e43db86b7a637fbf500eef1`

## Guardrail statement
Positive-sum and norm aggregations are numerically large enough, but P2303 does not promote them because no strict provider-to-margin aggregation theorem is exported. G1/G3 remain open; no QW-2191, selector, legacy-kernel, or ToE closure is claimed.

## Next honest step
Either derive a strict norm-to-margin aggregation theorem from the P2300 operator basis, or lower the P2302 policy-lock candidate to a bridge bound that the strict direct coefficient floor can certify; do not update G1/G3 before that.
