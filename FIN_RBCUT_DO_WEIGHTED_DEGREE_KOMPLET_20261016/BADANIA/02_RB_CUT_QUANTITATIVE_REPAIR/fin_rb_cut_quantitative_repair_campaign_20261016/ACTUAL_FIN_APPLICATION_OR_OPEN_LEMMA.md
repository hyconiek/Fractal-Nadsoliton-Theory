# Frozen q12/A7 graph: exact consequences and remaining all-N obstruction

**Status: FINITE_THEOREM_APPLIES; UNIFORM_NATIVE_RELATIVE_CUT_NOT_PROVED; ACTUAL_H_IDEAL_SIGN_OPEN.**

For every genuine MT `D_N`, let `H_N` be the *unchanged* leaky support (within-phase edges plus cross-phase edges touching the actual departure layer), and use the exact reflected Gibbs/LOO source. Since all legal native opposite-parity count conductances are strictly positive, the new cut theorem applies verbatim.

Define the *actual*, fully typed finite-N statistic

`K_N^*= sup_{S subset D_N} [sum_{x in S}(N/2)W_x L^0_D m(x)]_+ / sum_{e=(x,y) in H_N; x in S,y outside S} W_x q^0(x,y)`,

where a positive numerator divided by zero is `+infinity`, `0/0` contributes zero, and the supremum also includes 0. No construction of a real-number value for this enormous statistic is claimed; this definition uses exact `g_c`, full multinomial weights, all boundary states, and the original **fixed** phase threshold.

**Consequence 1.** `K_N^*<infinity` if and only if every H-lower ideal has nonnegative source (old RB-CUT question). This gives an exact flow construction if the sign holds for that N.

**Consequence 2.** If `sup_(sufficiently large N in 12*N) K_N^* < infinity`, the constructive flow produces actual nonnegative reversible count conductances with exact native statewise parity drift and one uniform `CN` cap. Hence this stronger theorem would fully unlock RB-003 (and its previously specified phase-leak IBA argument), subject to that prior premise's own independent review.

**Consequence 3.** A finite-N cut S with positive demand and no allowed outgoing H edge is an exact ideal obstruction at that N. A sequence with the required unbounded-N quantifiers would refute the frozen graph route. If instead `K_N^*` diverges while each is finite, the *native-edge-dominated* constructive subclass fails to enjoy one fixed bound; **the original vertex-capped q12 family is not thereby excluded**.

**Consequence 4.** Every whole parity level cutoff S has the exactly computable necessary factor `1+C_F^+(S)/C_H^+(S)` as long as its denominator is positive. The existing RB theorem establishes that the source on every such cutoff is positive, but not that these ratios or irregular-subset ratios are uniformly bounded.

**Consequence 5.** The `RB-CUT-SHELL` small stationary mass estimate for F gives no lower bound on the H outgoing cut capacity for an *arbitrary irregular subset* and cannot alone bound K* or force any exact sign. It is consistent with an unfavorable cut concentrating on rare shell vertices.

## Exact remaining problem

For the actual q12/A7 MT domains at `lambda=0` and exact `g_c`, prove `sup_N K_N^*<infinity`, or show that this sufficient subclass does not cover some N via an exact cut, and then try the broader vertex-cap transshipment directly. Alternatively find a genuine H-lower ideal with negative full native source on the required asymptotic N sequence.

This handoff **does not** certify any of these alternatives. No full q12 count lattice, current remote Git checkout, original RB review or external apparatus was available. No population-size-specific computational result has been misreported as an all-N theorem.
