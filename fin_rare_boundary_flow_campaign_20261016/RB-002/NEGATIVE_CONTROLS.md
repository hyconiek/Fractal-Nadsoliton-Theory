# RB-002 necessary negative controls

1. A correct native heat-bath Gibbs measure and drift do not alone provide a flow on a smaller edge support.
2. Global parity-level cuts are automatic because the given drift is divergence of positive conductances on the original *full* graph. They cannot detect deletions of interior cross-phase edges.
3. Phase-total balance is guaranteed by symmetry; it does not ensure finer within-phase ideal signs.
4. In an exact 12-vertex rational fixture, every ordinary within-phase level-prefix cut is nonnegative, every vertex with positive/negative source has an allowed up/down neighbour, and both phase sources balance. Yet an irregular H-lower ideal has negative source (-9/10). See verify_rb.py/EXACT_DIAGNOSTICS.json.
5. The toy is NOT an actual q12/A7/MT sequence and is not a universal-IBA counterexample.
6. The absence of a feasible capacity-free flow proof is NOT evidence that the actual frozen partition is impossible.
7. If a capacity-free flow were established, it would NOT automatically obey a uniform `CN` jump rate cap.
