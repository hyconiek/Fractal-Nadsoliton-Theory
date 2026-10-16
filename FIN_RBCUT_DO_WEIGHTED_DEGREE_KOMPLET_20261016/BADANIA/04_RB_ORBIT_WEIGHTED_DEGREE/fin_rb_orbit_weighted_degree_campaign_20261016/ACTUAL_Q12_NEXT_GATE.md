# One exact still-unpaid native q12 weighted-cut estimate

**Scientific status:** new FIN-specific equivalence has been proved, but the **asymptotic cut bound is OPEN**. This file is the next atom, not a claim that the counterfamily exists.

The actual and only relevant diagnostic class is the one from RB/FP: fixed exact `g_c=12/Lambda6`, actual MT connected domain and departure set, unchanged phase `A_N={R_N<=T(Z_N)}`, no changed energy/clock, legal count edges H between same-phase states or touching the original MT departure layer, stationary `W_N`, exact source `s_x=(N/2)W_N(x)L0_D m(x)` at all states, and nonnegative reversible conductances matching the *pointwise* native drift.

For a parity-preserving G_+-orbit representative `x_O` and any set J of genuine vertex orbits let

`b_O(J)=#{(a,b): a odd, b even, n_a>0, y=n-e_a+e_b in D_N, xy in H_N, [y] not in J}`.

This number lies in `{0,...,36}` and requires exactly the **actual** MT reflected-edge and frozen phase predicates. Set

`D_N(J) = sum_(O in J) |O| W_N(x_O) b_O(J)`;
`S_N(J) = sum_(O in J) |O| (N/2)W_N(x_O) L0_D m(x_O)`;
`Q_N = sup_J [S_N(J)]_+/(N D_N(J))`, with positive/zero infinity.

**Precisely missing theorem or obstruction:** prove `sup_(large N in12N) Q_N<infinity`; or provide a correct actual q12 sequence of orbit-union cuts J_N whose quotient goes to infinity, including exact population quantifiers and membership/reflection/weight proofs. Either outcome would determine the original *uniform vertex-rate feasibility* question: by the newly proved result `Q_N <= C_N^vert <= Q_N/a_0`, where `a_0=(103/15000)e^(-24)` is universal for this local domain. A single finite-N bad cut settles that N only, not all-N alternatives.

A negative ideal is the `Q_N=+infinity` case. A divergence with positive denominator would show **no** uniform cap for that sequence, even if every finite-N flow exists; this is stronger than failure of the former purely native-dominated sufficient *method* because the FIN-specific lower bound proves equivalence of the two caps in the real MT cylinder.

The remaining inequality is still globally nontrivial: weighted directed edge expansion for irregular phase-shell orbit cuts. It is not proved by positive level cuts, rare-shell stationary probability, the continuum phase profile, or toy rational flow checks. An exact comparison must retain the signed native source and the actual departure-layer edges.

**No automatic continuation:** One analytic monotonicity/weighted-isoperimetry argument or one real q12 countercut is the next research target. Do not repeat toy enumerations, alter the phase threshold, approximate `g_c` to classify ties, invent a source provider, or treat `Q_N` as a probability/statistical estimator.
