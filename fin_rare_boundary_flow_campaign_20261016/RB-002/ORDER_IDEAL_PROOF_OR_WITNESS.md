# RB-002 — capacity-free order-ideal criterion, new native defect reduction, and first open lemma

**Scientific status: `OPEN_ACTUAL_Q12_ALL_IDEAL_SIGN`.** Neither a legal q12/A7 countercertificate nor an all-N nonnegative flow is constructed.

## 1. Exact finite equivalence

Use the exact RB-001 `s_n=(N/2)W_n d0(n)` and directed opposite-parity graph `H_N`. Seek `c_e>=0` with `A_H c=s`, where `+1` denotes the lower endpoint. Total source is zero, by native reversibility. A finite max-flow/min-cut argument (super-source at `s_x>0`, super-sink at `s_x<0`, unbounded capacities on H-edges) gives:

`exists c>=0 with A_H c=s` **if and only if** `sum_{x in L}s_x>=0` for **every** lower ideal `L` of the directed H graph.

Necessity: sum the incidence equations over a lower ideal; every contributing H crossing edge leaves the ideal and has nonnegative `c`. Sufficiency: in an augmented source/sink network, the only finite cuts have the graph-side source set upward closed. Their min-cut demand condition is exactly the complementary lower-ideal inequality; max-flow/min-cut provides a feasible nonnegative flow. The finite statement makes no rate-cap claim and requires no irreducibility.

An explicit *capacity-independent* Farkas witness of failure is `phi=1_L`, `beta=0`, whenever `L` is a lower ideal with `sum_L s<0`; on all oriented H-edges, `phi_lower-phi_upper>=0` and `sum phi_x s_x<0`. A genuine all-N result must identify such an ideal for the true H_N on a sufficiently strong quantifier, or prove every ideal nonnegative. No assumed finite-size LP output can replace that.

## 2. Native-defect form and a new sufficient route

From RB-001, for every H-lower ideal `L`, the sign reduces to the actual native crossing weights:

`cut_down(F_N,L) <= cut_up(H_N,L)+cut_up(F_N,L)`.

This is a **fully model-specific exact** finite-N test: source and edge weights are the known Gibbs/LOO heat-bath objects; the prospective diagnostic rates do not enter. It gives one concrete possible constructive proof route: establish this inequality for every native order ideal, possibly using a valid monotone injection/weighted pairing from deleted downward-crossing traffic into upward-crossing traffic. Such an injection has *not* been proved. A native globally stationary cut identity alone does not provide it.

If that sign theorem succeeds, a feasible nonnegative count flow exists at that N, but no bound `sum_e c_e<=C N W_n` follows from max-flow/min-cut. The uniform cap remains an **independent** quantitative theorem gate.

## 3. Exact rational negative control — intentionally NOT q12

`verify_rb.py` builds a twelve-vertex layered reversible toy graph with three chains `a,b,c` over levels `0..3`, phase A containing a/c chains and phase B the b chain. The departure layer is `{a3,b3}`. Every edge has a strictly positive rational native conductance. One positive native edge `b0->a1` crosses phases with neither endpoint in the departure layer, so is deleted from H; the two allowed cross-phase edges `a2->b3` and `b2->a3` touch that layer. An additional `c` branch carries sufficient A-level flow to keep **every** ordinary within-phase parity prefix nonnegative. Phase-total source balances both vanish.

Nevertheless `L={a0,a1}` is an H-lower ideal with `sum_L s=-9/10`, despite local direction availability and all global *and per-phase* parity prefix cut tests passing. Its forbidden native incoming edge has weight `1`, while its permitted native outgoing edge has weight `1/10`. The exact Farkas witness `phi=1_L, beta=0` proves the toy's capacity-free impossibility. This independently establishes that coarse phase/level and singleton tests are **insufficient in general**; it does not establish a failing ideal in the actual q12 MT graph.

The toy uses positive rational conductances, constant positive reference weights and a stipulated departure set, **not** actual q12/A7 Gibbs weights, not the MT energy boundary, and no large-N family. It must never be relabelled an admissible FIN countermodel.

## 4. Actual attempted native decision and stopping point

The frozen q12 predicate is transcendental and the actual D_N contains an enormous count graph at the relevant asymptotic N. The accessible original MT source and FP frozen partition suffice to define every vertex/source/edge symbolically. They do **not** yield a proof of the native ideal inequality for *all* irregular ideals. The specifically mentioned newest review seed and current complete repository checkout were not available as verifiable original bytes. No huge enumeration or threshold retuning was attempted.

**First open lemma (RB-CUT):** for all sufficiently large `N in 12*N` at exact `g_c`, prove or refute

`for every lower ideal L of the actual H_N:  sum_{n in L}(N/2)W_n L0_N m(n) >= 0`.

A negative ideal supported by exact Gibbs bounds would be a legitimate finite obstruction at those particular N. An all-N positive proof would open RB-003's uniform-cap obligation. We obtained neither. This is a precise open model-specific lemma, *not* another name for generic IBA.
