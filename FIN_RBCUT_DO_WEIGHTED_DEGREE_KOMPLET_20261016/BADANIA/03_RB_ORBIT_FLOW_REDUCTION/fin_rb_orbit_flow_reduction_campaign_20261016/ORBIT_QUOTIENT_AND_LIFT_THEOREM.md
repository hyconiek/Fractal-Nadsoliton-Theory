# FIN RB-CUT continuation: exact parity-preserving orbit quotient and lift

Date: 2026-10-16. Scientific status: **PROVED finite graph theorem; conditionally applicable to every genuine finite MT q12/A7 domain; no all-N estimate obtained**.

## 1. Frozen objects and exact group

For each admissible `N in 12*N`, retain the exact original critical `g_c=12/Lambda6`, genuine MT connected count domain `D_N`, departure set `E_N`, `A_N={R_N<=T(Z_N)}`, the existing leaky upward graph `H_N`, full count weights `W_x`, exact reflected native heat-bath generator `L0`, and native conductances `c0_e=W_x q0(x,y)`.

Every `H_N` edge is oriented in increasing parity `m`, with increment `2/N`. Define `s_x=(N/2) W_x L0 m(x)` and native relative capacity `u_e=c0_e>0` on the retained edges. The candidate diagnostic reversible parity-only conductances are `c_e>=0` with `A_H c=s`; all same-parity rates may be set to zero (or separately retained at an extra native `N` attempt-rate cost).

Define the 12-element parity-preserving dihedral subgroup

`G_+={j -> j+2k, j -> 2k-j (mod 12): k=0,...,5}`.

This is the dihedral group of the six even indices. Every group element permutes count coordinates and fixes parity `m`. The A7 Fourier-projector/circulant interaction and native LOO generator are equivariant under this group; `F_g`, the uniform-centred connected `D_N`, `E_N`, `W`, and the radius-threshold phase partition are invariant. Thus every `H` edge is carried to another **upward** `H` edge, `u_{ge}=u_e`, `W_{gx}=W_x`, `s_{gx}=s_x`.

These premises are part of the frozen FIN algebra and actual MT-domain contract. The theorem below is elementary finite graph mathematics and applies if they are fulfilled; the current response does not claim new byte-level custody of the latest MT review or live HEAD.

## 2. Orbit network

Let `O=V/G_+` be the vertex orbits and for each orbit define

`W_O=sum_(x in O)W_x`, `s_O=sum_(x in O)s_x`.

Every original directed edge joins two **different** vertex orbits because it changes parity while parity is constant on each orbit. For orbit pairs `(O,P)` define aggregated reference capacity

`u_OP=sum_(e=(x,y) in H, x in O, y in P)u_e`,

and include an orbit-network arc iff `u_OP>0`. Multiple original arcs between an orbit pair are aggregated. Let `Abar` denote the incidence of this directed orbit graph.

**Theorem Q1: exact reduction for native-relative caps.** For every fixed finite `N` and `K>=0` the following are equivalent:

1. `exists 0<=c_e<=K u_e` on actual `H_N`, with `A_H c=s`.
2. `exists 0<=f_OP<=K u_OP` on the orbit graph, with `Abar f=sbar`.

Furthermore the sharp constant satisfies

`K_N^*=max(0, sup_(J subset O) [sum_(O in J)s_O]_+ / sum_(O in J,P outside J)u_OP)`,

with positive/zero interpreted as infinity and zero/zero as zero. In particular, **some union-of-orbits cut realizes any failing sharp cut inequality**, including the zero-outgoing-capacity/lower-ideal obstruction. No arbitrary nonsymmetric cuts need be searched to decide `K_N^*`.

**Proof: original to quotient.** Average a feasible `c` over the 12 action maps on original edges. Invariance of `s,u,H` ensures the resulting `c^G` is nonnegative, solves the same divergence and keeps the edgewise upper bounds. Set `f_OP=sum_{e: O->P} c^G_e`. Then `Abar f=sbar` by summing vertex divergence within each orbit, and `f_OP<=K u_OP`.

**Proof: quotient to original (constructive, exact).** Given a feasible `f` define on every original upward edge `e=(x,y)` with `x in O, y in P`

`c_e=(f_OP/u_OP) u_e`.

For a fixed `x in O`, transitivity plus reference-capacity invariance gives `sum_(y in P)u_xy=u_OP/|O|`. Similarly at `y in P`, `sum_(x in O)u_xy=u_OP/|P|`. The lift therefore has per-vertex divergence `(Abar f)_O/|O|=s_O/|O|=s_x`, is nonnegative, and satisfies `c_e<=Ku_e`. This is a fully explicit lift, not merely an averaged existence argument.

The sharp cut formula follows from the ordinary finite maxflow/mincut theorem *on the orbit graph*. Its equality with the full sharp factor follows from the two-sided feasibility equivalence. Even if a nonsymmetric violating full cut exists, a violating orbit-union cut must therefore exist. The same conclusion follows directly by group averaging a cut's indicator and coarea decomposition of its positive-variation cut capacity.

## 3. Exact reduction for the ORIGINAL uniform vertex-rate cap

**Theorem Q2: orbit quotient preserves unrestricted local vertex capacities.** For each `C>=0`, there exists a nonnegative reversible parity-changing count generator on the full actual `H_N` with exact native drift and

`sum_{e incident x} c_e <= C N W_x` at **every** vertex

if and only if there exists nonnegative orbit flow `f` with

`Abar f=sbar`, `sum_(P adjacent to O)(f_OP+f_PO)<= C N W_O` at **every** vertex orbit `O`.

**Proof.** Average any original feasible flow over `G_+`; the per-vertex degree bounds are preserved because `W` is invariant. Aggregate its edge flows. Conversely lift arbitrary quotient `f` with the explicit ratio formula above. By transitivity, each original vertex in orbit `O` has total incident flow exactly `sum_P(f_OP+f_PO)/|O|`, and `W_x=W_O/|O|`. Its exit rate is therefore at most `CN`. The cap equivalence holds without any positive lower bound on an edge relative to native `u_e`, so it reaches the **original RB/FP diagnostic class**, not only the stronger native-dominated subclass.

The minimal cap `C_N^vert` is the exact optimum of this finite quotient linear programme. The practical difficulty remains the number and geometry of genuine MT count orbits; no polynomial-in-N state count is claimed.

## 4. Complete finite dual, reduced to the orbit graph

Let `Bbar` be the unoriented vertex-edge incidence matrix on the orbit graph (each arc contributes `+1` to both endpoint orbit degrees). For a fixed C, feasibility is equivalent to the infinite but finite-dimensional exact family of inequalities

`sum_O phi_O s_O <= C N sum_O beta_O W_O`

for all real potentials `phi_O`, nonnegative `beta_O>=0` obeying

`phi_O-phi_P <= beta_O+beta_P` on **each upward** quotient edge `O->P`.

Proof of necessity: multiply `Abar f=sbar` by phi, use the displayed edge inequalities and the vertex-degree bounds. Sufficiency: Farkas' lemma / finite LP separating hyperplane theorem for `Abar f=sbar, Bbar f<=C N Wbar, f>=0`. The dual certificate is a genuine vertex-cap obstruction, not merely an ideal sign test or a bound on the stronger K*. All non-invariant full-graph certificates can be averaged to orbit-invariant ones by convexity and symmetry.

For `beta=0`, the inequalities demand a non-increasing potential along upward edges; its level sets recover capacity-free ideal obstructions. For general beta the dual detects bottlenecks that need more than the allowed `CN` rate even when all ideal signs are nonnegative.

## 5. Exact distinction: conditional finite theorem versus actual all-N q12 research

No actual huge `D_N` has been enumerated, no actual `K_N^*` or `C_N^vert` evaluated, and neither `sup_N K_N^* < infinity` nor `sup_N C_N^vert < infinity` has been proved. One actual orbit-ideal with negative source, correctly bound to exact `g_c`, would obstruct the leaky graph at that N; the required asymptotic quantifiers remain separate. An unbounded `K_N^*` may coexist with bounded `C_N^vert`.

The reduction exploits the symmetry already present in the original equations. It neither supplies a new kinetic source nor proves physical locality, QM, GR or ToE. It does not relabel a synthetic D6 fixture as the MT q12 count graph.
