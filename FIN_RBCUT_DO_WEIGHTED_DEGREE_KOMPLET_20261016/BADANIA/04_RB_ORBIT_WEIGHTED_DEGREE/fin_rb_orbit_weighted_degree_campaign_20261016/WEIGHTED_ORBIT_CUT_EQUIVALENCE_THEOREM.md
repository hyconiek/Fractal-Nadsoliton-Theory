# Exact weighted-incidence reduction of the original RB vertex-rate problem

**New conditional-to-actual structural theorem** (2026-10-16). All finite assertions below are analytical and exact under the frozen FIN q12/A7, MT local-domain and native-HB contracts. The remaining cross-N weighted cut estimate is **not proved**.

## Frozen finite problem

For each `N in 12*N`, use the genuine connected MT domain `D_N`, the unchanged phase partition `A_N/B_N` from `R_N<=T(Z_N)`, and the RB leaky edge support `H_N`: all legal within-phase one-copy count transitions and phase-crossing transitions with an endpoint in the **actual** native departure layer `E_N`. Restrict here to opposite-parity edges, directed from lower to higher `m`. Same-parity rates can be set to zero, or retained at an additional native clock cost; they never enter `L m`.

The exact count weight is `W_x=N!/(prod n_j!) exp[(g_c/(2N))n^T A7 n]`. Put `s_x=(N/2) W_x (L0_D m)(x)` using the FULL reflected leave-one-out native generator, including boundary states. By reversibility, `s=A_U u` with positive native edge conductances `u_e=W_x q0(x,y)`. Let `G_+` be the admitted 12-element parity-preserving dihedral subgroup and let `O_N=D_N/G_+`. The entire graph, source and weights are invariant under this action.

For any union-of-orbits set `S_J=union_(O in J) O`, define:

`D_H^+(x,S_J^c) = #{ y in D_N\S_J : x->y is an upward H_N edge }`,

`B_N(J) = sum_(x in S_J) W_x D_H^+(x,S_J^c)`,

`S_N(J) = sum_(x in S_J) s_x`.

An upward count replacement necessarily removes one odd label and adds one even label, so each `D_H^+` is an integer between 0 and 36. Define the exact weighted-incidence score

`Q_N := max_(J subset O_N) [S_N(J)]_+ /(N B_N(J))`,

with `0/0:=0` and `(positive)/0 := +infinity`. The empty/all-orbit cuts have zero numerator and contribute zero. **No Gibbs partition function is needed**; an overall weight normalization cancels.

Let `C_N^vert` be the infimum of `max_x sum_(e incident x)c_e/(N W_x)` over all nonnegative conductances `c_e` on H solving the exact statewise drift `A_H c=s`. If infeasible, set `C_N^vert=+infinity`. Let `K_N^*` be the sharp native-relative cap from the predecessor, `inf max_e c_e/u_e`, again `+infinity` if infeasible. Write `a_*=(103/15000) exp(-g_c Delta_A)>0` as proved in `UNIFORM_NATIVE_EDGE_ELLIPTICITY.md`. The explicitly computable conservative constant `a_0=(103/15000)exp(-24)` may replace `a_*` throughout because `||A7||_op=Lambda6` and `g_c Lambda6=12`.

## Theorem: weighted-degree / original-cap comparability

For **every** admissible finite N, including the `+infinity` feasibility convention,

`Q_N <= C_N^vert <= K_N^* <= Q_N/a_*`.

The supplementary bound `a_* K_N^* <= C_N^vert <= K_N^*` also holds. Hence, for the **actual frozen MT family**,

`sup_(N in 12*N sufficiently large) Q_N < +infinity`

**if and only if** `sup_N C_N^vert < +infinity` **if and only if** `sup_N K_N^* < +infinity`.

These equivalences cover capacity-free obstructions (`+infinity`) as well as the sequence-uniform finite kinetic cap. In particular, the artificial-network caveat “unbounded K may coexist with bounded C” is inapplicable to this actual uniformly elliptic local heat-bath reference; it remains a valid counterexample in the more general graph category.

### Proof of lower inequality Q <= C

For any feasible c with vertex cap `C N W_x`, sum its incidence equation over an arbitrary `S_J`:

`S_N(J)=sum_(edges leaving upward S_J)c_e -sum_(edges entering upward S_J)c_e`
`<=sum_(x in S_J) sum_(y outside S_J, x->y in H)c_xy`
`<= C N sum_(x in S_J) W_x D_H^+(x,S_J^c)=C N B_N(J)`.

The second inequality only uses `c_xy <= sum_(e incident x)c_e <= C N W_x` at every tail. If `B_N(J)=0` and `S_N(J)>0`, no feasible c exists. Take the supremum over J and infimum over c. This proves `Q<=C` including infinity.

### Proof of upper inequality K <= Q/a_*

For every upward H edge with tail x, the native reference conductance satisfies `u_e >= a_* N W_x` by the genuine MT cylinder / LOO rule. Therefore the native reference capacity of any orbit-union cut is

`U_H^+(J) = sum_(x in S_J, y outside, x->y in H) u_xy >= a_* N B_N(J)`.

If `Q_N<infinity`, then by its definition `S_N(J) <= N Q_N B_N(J) <=(Q_N/a_*) U_H^+(J)` for **every orbit-union cut**. The admitted exact G_+-orbit maxflow theorem (the preceding handoff) says these inequalities are necessary and sufficient for a full exact-divergence flow with `0<=c_e<=K u_e`, `K=Q_N/a_*`. Thus `K_N^* <= Q_N/a_*`. This argument does not require positive rates on edges outside H, nor a physical mediator.

For a flow dominated by `K u`, its vertex rate is bounded by `K times` the original native total-attempt rate `<= K N`, yielding `C_N^vert <= K_N^*`.

Conversely, any vertex-capped flow has `c_e<=C N W_x <=(C/a_*)u_e` at the tail, so `K_N^*<=C_N^vert/a_*` and therefore `a_* K_N^*<=C_N^vert`.

This completes the exact finite theorem. No step uses an unproved *uniform all-N* cut inequality; such an estimate is precisely what remains open.

## Closed-form orbit representatives (no full edge enumeration)

For an orbit representative `x` with orbit size `|O|`, the exact source and capacity can be computed symbolically using no more than the 36 upward replacement pairs `(a odd,b even)`:

`W_O=|O| W_x`, `s_O=|O|(N/2)W_x L0_D m(x)`,

`B_N(J)=sum_(O in J)|O| W_x sum_(a odd,b even) 1_{n_a>0} 1_{y=n-e_a+e_b in D_N} 1_{edge xy in H_N} 1_{[y] not in J}`.

`L0_D m` still uses ALL permitted native outgoing replacement pairs in both parity directions and their full twelve-label softmax, not a cubic asymptotic replacement. Each orbit image has the same outgoing *edge count* to each target orbit, because `G_+` acts transitively and preserves H. The expressions avoid enumerating all orbit member edges; they do **not** bound the number of distinct orbits or solve their large-N cut problem.

## What is NOT yet established

No independent proof of `sup_N Q_N<infinity`, no actual `q12/A7` negative orbit-union cut, no all-N flow construction and no temporal IBA counterfamily. The new result **identifies an equivalent weighted combinatorial inequality** rather than a mere stronger sufficient test. Exact support, the true MT domain and departure set, all boundary sources, and fixed phase threshold remain essential. An arbitrary synthetic graph with degenerating q0/N need not satisfy this theorem with an N-independent `a_*`.
