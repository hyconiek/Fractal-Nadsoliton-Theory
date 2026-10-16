# RB-CUT continuation: sharp native-dominated transport and signed repair

Date: 2026-10-16. **NEW finite-network theorem**, applied only *conditionally* to the frozen q12/A7 MT graph. This is not a q12 all-N feasibility proof.

## 1. Exact objects (unchanged primary FIN)

Fix one `N in 12*N`. The state space `D_N`, true departure layer `E_N`, partition `A_N/B_N`, parity step `d=2/N`, and directed allowed opposite-parity graph `H_N` are exactly those of `fin_rare_boundary_flow_campaign_20261016/RB-001/EXACT_SOURCE_AND_CUTS.md`. The native opposite-parity graph `U_N` decomposes into `H_N` and deleted interior cross-phase edges `F_N`. Write

- `W_x>0` for the *full multiplicity-weighted* MT Gibbs count weights;
- `c^0_e=W_x q^0(x,y)=W_y q^0(y,x)>0` for native conductances;
- `s_x=(N/2)W_x L^0_D m(x)=(A_U c^0_U)_x`;
- `b_x=(A_F c^0_F)_x` for the deleted-edge divergence, so `s=A_H c^0_H+b`.

Every oriented edge `e=(x,y)` has `m(y)-m(x)=2/N`; incidence has `+1` at x (tail) and `-1` at y (head). All definitions include MT reflection and **every** count state in `D_N`. Same-parity rates may be separately retained from native HB because they do not contribute to the parity drift.

## 2. A sharp finite, simultaneously constructive cut theorem

For any finite oriented H with positive reference conductances `c^0_e`, a signed demand `s` with `sum_x s_x=0`, define for every subset `S`:

`S_s = sum_(x in S) s_x`, `C_H^+(S)=sum_(e=(x,y) in H, x in S, y outside S) c^0_e`.

Define the **sharp native-relative cap**

`K_star(H,s;c^0) = inf {K>=0: exists 0<=c_e<=K c^0_e, A_H c=s}`.

Then, exactly,

`K_star = max(0, sup_(S subset D, C_H^+(S)>0) S_s/C_H^+(S))`,

with `K_star=+infinity` if any S has `S_s>0` but `C_H^+(S)=0`. The max is attained if finite, and a capacity-respecting rational flow can be constructed by augmenting H with a super-source/super-sink and invoking ordinary max-flow. **No nonnegative correction-only premise is used.**

**Proof.** Necessity: `S_s=out_c(S)-in_c(S)<=out_c(S)<=K C_H^+(S)`. For sufficiency, add super-source arcs of capacities `s_x` to vertices with `s_x>0`, sink arcs of capacities `-s_x` from vertices with `s_x<0`, and arcs of capacities `Kc^0_e` for allowed H edges. A source/sink cut indexed by S has capacity equal to total positive supply plus `K C_H^+(S)-S_s`. Therefore all cuts have capacity at least the total supply iff all the displayed inequalities hold. Max-flow saturates every source and sink arc, and its restriction to H has incidence `s`. This proof is valid for **arbitrary disconnected finite H** and does not assume that the count parity alone is a Markov state.

`K_star<infinity` iff all H-lower ideals have nonnegative `sum_L s`; this is RB's older capacity-free ideal criterion recovered quantitatively: if `C_H^+(S)=0`, the complement `L=D\S` is a lower ideal and `s(S)=-s(L)`. Thus this new theorem sharpens the old *existence* test to the **least edgewise inflation factor relative to native HB**.

### Immediate FIN corollary (conditional on a uniform cut estimate)

For every N, if `K_star(N)<=K` for one N-independent constant K, the max-flow construction gives exact nonnegative reversible count conductances on the real `H_N` with every native parity drift and stationary `pi_D,N`. Their total count exit rates obey

`sum_y q_tilde(x,y) = sum_(e incident x) c_e/W_x <= K sum_(e incident x) c^0_e/W_x <= K N`.

The same-parity rates may instead be retained native, in which case the total bound becomes `(K+1)N`, with no effect on parity drift. **This would settle both RB capacity-free and uniform-rate gates at once**, but only after the specific q12 `sup_N K_star(N)<infinity` has been paid. No such all-N estimate is proved here.

### Exact native level-cut lower bound

For each parity level `ell`, put `S_ell={x:m(x)<=ell}`. Then the native upper-parity cut has no incoming U edges, and

`S_s(ell)= C_U^+(S_ell)= C_H^+(S_ell)+C_F^+(S_ell)`.

Where `C_H^+(S_ell)>0`,

`K_star(N) >= 1+C_F^+(S_ell)/C_H^+(S_ell)`.

A fully blocked cut with nonzero native traffic certifies `K_star=infinity`. This is exact at every finite N (and every MT boundary count); it does not establish an upper bound or control the irregular ideals. The shell result controls *total deleted mass* but cannot bound the denominator for every irregular S.

## 3. Exact separate correction-only criterion: useful but too restrictive

Suppose one demands **purely additive repair**, `c_e=c^0_e+delta_e` with `0<=delta_e<=D c^0_e`. This solves the exact drift iff `A_H delta=b`. By the same cut argument, such a repair exists iff, for every S,

`sum_(x in S)b_x <= D C_H^+(S)`.

For fixed N the minimal D is the analogous maximum cut ratio with demand b. If `D` is uniformly bounded, full rates are bounded by `(1+D)N`. But signed repairs (`-c^0_e<=delta_e` allowed) can work where **no nonnegative additive repair exists**. Consequently additive repair failure is **not** a q12 feasibility certificate. The 8-state exact fixture in `EXACT_FIXTURES.md` demonstrates this strict distinction.

## 4. Local monotone-detour obstruction

Every upward edge changes m by exactly `2/N`. Given any removed upward edge `x->y` in F, a directed H-path x->y, if it existed, would necessarily have exactly one edge (total m increment is one step). But the edge `{x,y}` is deleted and H has no parallel count transition of the same endpoints. Therefore **there is no positive upward H-path rerouting that single removed edge to its same endpoint**.

This does *not* obstruct collective compensation of several removed edges or signed redistributions of retained conductances. It excludes a tempting but invalid edge-by-edge construction strategy without making a claim about the existence of some other repair.

## 5. Strong separation from the actual uniform vertex-rate problem

A uniform bound of `K_star(N)` is sufficient, not necessary, for the KR/RB cap. Native-edge domination can be too restrictive if individual `c^0_H` become small. The four-state family with allowed native reference edge weights `1/N` and forbidden edge weights `1` has `K_star(N)=N+1` but a feasible flow on H of weight `1+1/N` with total rate <=2 when `W=1`. Thus one can have an artificial graded family with uniformly bounded *absolute* rates and diverging **relative** cap. The fixture is **not** the genuine q12 Gibbs family.

Failure of a uniform `K_star` upper bound alone must NEVER be promoted to a negative verdict for the original RB graph's unrestricted `CN` rate cap. It remains a restricted constructive criterion.

## 6. Exact scientific status

Proved: finite-network maxflow/cut formula, sharp reference-relative cap, additive-versus-signed repair separation, no single-edge upward detour, parity-level necessary ratios.

Unpaid: evaluate/upper-bound `K_star(N)` on all actual q12 `H_N` or find genuine negative H-ideal; alternatively establish actual vertex-capped flows not relatively dominated by native conductances. Nothing here derives a physical kinetic source, an actual all-N closed/leaky-phase family, a counterexample to IBA, or any ToE/QM/GR statement.
