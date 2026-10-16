# Combined RB / RBCUT proof audit and corrections

Date: 2026-10-16. Analytical audit under the inherited exact A7/MT premises;
not mechanized proof, independent peer review or physical evidence.
Historical source bytes are preserved. `degree`, `orbit`, `repair` and
`shell` below refer to the four campaigns inside the supplied RBCUT bundle.

## 1. What the five submissions actually pay

RB pays the exact native-source identity for lower ideals and shows, on a
synthetic twelve-vertex graph, why phase totals, parity prefixes and local
directional availability do not suffice. Twelve graph vertices do not make
that example a q12 FIN count process. No actual asymptotic bad cut is supplied.

Shell pays the exact integer Fourier step and the moving-shell inclusion
for both endpoints of every forbidden opposite-parity edge. MT convergence
and the null limiting boundary imply only `pi_D(Shell_N)->0`, without a
local anti-concentration rate. An `o(1)` normalized negative cut deficit is
not an exact nonnegative cut inequality. The Fourier interface shell is
distinct from the exponentially rare **energy-domain departure set E_N**.

Repair pays the finite maxflow identity

    K_N^*=max_S [s(S)]_+/U_H^+(S),

with positive/zero interpreted as infinity and zero/zero as zero. Its
collective signed repair can succeed where a nonnegative additive repair
cannot. An upward detour cannot replace one missing nearest-level edge
between the same two endpoints. None of this proves an actual all-N flow.
For real input capacities/demands the theorem gives a real flow. Rational
flows follow for rational inputs, not automatically for transcendental
Gibbs weights or unresolved numerical phase comparisons.

Orbit pays exact averaging and lifting for the parity-preserving subgroup
`G_+={j->j+2k, j->2k-j}`, of order twelve. It fixes m and preserves arrows.
For vertex orbits O,P and aggregate capacities u_OP, a quotient flow lifts as

    c_xy=(f_OP/u_OP)u_xy.

The row sum from each x in O is u_OP/|O| by transitivity. This proves
statewise divergence and vertex-cap preservation even when stabilizers
give different orbit sizes. Orbit-union cuts suffice for K*, by the finite
quotient maxflow theorem. At most 36 upward replacements per representative
does **not** bound the number of orbits; the symmetry reduction alone is
at most a factor of twelve in vertex count.

Degree supplies a genuinely model-specific advance. The actual MT cylinder
gives `103/1250<p_j<158/1875`; full twelve-label LOO proposals, followed by
departure rejection, give on each retained native edge

    a0 N <= q0(x,y) <= p_max N,
    a0=(103/15000) exp(-24)>2.59e-13, p_max=158/1875.

Indeed score differences are at most g_c Delta_A, while
`Delta_A<=2||A7||=2 Lambda6` and `g_c Lambda6=12`. The independent checker
pays the displayed coarse decimal lower bound using a positive rational
exponential series and a geometric remainder, not floating exponentiation.
This is ellipticity of the **given native reference**, not a lower-rate
premise imposed on candidate diagnostic kinetics or an activity source.

For orbit-union cuts, let `B_H(S)=sum_(x in S)W_x #H_up(x,S^c)` and
`Q_N=max_S [s(S)]_+/(N B_H(S))`. The finite comparison is

    Q_N <= C_N^vert <= K_N^* <= Q_N/a0,
    a0 K_N^* <= C_N^vert.

Here C_N^vert is the infimum of the **total incident conductance** cap,
not only outgoing upward flow. Every feasible c satisfies Q<=C by cut
balance. Native ellipticity and quotient maxflow give K<=Q/a0; native
total rate <=N gives C<=K. Conversely c_e<=CNW_x<=C u_e/a0 gives a0K<=C.
The comparisons include infeasibility/infinity. Thus uniform K is now an
equivalent target for the actual MT family, not merely a stronger sufficient
condition. This corrects application of the earlier generic-graph warning,
without invalidating its deliberately degenerating reference example.

## 2. A factor of two in the generic negative control

In the repair/degree four-vertex example, the retained chain is 0->2->3.
The two final conductances both equal `1+1/N`, with W=1. A reversible
generator has two departures at the middle vertex. Therefore

    max_x sum_(e incident x)c_e = 2(1+1/N),
    C_feasible = 2(1+1/N)/N.

The source degree JSON omits the factor two in C_feasible. The corresponding
repair assertion bounds a source/divergence, not the middle vertex's total
rate. For N>=2 the correct absolute total-rate bound is <=3, and the
normalized cap is <=3/2, rather than the stated <=2 and <=1 respectively.
The qualitative example remains valid: K*=N+1 diverges, while the actual
vertex cap is uniformly bounded. Its q0/N degenerates, so it is outside the
genuine MT ellipticity theorem. The separate 36-vertex example's C*=25/72,
K*=13/9 and its valid primal/dual certificate are unaffected.

## 3. Dual sign convention

The orbit dual formula `phi_x-phi_y<=beta_x+beta_y` on upward edges is
correct. At beta=0, phi is **nondecreasing** upward, not non-increasing as
one explanatory sentence says. A lower ideal L has the witness `phi=-1_L`
in this convention. In the earlier Farkas convention `alpha=-phi`, the
witness is `alpha=1_L`, with the earlier nonincreasing convention. A negative
ideal source violates the same primal gate in both notations. The displayed
dual formulas, unlike that sentence, do not need a sign reversal.

## 4. Same-parity edges need the phase permission test too

They do not affect Lm or Gamma(m,m), but they can change C,S and phase.
Consequently “retain the native same-parity edges” is licensed in the
rare-leak construction **only for edges in the permitted phase support**.
Setting all same-parity rates to zero is sufficient. Adding all native
same-parity edges may introduce ordinary interior phase changes, invalidating
the rare-leak implication even though every scalar drift equation survives.

This is not just a hypothetical support issue. Let n* be the balanced type,
T0=T(0)>0 and k=floor(sqrt(N T0/8)). For sufficiently large N in 12N put

    n^k=n*+k(e0-e2), m=0, C=2k, S=0, R=8k^2/N.

The even-to-even move 2->0 takes n^k to n^(k+1). The first is in phase A
(inclusive ties), the second in B. Both types and all their neighbors have
positive populations, F=O(1/N), and eventually lie in the central MT component
strictly below its fixed boundary energy. Hence neither endpoint belongs to
E_N. The crossing is genuinely possible in the actual count space.
No floating decision about equality to T0 is needed for this floor argument.

## 5. What remains unpaid

None of the five packages bounds Q_N uniformly on the actual all-N graph,
constructs the required exact-drift family or gives a genuine asymptotic
bad-cut sequence. RB-003 and FP-003 therefore remain NOT_ACTIVATED.
Universal q12 IBA without NC/IBA assumptions remains neither proved nor
refuted. The new review derives partial cut estimates and a deletion
obstruction separately; it does not promote either into a completed flow.
