# AA: focused audit of the conditional critical-window action law

Date: 2026-11-07. Audited head: 789c4ebc5564cafe32033217d2e4571028f36053.
Disposition: ACCEPT_AT_STATED_DEPENDENCY_BOUNDARY.
This is a same-author focused audit, not independent review or formal verification.

## 1. Exact scope

The target is AW1 and its passage to AW10: residual radiation work is
O_S(epsilon^3[1+log(1/epsilon)]), and the actual canonical actions track
CD's three-channel ODE on 0<=t<=S epsilon^-2. Take fixed S and sufficiently
small positive epsilon; constants may depend on S and fixed signed model
coefficients, uniformly in w in [5/4,7/4].

The audit reconstructs the source ledger, adjoint estimate, projection
order, envelope budget and ODE comparison. It inspects how the cited
critical local-energy estimate is obtained. It does not independently
reprove the entire historical SP existence/coercivity theorem, VD virial
inequality, CM canonical chart, LE propagation theorem or CD spectral
coefficients. Those are explicit inherited inputs, not newly certified
facts. No blocking defect was found within this scope.

## 2. Dependency and circularity audit

The actual logical order is:

1. Global smallness from SP gives phase norm O(epsilon), without a damping
   law or a persistent shape-mode size O(epsilon^2).
2. QP supplies all three free-mass symbols and localization identities.
   CN retains actual continuum inputs and the complete forcing ledger.
3. ML proves the local/global estimates and absorbs localized feedback
   using the inherited VD/LF virial inequality. Its exponential time
   weight gives a constant depending on fixed S, not a uniform all-time
   estimate. CS separates the cubic source from higher terms.
4. CE supplies L1 preimages of the converted cubic source. TP2 uses
   uniform spatial adjoint kernels and a high-frequency trace bound to
   obtain integral E^2 <= C_S epsilon^2. It uses global smallness and
   source structure, not MA, RW or AW.
5. MA obtains the action remainder and envelope derivative bounds from
   that local-energy estimate and CM. RW transfers the causal quadratic
   response to the fixed CD polynomial. AW estimates its remaining work.

There is no circular appeal to the resulting action law in these inspected
steps. This conclusion concerns the displayed chain, conditional on its
historical input theorems. In particular the decay of the reference ODE
must not be substituted back as a premise for TP2.

## 3. Reconstructed adjoint estimate

For each work profile J=P_c B, Wcal-adjoint J is Schwartz. The explicit
Darboux eigenfunction and its first two spatial derivatives grow at most
polynomially in momentum, uniformly in position. Multiplication by the
Schwartz source transform therefore leaves an integrable amplitude.

For an open carrier, the roots of nu +/- sqrt(m_A^2+k^2) are simple.
The carrier set is finite and uniformly separated from threshold on the
compact mass band. A change of variable near each root gives the bound

    integral_(-d)^d min(t,2/|y|)dy <= C_d log(2+t).

For t>=2/d its value is exactly 4+4 log(dt/2); smaller t is bounded.
The complementary momentum region has a bounded denominator away from
zero, with integrable tails. Closed carriers have no roots. This proves
the spatially uniform primitive estimate, including derivatives through
order two, without estimating a stationary wave in global energy.

For a varying envelope A, the finite-time integration-by-parts identity is

    integral_0^L A(s+r) K_nu(r)dr
      =A(s+L) F_nu(L)-integral_0^L A'(s+r) F_nu(r)dr,
    F_nu(r)=integral_0^r K_nu(q)dq.

There is no lower endpoint contribution because F_nu(0)=0. The upper
endpoint is essential. MA gives sup|A|+integral|A'| <= C_S epsilon^2.
Thus all backward tests are O_S(epsilon^2 log(2+T)). Negative carriers
are included by conjugation. The work has no DC self-carrier: the zero
terms cancel exactly in PZ and RC. Difference carriers in PC and RZ
are retained even though their susceptibilities are closed.

## 4. Projection and physical momentum audit

In physical phase the projection is diag(P_c,P_c). With
J(r,p)=(r+i Omega^-1 p,r-i Omega^-1 p), its half-wave representation is
J diag(P_c,P_c) J^-1. It is not diag(P_c,P_c) in general.

For a continuous backward position test S, the momentum-source pairing is

    <S,P_c Omega B_0>=<S,Omega B_0>=<Omega S,B_0>.

The first equality uses P_c S=S. It does not commute P_c and Omega.
The correct bounded test is Omega S, estimated through
Omega^-1(-S''+m_A^2 S). It need not remain continuous.

The executable exact matrix countercontrol uses a projection not commuting
with a positive free mass matrix. The correct pairing is 2, whereas the
illegally reordered pairing is 5/2. It also checks that the conjugated
phase projection is idempotent and differs from componentwise projection.
This is an algebraic countercontrol, not a discretization of the RP PDE.

## 5. Odd channels and residual ledger

For output i, the input triples and coefficients are

| Output | Inputs | Coefficient a | Free masses |
| --- | --- | --- | --- |
| phi | phi,phi | -3/2 | 2;2,2 |
| phi | psi,psi | -3/2 | 2;m,m |
| psi | phi,psi | -3 | m;2,m |

The factor 1/4 comes from the two half-wave position reconstructions.
For every branch, (partial_t+sigma i Omega_i)B_b equals
 i Omega_i^-1(h_j h_k)+B_b(F_j,h_k)+B_b(h_j,F_k).
Multiplying by sigma a K and subtracting the correction gives exactly
AW's commutators, potential term and input-forcing terms. The potential
term is applied once to the sum of channel corrections. Both occurrences
of F in a same-field channel are required.

The response being transformed contains G_i only: N_X=0 in CN9.
The actual input F_C is still retained. Direct physical cubic forces and
F_C-induced cascades are different terms of that exact identity; retaining
both is not double counting. Initial continuum inputs vanish, so this
correction vanishes at zero. QP's whole-state initial shift is a different
map and must not be used here.

ML's weighted identity uses derivatives in an independent input momentum.
Its polynomial weight argument and localized-kernel estimates do not
single out the even output. QP's uniform derivative bounds cover all three
triples. CE's L1 kernel proof likewise uses bounded derivatives through
two and therefore extends to both odd channels. The output multiplier in
each physical conversion is the mass of that output, never a bound-mode
frequency.

| Exact residual class | Instantaneous bound | Integrated budget on T=S epsilon^-2 |
| --- | --- | --- |
| Bound chart reconstruction | L1 O(epsilon^3) | O_S(epsilon) |
| Bound/continuum quadratic terms | L1 O(epsilon E) | O_S(epsilon) |
| Full physical cubic/quintic forces | L1 O(epsilon^3) | O_S(epsilon) |
| Potential, commutators, F_V cascades | Positive-weight phase O(epsilon E) | Component L1 O_S(epsilon) |
| F_L cascades | Positive-weight phase O(epsilon^3) | Component L1 O_S(epsilon) |
| F_G and F_C cascades | L1 preimages O(epsilon^3) | Preimage O_S(epsilon) |
| Direct normal-form correction | Local pairing O(epsilon E) | Work O_S(epsilon^3) |

For the correction only pairings with exponentially localized work
profiles are needed. ML directly provides an H1_-4 half-wave bound;
there is no need to strengthen it to the H1_-3 norm defining E. This
clarifies AW's phrase 'local phase size': weighted duality suffices.
Every history row then costs epsilon^2 log(2+T) times O_S(epsilon).
Projection is removed against the already continuous adjoint tests only.

## 6. Final comparison and limits

RW's smaller O(epsilon^3 sqrt(log(1/epsilon))) error and MA's O(epsilon^3)
error fit inside AW's budget. The action polynomial is Lipschitz with
constant C epsilon^2 on the nonnegative O(epsilon^2) action box. Its
Gronwall factor is exp(C epsilon^2 T)=exp(CS). This establishes the
stated finite-S comparison. It does not establish uniformity as S tends
to infinity, pointwise monotonicity of actual actions or relative accuracy
for J when J is smaller than the absolute error.

In particular J(0)=O(epsilon^4) and the proved J bound
O(epsilon^3 log(1/epsilon)) do not resolve its leading evolution.
A vanishing leading shape chart shift at w=3/2 does not imply that the
physical shape mode stays zero.

The focused audit accepts AW at its explicit dependency boundary. Full
independent proof review remains OPEN. The next recommended mathematical
step is to audit the inherited VD/LF virial estimate and its absorption
into ML/CS, including boundary terms and all dependence on S. That is the
main upstream input not rederived here; another numerical PASS cannot
replace it. Empirical FIN and natural-source claims remain unchanged.
