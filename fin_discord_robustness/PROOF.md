# New goal: microscopic robustness of the discord obstruction

The previous separability/discord goal is complete. This study attacks its
explicit unresolved exchange shift a=-5/12. No previous result is counted
again as a new discovery. The supplied strict reference is C=I/12+W/20.
All operators act on two twelve-level factors unless a different scope is
explicitly stated. No PDF or physical apparatus is generated.

## 1. The old singular block does not mark a stationary-state exception

The complete reciprocal class with the same mixed-projector Hartree flow
is H_a=V+aS+bI. Dropping the irrelevant scalar b, its four bands are

    Omega: e_O=a+(n-1)/2;
    off-diagonal symmetric: e_+=a+1/2;
    traceless diagonal: e_D=a-1/2;
    antisymmetric: e_A=-a-1/2.

The canonical V has respective eigenvalues (n-1)/2,1/2,-1/2,-1/2.
Let p_+ and p_O be the Lagrange interpolation polynomials for the first
two distinct target projectors. Then

    V=-I/2+p_+(H_a)+(n/2)p_O(H_a).

Their denominators vanish only at a=-1/2 or a=-n/4. The coincidence
at a=0 is harmless because the colliding bands have the same V value.
Thus [H_a,R]=0 implies [V,R]=0 outside those two values. In particular
a=-(n-2)/(2n)=-5/12 is admitted by the exact polynomial certificate.
The earlier invertible-block proof failed there, not the conclusion.

The two actual transfer exceptions must be treated separately. They admit
positive stationary states with cross-band coherence that do NOT commute
with V, so one must not silently extend the polynomial license to them.
However the u,v cross block has singular values kappa=(n-2)/(2n),
|a+kappa| and 1/n. At both exceptional values, n=12 gives inverse norm
at most 12, while ||H_a|| is respectively 5 and 7/2, both no greater
than the canonical bound 11/2. The old quantitative block proof therefore
applies with the same conservative constant.

**Uniform mixed-flow conclusion.** Every stationary density for EVERY
V+aS+bI, whose first marginal is the stated real strict C, has distance
from the one-sided classical-quantum set at least

    d0 = d_min_lower*(c0-c6)_lower/7392 > 3.45e-8.

This closes the exceptional-value gap for the full two-parameter family.
The bound is not an entropic discord value and does not assert existence
of a stationary density with those marginals for every a.

## 2. The pure-flow family: exact state symmetrization identity

The larger reciprocal class that agrees only on pure-projector flow is

    H=V+cP_s+B_a,  B_a=P_a B_a P_a Hermitian.

For any matrix R, set Rbar=(R+SRS)/2. The exact identity is

    [V,Rbar] = P_s[H,R]P_s.                                (1)

Indeed Rbar is block diagonal under S. On its symmetric block H differs
from V only by cI; V is scalar on the antisymmetric block. This proves
(1) without stationarity or positivity assumptions.

If R is stationary, Rbar is canonical-stationary. Its two marginals are
the common average (R_A+R_B)/2. Therefore, when this average equals C,
the canonical lower bound applies to Rbar. If R itself is exchange
symmetric, the same lower bound applies directly to R across all 4357
pure-flow-invisible interaction parameters.

One cannot instead assume that symmetrizing preserves classical-quantum
form. A positive CQ state with equal marginals C is

    R=C tensor C+epsilon(P_u-P_v) tensor (|u><v|+|v><u|),

for sufficiently small epsilon>0. It commutes with C tensor I and is
block diagonal in an eigenbasis of C, but its swap average does not commute
with that marginal operator. It is an explicit counterexample to that
otherwise tempting shortcut; it is NOT asserted stationary.

## 3. A uniform three-way inequality

Let d_A(R)=inf_(Q in CQ_A) D_tr(R,Q), let
Asym(R)=D_tr(R,SRS), and epsilon=||[H,R]||_1. Assume only that
(R_A+R_B)/2=C and H is in the pure-flow family above.

Pinch Rbar in the distinct first-factor eigenspaces of C to obtain Rtilde.
The invertible canonical u,v block gives

    c0-c6 <= n epsilon+2n||V|| ||Rbar-Rtilde||_1.

Equation (1) paid the first term without any bound on B_a. The spectral
gap and the Frobenius/trace norm comparison give

    ||Rbar-Rtilde||_1 <= n/d_min ||[Rbar,C tensor I]||_1
                     <= 4n(1+||C||)/d_min * d_A(Rbar).

The last inequality follows by comparison with an arbitrary CQ state and
partial-trace contraction, as in the preceding report. Distance to any
fixed set is 1-Lipschitz, even when the set is not convex. Consequently

    d_A(Rbar) <= d_A(R)+D_tr(Rbar,R)
              =d_A(R)+Asym(R)/2.

Using n=12, ||V||=11/2 and ||C||<1/6 proves the exact, outward-certified
inequality

    d_A(R)+Asym(R)/2+(d_min_lower/616)||[H,R]||_1 >= d0.    (2)

It also holds with A and B interchanged. For stationary R, the residual
term vanishes. Thus a stationary CQ state in the entire pure-flow class
would have to carry Asym(R)>=2d0. Section 5 gives an actual asymmetric CQ
branch for the AVERAGE-marginal premise. Saturation of the bound, or a CQ
equilibrium with BOTH individual marginals equal to C, is not thereby proved.
This does not construct an internal orientation source.

For an overall Hamiltonian scale g!=0 the commutator residual in (2) must
be divided by |g|. There is no absolute clock or dimensional calibration.
The state-exchange asymmetry refers to the two factors, not to directed
Z12 generators, parity of a particle species, or QW-2191 selection.

## 4. Scope and next obligation

The result is two-body. Permutation averaging on an arbitrary many-body
system is not silently substituted for (1). Equal one-body marginals alone
do not imply exchange symmetry. Matching pure flow is not matching the
full source T or its finite-speed controller.

The remaining distinction is precise: all-mixed-state matching excludes a
CQ equilibrium with its first marginal equal to C. The weaker average-
marginal premise permits the asymmetric construction below, even in the
mixed-flow family. Whether the general pure-only class admits a CQ equilibrium
with BOTH individual marginals C, or an internal source can prepare the
required resources, remains separate. No laboratory, kernel
bridge, selector, dimensional action, SM/GR or ToE closure follows.

## 5. An exact asymmetric CQ equilibrium and what is actually observed

Let C=I/n+gamma W, gamma=1/20, n=12, and supply
C'=I/n+[2gamma(n-1)/(n-2)]W=I/12+11W/100. Its positivity is already
paid by the exact strict bounds. For Q_i=I-|i><i| set

    tau_i=n/(n-1) Q_i C' Q_i,
    R_CQ=(1/n) sum_i |i><i| tensor tau_i.

This is explicitly classical-quantum on the first side, with positive
trace-one conditional states. Its first marginal is I/n, its second is
I/n+2gamma W, and their average is C. It is NOT a counterexample to
the theorem with both individual marginals C.

Every branch has disjoint coordinate support, so D R_CQ=0. The shifted
Hamiltonian H_-=V-S/2=(|Omega><Omega|-2D)/2 annihilates that support.
Thus R_CQ and each conditional branch are stationary for H_-.
Prepare a uniformly chosen classical label i on the first factor and
postselect Q_i on ONE supplied C' copy on the second. The success
probability is exactly (n-1)/n=11/12. The program and heralding remain
supplied resources in this particular implementation. Heralding is not
necessary in a different implementation: the n Kraus operators
Q_i/sqrt(n-1) form a complete Lüders instrument. Recording outcome i on
the first factor prepares R_CQ deterministically from one supplied C'.
This instrument is explicitly specified, not derived as an intrinsic FIN law.

Its swap average Rbar has both marginals C, is stationary for canonical V,
and has nonzero one-sided discord by the earlier theorem. R_CQ is generally
not classical on the second side; zero discord on one specified side is
not a fully classical two-sided ontology.

For every effect O with [O,S]=0, Tr(O R_CQ)=Tr(O Rbar). Moreover

    Pi_swap U_a(t) = U_0(t) Pi_swap,

where U_a denotes the UNCONTROLLED density channel generated by V+aS.
Any common fixed instrument whose individual outcome maps are swap-covariant
commutes with Pi_swap. Induction therefore proves equality of all finite
adaptive outcome records for the descriptions (H_-,R_CQ) and (V,Rbar),
within that explicitly restricted operational class.

This is not equivalence under every physically conceivable symmetric probe.
A labelled local measurement distinguishes the different individual marginals.
Coherent control of the unknown Hamiltonian is also a stronger oracle than
its uncontrolled density channel: with a |+> control, the X expectation is
Re Tr(U(t)R). Here it is cos(t/2) for (V,Rbar) but 1 for (H_-,R_CQ).
At t=pi they are 0 and 1. Both controlled Hamiltonians can respect swap;
thus it would be false to claim that symmetry alone excludes this access.
Absolute energy queries likewise are not included in the passive-channel
equivalence. The theory must state which records and control resources its
observer actually possesses.

## 6. Universal stationary broadcasters with identical local channels

Write P_+=(I+S)/2-D for the off-diagonal symmetric band, and P_-=(I-S)/2.
For every input density rho define

    E_+(rho)=2/(n-1) P_+(rho tensor I)P_+,
    E_-(rho)=2/(n-1) P_-(rho tensor I)P_-.

These are CPTP: their Kraus operators are
sqrt(2/(n-1)) P_± (I tensor |j>), and completeness follows from
Tr_2 P_±=(n-1)I/2. Their outputs lie in the +1/2 and -1/2 energy bands
of canonical V, so both are stationary for every input. Both marginals
of BOTH channels are exactly

    B(rho)=[I+(n-2)rho]/[2(n-1)].

Consequently complete local input-output tomography, even with a reference,
cannot distinguish these joint preparation channels. They are noisy
broadcasters, not exact cloning or broadcasting of an unknown state.

Their balanced mixture is the swap average of the deterministic CQ
instrument from section 5, for arbitrary input rho, not just the strict
program. It is therefore separable for every input. More strongly, for

    E_r=r E_++(1-r)E_-,  0<=r<=1,

the output is separable iff r=1/2, for EVERY input density. Indeed all
outputs have zero coordinate-coincidence diagonal entries. The partial-
transpose principal block on |ii>,|jj> has diagonal zero and off-diagonal

    (2r-1)(rho_ii+rho_jj)/(2(n-1)).

At least one such pair has rho_ii+rho_jj>0. Thus every unbalanced mixture
has a negative principal minor and is NPT. For uniform input diagonal,
compression onto the whole coordinate-coincidence space gives negativity
at least |2r-1|/n. Swap expectation is 2r-1, so a joint, exchange-invariant
observable DOES distinguish r. This is not the passive indistinguishability
of section 5: only the LOCAL channels are identical in this family.

For the strict program C', B(C')=C. The symmetric-band channel E_+ is
stationary even for every pure-flow-equivalent interaction V+cP_s+B_a,
because its entire output is in the fixed symmetric energy band. Thus it
also supplies a non-vacuous stationary strict completion uniformly in the
arbitrary antisymmetric microscopic block. It stores the prescribed data;
it does not itself derive the strict propagator or the program source.

Resource caution: the balanced channel is a deterministic GLOBAL CPTP
channel from one input copy. It must not be called an LOCC broadcast from
one program initially held at a fixed party. Its marginal depolarizing
factor is (n-2)/(2(n-1))=5/11 at n=12. The corresponding normalized Choi
state has partial-transpose eigenvalue -3/88, so that marginal channel is
not entanglement breaking. Producing it at the other party requires quantum
routing/communication or an appropriate prior quantum resource. Randomly
swapping a record and a quantum output is not automatically a free classical
operation. The earlier two-program-copy LOCC construction remains a
different, valid resource model.
This limitation is for the universal channel on an unknown input, possibly
entangled with a reference. A fixed known separable target can still be
locally prepared given its classical description; that is a different input
resource and does not contradict the channel obstruction.

The general projector-channel and noisy-broadcasting mechanisms are related
to established quantum cloning theory; no universal cloning optimum or
global mathematical priority is claimed here. The new FIN consequence is
an exact stationary family with identical complete local channels but a
fully classified and different joint entanglement resource.
