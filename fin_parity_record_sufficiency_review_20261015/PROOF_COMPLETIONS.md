# LI independent proof audit and scope completion

Date: 2026-10-15. Status: analytical model audit under inherited exact spectral
and MT/UF/OBS premises. Not a formal proof-assistant or peer-reviewed result.

## 1. The inherited inputs are actually available here

The source correctly reports that several later review files were unavailable
to its executor. They are now present locally, and their actual manifests and
scientific input hashes are checked in verification.json. In particular the
MT nonuniform-mixing/lattice completion, OBS labelled path transfer and RI
feedback lemma supply the named premises. This binds the current proof;
it does not attest what the historical executor read.

Do not interpret SOURCE_INPUTS.sha256 as eight ordinary package paths: it
lists named external references. Seven match local producer files. The
original uploaded LI-plan bytes are unavailable; the actual repository plan
is separately hash-bound, without an invented formatting-equivalence claim.

## 2. Initial reference and future target

The constant-diagonal Gibbs identity gives the exact softmax conditional
tag law on backgrounds whose twelve completions remain in D_N. If one
completion is allowed and another forbidden, replacing the actual tagged
label by the forbidden one is a genuine count move. Its departure state
lies in MT's outgoing layer. The local partition bound and F>=b there give
the required small bad-completion mass.

Thus the source's initial joint-to-product parity-TV bound follows from MT
soft/fast first moments, retaining the boundary event. The background marginal
is unchanged by this product replacement. It is not an N-1 Gibbs start.
The full-label feedback comparison covers time sqrt(N), including the future
target, before any tag projection. After feedback deletion, background/parity
is exactly Markov and full q12 tag labels may be forgotten.

The full-to-background sign change occurs only when |Z_bg,N(1)|<=N^(-3/4).
Critical weak convergence and a continuous limiting density make this event
o(1): first compare to |Z|<=eta, then send eta down to zero. This is not a
small-N rate or permission to delete finite ZERO outcomes.

## 3. Pay occupation events without moving moments through TV

Let P_D be the stationary reflected labelled comparison, and let delta_N be
the sum of labelled initial error, its local exit error and feedback path TV.
For every background-path event E_N, even an N-dependent one,

    Q_bg,N(E_N) <= P_D(E_N)+delta_N,   delta_N->0.

Product replacement does not change the background marginal. Under P_D,
one-tag deletion changes m and w by O(1/N). Stationarity and MT polynomial
moments give, on T=sqrt(N)/2,

    E_D integral w_bg^2 dt = O(N^(-1/2)),
    E_D integral m_bg^6 dt = O(N^(-1)),
    E_D integral m_bg^2 dt = O(1).

For any fixed eta,K>0, apply Markov first and transfer only probabilities:

    Q_bg,N(integral w_bg^2>eta) <= C N^(-1/2)/eta+delta_N,
    Q_bg,N(integral m_bg^6>eta) <= C N^(-1)/eta+delta_N,
    Q_bg,N(integral m_bg^2>K) <= C/K+delta_N.

No factor sqrt(N) multiplies delta_N. The source's bound
`|q-m_bg|<=C(||w_bg||+|m_bg|^3)` then pays the time-integrated squared error.
UF/OBS soft C-tightness and fast-tube control separately pay sup|q|->0.
These arguments use the actual parent via path coupling; no reflected
controller is added to the physical model or the original ensemble.

## 4. Compensated likelihood errors need only probability bounds here

Under Q the background is independent of the fair tag. Conditional on the
entire background path, every coefficient a(t) below is deterministic and

    E_Q[(integral a dM_S)^2 | background]=2 integral a^2 dt.

Consequently, if integral a^2->0 in probability, conditional Chebyshev gives
an unconditional error bound by
`E[min(1,2 integral a^2/eta^2)]`, which tends to zero. This is enough; an
unpaid unconditional expectation bound under Q is unnecessary.

Apply this to a=q-m_bg. The analogous M_J isometry uses coefficient 1/2.
Since `integral q^4 <= sup(q^2) integral q^2=o_P(1)`, the quadratic martingale
remainder vanishes. Moreover `integral q^2 dJ` is tight after compensation,
so the cubic Taylor error also vanishes. For |q|<=1/2, the elementary series
bound is `|log(1-Sq)+Sq+q^2/2| <= (2/3)|q|^3`. Its complement has probability
o(1). This pays the likelihood expansion without observing attempt marks.

## 5. Make the stochastic-integral localization explicit

The source's partition/isometry route is valid with this precise convention.
Clip the *analysis integrand* Z_bg,N to [-R,R]. This does not change a process
law. On sup|Z_bg,N|<=R it is the original score. Its clipped process is bounded
and C-tight, so at a fixed partition its left-step approximation has mean
integrated squared error tending to zero as mesh->0 after N->infinity.
Use the conditional fair-tag isometry to transfer this error to the score.

For every fixed partition, joint background and tagged-martingale convergence
includes Z_bg,N(1), although observations end at1/2. Background and tag are
independent under Q; their actual jumps have zero cross variation. The
record and N^(-1/4)M_S differ uniformly by at most2N^(-1/4). At fixed clipped
partition this also pays the endpoint replacement in the score sums.

Then let R->infinity. Background weak convergence to a nonexplosive continuous
path makes the large-radius exception vanish. This order of localization
avoids silently taking expectations of unbounded scaled paths. It establishes
the joint terminal likelihood/target/path convergence used in LI-004, not
an all-time likelihood-process TV theorem.

## 6. Unit mass and UI are paid, not inferred from the amplitude CLT

For every N the conditional finite likelihood has mean one given background.
In the limit, conditional on continuous Z, the stochastic integral against
the independent variance2 Brownian record is Gaussian and

    E[L|Z]=exp[((1/2)^2*2/2-1/4) integral Z^2 ds]=1.

The integral is finite on each continuous compact path. Nonnegative likelihood
convergence with both expectations equal to one therefore supplies UI by
the source's min(L,K) truncation argument. No uniform exponential moment
transport through TV is needed. Weighted bad-background mass equals its
background probability by the same conditional mean-one identity.

## 7. The finite-cylinder argument really closes this bounded risk

The limiting likelihood has a Borel version as a function of (Z,Y0) paths:
use a convergent subsequence of stochastic-integral step sums. Finite rational
time coordinates generate their continuous-path sigma-fields. Truncate the
integrable target-weighted likelihood, condition on refining finite grids,
then approximate finite-dimensional measurable functions in L1 by bounded
continuous functions. This justifies the source's h_c approximation.

Joint convergence with likelihood and the discrete target, plus UI, transfers
its L1 error to finite N. Conditional expectation is L1-contractive. The
**paid product reference**, not merely joint weak convergence, permits hidden
background integration. Background tightness and a finite net on compact
observed-coordinate sets give the stated locally uniform convergence of H_c,N.

Hence the expectations of the maxima of the three target subdensities
converge. No posterior normalization is missing: these are densities relative
to the reference observation law. Future background coordinates used inside
h_c disappear under that integration; the decision sees observations only
through1/2. Initial parity information is in the full finite-N sigma-field
and is covered by the same projection bound, not silently excluded.

The argument is valid for the particular three-valued forecast and bounded
0-1 loss. It does not prove a sufficient statistic for all microscopic
quantities, full observation-law TV, all filters, entropy or mutual information.

## 8. Accepted terminal scope

Combining the paid initial/target comparisons, score, mass/UI and conditional
cylinders gives

    r_N(full one-tag parity through1/2 -> collective sign at1)
       -> r_*(continuous Y through1/2 -> sign Z1).

This closes the earlier RI-LIKELIHOOD-UI gate **only for this target, gain,
kinetic rule, preparation and finite scaled horizon**. MT's existential
burn-in and inherited exact spectral premises remain. No quantitative
finite-N rate, useful N0, equality of all information, physical observer,
natural source, units, locality or quantum/gravity conclusion is admitted.
