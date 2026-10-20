# CF-003: a paid causal-defect witness and a necessary countercontrol

Date: 2026-10-20. Source is the **existing ES diagnostic lift**, not a new
physical mediator or an independently sourced FIN law. Null refreshes stay
hidden and are included by the exact filter, not granted as observed marks.

## 1. Derive the full source before computing the comparison

The hidden state is z=(w,k), with w in E^r and k=0,...,r-1. At rate lambda,
k advances within the pre-sampled block; at wrap a fresh independent block
is sampled from the already frozen ES law T_r. Readout is R(w,k)=w_k and
the supplied stationary law is mu(w,k)=T_r(w)/r.

For i=w_k, the full-source observable rates are

    a_ij(w,k)=lambda 1_{w_(k+1)=j}, k<r-1, j!=i,
    a_ij(w,r-1)=lambda p_j, j!=i.                    (1)

The source filter in CF-001 sums over hidden phase/block and retains hidden
transitions with unchanged output. Its rates alpha, not merely stationary
averages of (1), are the actual observed-history intensities.

## 2. Why a raw source-fibre test is not necessary

Adjacent labels in a block have iid p marginals for all epsilon. Averaging
the absolute full-source rate defect against q_ij=lambda p_j gives exactly

    D_src(T)/T=2lambda(r-1)/r sum_j p_j(1-p_j)^2>0.    (2)

For each i and candidate j!=i, the mean of
`|lambda 1_{next=j}-lambda p_j|` over the next mark is
`2lambda p_j(1-p_j)`; then average over i and the r-1 nonwrap phases.

At epsilon=0, however, all block marks are iid p. The observed process is
**exactly** the native reset CTMC, so D_vis=C_vis=0, despite (2). Its full
source contains extra knowledge of pre-sampled future marks. This is an
exact preparation/filtration countercontrol against treating full-source
strong lumpability or D_src=0 as a necessary condition for this one law.
For p=(1/2,1/3,1/6), lambda=2, r=3, (2) equals `28/27` per time unit.

## 3. A positive visible defect for the actual ES memory family

For r=3, v=E_p h^2 and nonzero epsilon, the ES three-time formula implies
the following exact history-sensitive drift witness at t>0:

    W(t)=E[h(D_0)h(D_t)
              {L(h o R)(Z_t)-(Qh)(D_t)}]
        =(epsilon/3)v^3 lambda^2 t exp(-lambda t).    (3)

To see this, differentiate the ES moment difference at its **future** gap
s=0. Pair data and repeated-index terms agree, and Qh=-lambda h.
Finite source rates justify differentiation. Equivalently, in the discrete
uniformization series the only nonzero coefficient is the one-refresh
coefficient `lambda epsilon v^3/3`; all others vanish by the ES proper-
block marginal identity. No full-history filtering approximation is used.

The left side equals the same expression with the source drift replaced
by its projection to the actual observed filtration. Since |h|<=1,
`|h(D_0)h(D_t)|<=1` and |h_j-h_i|<=2. Therefore

    |W(t)|<=2 E sum_(j!=D_t)|alpha_ij(t)-q_ij|.

Integrating yields an actual **visible** lower bound, not just a large
sufficient source bound:

    D_vis(T)>=|epsilon|v^3/6
                        [1-(1+lambda T)exp(-lambda T)]>0. (4)

Thus the pair-conditioned covariance agreement cannot pay causal closure.
For the rational ES control, epsilon=1/2, v=1/4, lambda=2, T=1, (4) is
`[1-3/e^2]/768 >1/3072`, using e>2. For higher r the already paid bounded
r-time contrast also gives `D_vis(T)>=|Delta E|/2` by the path coupling
bound. The claim concerns each fixed chosen r, not simultaneous equality
of all finite-order history laws by one different process.

The entropy criterion provides another obstruction. With the same initial
law, CF-002 gives `C_vis=D_KL(P||Q_path)`. A bounded snapshot contrast d has
`TV>=|d|/2`, hence `C_vis>=d^2/2`. For the three-state rational control at
times 0,1/2,1, the ES `d>1/3456` gives `C_vis>1/(2*3456^2)`.
These are mathematical diagnostic constants, not physical noise or sample
requirements, and do not refute native HB's supplied Markov law.

## 4. A waiting-time law that pair endpoints do not determine

Conditional on D_0=i, the probability of no **observed** change up to t is
`eta_0^i exp(t L_ii)1`, not the endpoint probability Pr(D_t=i|D_0=i).
Uniformization gives

    survival_i(t)=exp(-lambda t) sum_(n>=0) (lambda t)^n/n!
                     Pr(Y_1=...=Y_n=i | Y_0=i).

For r=3, n=0,1 match iid marks, while n=2 has correction
`epsilon p_i^2 h_i^3/3`. The remaining terms are O(t^3), bounded by the
Poisson tail because all conditional probabilities are in [0,1]. Therefore

    survival_ES,i(t)-exp[-lambda(1-p_i)t]
       =epsilon p_i^2 h_i^3 lambda^2 t^2/6+O(t^3).    (5)

In the rational control at i=0, the quadratic coefficient is exactly 1/96.
This exposes memory in an unmarked holding-time event even though every
two-time endpoint law agrees. It is a statement about this ES family,
not a general test that any matching waiting-time law implies Markovity.
Unobserved null refreshes must remain inside L_ii to obtain the right law.

## 5. What has actually been derived

The source's current conditional criterion now has: an explicit causal
filter from a specified generator/readout/preparation; full-path error
upper bounds; and a nonzero lower-defect witness for the ES counterfamily.
The epsilon=0 countercontrol prevents a false necessity claim. These are
newly paid conditional bridges, not source selection. L, R, mu, p, lambda
and epsilon remain premises of this diagnostic, not a natural nadsoliton
mechanism. No physical field, bath, clock or observer is introduced.
