# BA audit: exact parity cuts, surviving microscopic memory and proof scope

Date:2026-10-16. Analytical finite/model proof, not formal mechanization
or external peer review. The inherited exact MT/UF/KR premises remain explicit.

## 1. The new cut theorem is correct

Let c_xy=pi_x q_xy=c_yx and d(x)=Lm(x), with parity jumps0 or +/-delta,
delta=2/N. The exact statewise constraint is

    sum_y c_xy(m_y-m_x)=pi_x d(x).

Sum over A_l={x:m_x<=a_l}. Internal edges cancel in pairs. All edges crossing
this cut change m upward by exactly delta, so

    delta K_l=sum_(x in A_l)pi_x d(x).

This requires the actual nearest-parity-level support, not general multiple
copy jumps. In the MT domain the permitted count graph is connected, so
occupied parity levels cannot skip an intermediate lattice value. Diagnostic
rates may nevertheless vanish; no diagnostic irreducibility follows.

With common pi and exact native d, all K_l are exactly native. Conditional
instantaneous up/down rates are K_l/p_l and K_(l-1)/p_l. Consequently all
bilinear Dirichlet pairings on parity functions agree, and

    E[G|m=a_l]=N delta^2(K_(l-1)+K_l)/p_l
               =E[G_native|m=a_l].

The factor N and both incident cuts are necessary. No division by a zero
cut is needed. The formula is a compression PLP, not a strong-lumpability
or observed-semigroup theorem. It does not determine individual microscopic
conductances or within-level motion.

## 2. The asymptotic projection improvement is paid

Native HB has E_piD|G_native-2|->0 by continuity at u, MT concentration,
bounded bracket and vanishing reflected departure mass. The exact conditional
identity and Jensen give

    ||P_m(G-EG)||_1<=||G_native-EG_native||_1 ->0.

The uniform cap on G implies the same projected L2 convergence. This
strengthens KR's smooth weak conditional test to exact finite-N conditional
rigidity plus a strong asymptotic bound. It is not NC for the full G.

The remaining noise lies in the fibre residual r_perp=(I-P_m)(G-EG).
Stationarity gives directly

    E sup_(t<=S)|integral_0^t P_m r(X_sqrt(N)s)ds|
       <=S||P_m r||_1 ->0.

Thus IBA for the whole residual is equivalent to IBA for r_perp, up to this
vanishing error. Neither P_m nor the true invariant projection P_kerL may
be silently replaced by the other.

## 3. Component and Poisson statements are correctly qualified

On an actual closed communicating component C, invariance of m^2 gives
E_C G=-2N E_C[m d]. The source formula for h=P_kerL(G-EG) is exact.
Having all native domain edges available is not the same as assigning
them positive diagnostic rates. If all rates are positive, connectedness
does give h=0, but no quantitative lower bound or Poisson estimate follows.

The finite forward/backward martingale bound is already admitted in the
KR intake. Reuse its deterministic/rational-time endpoint convention for
jump reversal, its componentwise inverse, and its sqrt(N) time factor.
Re-deriving it is useful verification, not a new proof of the unpaid class
estimate. A large Poisson cost alone cannot refute IBA; a small full gap is
not required when the relevant residual has small slow-mode overlap.

## 4. Independent validation of the six-state coefficient

The source has a six-node cycle, not an all-N q12 count construction.
Its common drift and parity cut traffic are preserved; at its disconnected
endpoint, the component residual means are genuinely nonzero.

We independently check the Poisson costs using cycle currents instead of
the producer's matrix inverse. For edge i->i+1 put J_i=c_i(u_i-u_(i+1))
and f_i=pi_i r_i. The Poisson equation is J_i-J_(i-1)=f_i. Thus

    F_i=sum_(k=0)^i f_k, J_i=j+F_i,
    j=-sum_i F_i/c_i / sum_i 1/c_i,
    R=sum_i J_i^2/c_i.

As tau->0, the small edges are0,1,4 with c_i=tau. All F_i converge;
j tends to -(F_0+F_1+F_4)/3, and large-edge contributions to tau R vanish.
The exact small-edge current squares sum to1064/126075. This supplies a
direct rational proof of the reported leading coefficient, not a fit to
three sample points or an unproved coarse-network analogy.

The producer's24 flags are not all strict exact-arithmetic decisions:
one monotonicity flag compares binary64 previews. The intake validates the
same ordering with exact Fractions and validates the limit independently.
Do not relabel the original implementation as24 independent exact proofs.

The toy has a nearest-level graph and arbitrary rational pi, not the actual
q12/A7 MT count measure, statewise native drift, or an all-N critical sequence.
Its obstruction therefore does not refute universal q12 IBA.

## 5. Current admission and genuinely new next question

Accept exact finite parity-compression rigidity, strong projected-noise
convergence and the correctly scoped finite negative control. Preserve
OPEN_Q12_UNCONDITIONAL_IBA and the KR sufficient process theorem.

The current full AC/rate/KR/MT/UF review chain is locally bound, not
retrospectively read by the source producer. No physical activity provider,
initialization law, clock/units, locality/conservation, EV or ToE is supplied.

The accompanying count-partition seed narrows a possible counterexample to
one precise finite conductance-feasibility question. Its static properties
and conditional obstruction are derived; admissible generators closing its
two parts have NOT been constructed. Do not promote that candidate or a
dual infeasibility result for it into a universal scientific verdict.
