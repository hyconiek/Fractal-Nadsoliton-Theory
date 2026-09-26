# ROBUST-QUARTIC-DRIFT-OPTIMUM-24 — sign-robust quartic probe against hidden target-convention drift

Status: **RIGOROUS GLOBAL CERTIFICATE IN THE CO-PHASED FOUR-SECTOR CLASS + FULL-7D NUMERICAL PHASE AUDIT**.

Scope: declared stationary heat-bath + ME7 lane at `g_eq=3.7183448981203875`.
No physical apparatus, activity source or clock is inferred.

## 1. Robustness functional

For a retained probe `phi`, write

    m = <phi^2>_u,
    h = P_H(phi^2),
    H = ||h||_u^2,
    C = C_0(phi)

for the intrinsic general-quartic coefficient.  FINITE-N-TARGET-CORRECTION-THEOREM-21 gives, for hidden row-sum drift `b=P_H(B1)`,

    Delta C_B = -3 g m <b,h>_u.

If `||b||_u <= rho`, the worst convention shift is at most

    3 g m rho ||h||_u.

Therefore the sign-preserving radius of a non-null positive probe is

    rho_crit(phi) = C / [3 g m ||h||_u].

This quantity is scale invariant: under `phi -> s phi`, both numerator and denominator scale as `s^4`.  Hence the projective optimal direction is independent of the chosen positive homogeneous normalization (coefficient norm, uniform-Fisher variance, etc.).

For the coefficient-sphere convention `sum c_k^2=1` in the orthonormal real sectors `k=3,4,5,6`, `m=1/12`, so

    rho_crit = 4 C / [g sqrt(H)].

## 2. Certified co-phased optimum

The co-phased four-sector problem was reduced to the simplex `y_k=c_k^2`, `sum y_k=1`.

A strict Krawczyk iteration for the stationary equations contracts to

    y3 in [0.063925514393433802, 0.063925514393436383]
    y4 in [0.22584354095884618,  0.22584354095885201]
    y5 in [0.43854229606623746,  0.43854229606624517]
    y6 in [0.27168864858146646,  0.27168864858148256]

or amplitudes

    c3 in [0.25283495484887725, 0.25283495484888235]
    c4 in [0.4752299874364476,  0.4752299874364538]
    c5 in [0.6622252608185809,  0.6622252608185867]
    c6 in [0.5212376124009725,  0.5212376124009880].

The same interval-AD box gives a strictly negative Hessian.  At radius `5e-4` around the center, the interval perturbation bound gives

    lambda_max(Hessian) < -0.4207,

so the certified stationary point is a strict local maximum throughout that box.

Iterated Krawczyk contraction yields

    rho_* in [0.71753326806113993, 0.71753326806121198].

Two independent outward-rounded global exclusions were then used:

1. interval branch-and-bound on the whole simplex excluding the broad box
   `y3 in [0.02,0.12]`, `y4 in [0.15,0.30]`, `y5 in [0.35,0.55]`;
2. interval mean-value branch-and-bound on that broad box excluding the certified local Krawczyk box.

Both queues exhaust completely.  Thus the value above is a global certificate **inside the co-phased four-sector class**.

## 3. Gain over pure k=5

For a pure unit `k=5` probe,

    C5 = (g lambda5 - 6)/144,
    H5 = 1/288,
    rho5 = sqrt(2) (g lambda5 - 6)/(3g)
         = 0.32290507949561775...

Therefore

    rho_*/rho5 in [2.222118243484856, 2.222118243485080].

The robust mixed probe tolerates more than twice the RMS hidden row-sum drift that can adversarially cancel the pure-k5 intrinsic quartic coefficient.

## 4. Normalization independence

Because `rho_crit` is projective, fixed uniform-Fisher variance selects the same physical ray.  In Fisher-normalized coordinates `w_k=sqrt(lambda_k)c_k`, the same direction has approximately

    w^2 = (0.05532, 0.21917, 0.44475, 0.28076),

and reproduces the same `rho_*`.

Thus the robustness optimum is not an artifact of coefficient-sphere versus Fisher-sphere normalization.

## 5. Full 7D phase audit

The co-phased reduction is rigorously sufficient for the earlier *signal-strength* objective, but it was not automatically sufficient for the ratio `rho_crit`, because phases affect both `C` and `H`.

A direct full-7D optimization was therefore run with independent cosine/sine coordinates for `k=3,4,5` plus the real `k=6` coordinate.  Five differential-evolution campaigns and 100 independent BFGS starts all converged to

    rho = 0.717533268061...,

with the same sector magnitudes as the interval-certified co-phased solution.

All 100 local starts converged to the same value to displayed precision.  At the optimum the sector phases form an arithmetic progression:

    theta4-theta3 = theta5-theta4 = theta6-theta5  (mod 2 pi),

so the pair contributions entering each hidden Fourier channel are co-phased.  The hidden vectors `h` and `P_H(phi A7 phi)` are channelwise phase-aligned.

This is strong full-7D numerical evidence, but a rigorous phase inequality proving that no non-co-phased configuration has larger `rho_crit` is still open.  Therefore the word **global** is presently proof-grade only for the co-phased class.

## 6. Methodological boundary

Promoted:
- exact scale-invariance of the drift-robustness functional;
- rigorous global optimum value and direction in the co-phased sector class;
- rigorous >2.222x robustness gain over pure k5 in that class;
- strong full-7D numerical confirmation of the same optimum.

Still open:
- an interval/SOS proof of the phase reduction for `rho_crit` on the full 7D sphere;
- any preferred physical normalization or apparatus;
- any physical clock, sourced activity, QW-2191, role transfer, `L_total`, SM/GR or ToE conclusion.
