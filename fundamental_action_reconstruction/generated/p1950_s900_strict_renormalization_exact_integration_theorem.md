# P1950 S900 Strict Renormalization Exact Integration Theorem

Target: `STRICT_B1_ONE_LOOP_DIVERGENCE_EXACT_ALGEBRAIC_CANCELLATION`

We define strict kernel

`K_strict(d)=cos(omega d+phi)/(1+beta d^eta)` with
`omega=0.18575, phi=0.1625, beta=1.0, eta=1.8`.

B1 proxy UV moments on `[epsilon,1]`, `epsilon=1e-6`:
- m0=5.8585640482543466e-01
- m1=2.2265838659866413e-01
- m2=1.2596702576507912e-01

Solve exact linear cancellation system `A a = b` for counterterm coefficients
`a=(a_R2,a_Ric2,a_Riem2,a_GB)`.

Solution:
- a_R2=9.4736842105263330e-01
- a_Ric2=2.1052631578947498e-01
- a_Riem2=-5.2631578947369473e-02
- a_GB=5.2631578947369320e-02

Residual norm `||A a - b||_2 = 6.592e-16`.
Termwise max residual `max_i |(A a - b)_i| = 5.551e-16`.
Residual tolerance `3.270e+00` is computed from backward-error and Simpson-refinement bounds.
Quadrature relative discretization bound `1.615e-01`
against stability tolerance `1.000e-06`.

Direct backend profile source:
`p1848_direct_gravity_operator_profiles_B1`.

Backend profile expressions:
- R2: `cos(0.18575*d + 0.1625)**2/(1.0*d**(9/5) + 1)**2`
- Ric2: `3.24*(-d**(4/5)*cos(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1)**2 - 0.103194444444444*sin(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1))**2`
- Riem2: `41.9904*(d**(8/5)*cos(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1)**3 + 0.103194444444444*d**(4/5)*sin(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1)**2 - 0.00532454668209876*cos(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1) - 0.222222222222222*cos(0.18575*d + 0.1625)/(d**(1/5)*(1.0*d**(9/5) + 1)**2))**2`
- GB: `-12.96*(-d**(4/5)*cos(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1)**2 - 0.103194444444444*sin(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1))**2 + 41.9904*(d**(8/5)*cos(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1)**3 + 0.103194444444444*d**(4/5)*sin(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1)**2 - 0.00532454668209876*cos(0.18575*d + 0.1625)/(1.0*d**(9/5) + 1) - 0.222222222222222*cos(0.18575*d + 0.1625)/(d**(1/5)*(1.0*d**(9/5) + 1)**2))**2 + cos(0.18575*d + 0.1625)**2/(1.0*d**(9/5) + 1)**2`

Gram rank: `3` of `4`.
Gram nullity: `1`.

Verdict: `OPEN_OBSTRUCTION_WITH_TRACE`.
Fail trace: `operator_profile_gram_rank=3 < 4; nullity=1; condition_number=rank_deficient`.

Interpretation: the direct tensor-profile export is now consumed by the same
witness+gate pipeline.  If the Gauss-Bonnet tensor identity makes the B1 scalar
profile Gram matrix rank-deficient, the local cancellation residual may still
be small, but the four-channel coefficient identification is not theorem-grade
unique on this B1 scalar projection.
