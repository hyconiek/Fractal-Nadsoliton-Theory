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
- a_R2=9.9999999999999989e-01
- a_Ric2=0.0000000000000000e+00
- a_Riem2=0.0000000000000000e+00
- a_GB=0.0000000000000000e+00

Residual norm `||A a - b||_2 = 1.388e-16`.
Termwise max residual `max_i |(A a - b)_i| = 1.110e-16`.

Therefore cancellation is numerically exact at machine precision under strict-lane assumptions.
