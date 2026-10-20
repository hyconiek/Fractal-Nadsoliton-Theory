# CPL-003: the same static law and compressed drift do not give the same record

Date: 2026-10-20. Ideal mathematical observation, not physical readout.

## 1. A bounded two-time contrast, without ideal quadratic-variation readout

Compare the stationary canonical quartic diffusion
L0=b partial_z+partial_z^2 with the scalar record of the full CPL process.
Both have marginal rho, the same actual statewise soft drift b, and the
same compressed operator E[L f(Z)|Z]=L0 f. Therefore all one-state means
and first Dirichlet pairings on Z agree. This is the same distinction
between compression and a semigroup already paid in BA/HM, reused here.

Freeze f(z)=exp(-z^2), and let C(t)=E[f(Z0)f(Zt)]. In both laws f is in
the L2 generator domain: f is bounded and the actual Lf=b f'+g_P f''/2
is bounded, so Dynkin's formula and strong continuity yield the L2
derivative. No assumption f in Dom(L^2) is needed for the following scalar
correlation. Self-adjointness and the spectral Taylor bound give

    C(t)=||f||^2+t<f,Lf>+(t^2/2)||Lf||^2+o(t^2).

Equivalently the second forward difference equals
||exp(tL)f-f||^2; this separately checks the domain/sign convention.
The residual (L-L0)f=(g_P-2)f''/2 has zero conditional mean given Z. Hence

    C_CPL''(0)-C_0''(0)
        =(1/4)E[(g_P-2)^2 f''(Z)^2]>0.              (1)

Strict positivity follows from g_A<2<g_B and f''=(4z^2-2)exp(-z^2),
nonzero on an open set of positive measure. Thus the BOUNDED observable
f(Z0)f(Zt), valued in [0,1], has a strictly different expectation for all
sufficiently small positive t. No numerical lag threshold or finite-N
error rate is claimed. The generic compression identity is not new; its
application to this fully realized critical-state process is the new test.

## 2. A separate ideal continuous-path singularity

For any fixed S>0 define quadratic variation by the limit of squared
increments on the dyadic grid. It exists almost surely under both laws.
On compact stops bounded drift/diffusion give a mean-square error O(2^-n) for
the realized-square sum; summability and localization prove this claim.
The canonical process has Q_S=2S. The CPL process has

    Q_S=integral_0^S g_P(Z_s)ds<2S on A, >2S on B.

The strict signs persist on every finite continuous path. Thus the event
Q_S=2S has probabilities 1 and 0, and the two ideal scalar PATH laws have
TV exactly one despite (4) and their common one-time Gibbs marginal.
This is not a finite-resolution cost theorem or an empirical discrimination
claim; (1) supplies a different bounded, finite-time contrast independently.
No full native count-path law has been compared or refuted here.
