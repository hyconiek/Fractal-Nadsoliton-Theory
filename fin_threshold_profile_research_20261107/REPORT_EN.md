# Conditional closure of the critical cubic trace budget

The full actual cubic trace now obeys the integrated squared bound
O(epsilon^2) on each fixed critical interval S epsilon^-2. The CS reduction
transfers this result to the whole local-continuum squared integral and
compact trace. This removes the logarithmic loss in CE's upper bound.

The argument uses the complete physical cubic phase source in the form
(A_C, Omega_psi B_0+C_D), with all three preimages bounded in L1 by
C epsilon^3. CE supplies the cascade bounds and TB supplies the direct
one. The low-frequency observation kernels are bounded uniformly in space
by C(1+t)^-1/2. The momentum kernel retains the exact order Omega_psi Wcal.
Its bound follows by applying Omega_psi^-1 to two spatial derivatives
of the distorted test function. Integration over the actual source history
then gives the required O(epsilon^2) squared trace. FH controls the fixed
high-frequency output using the global phase norm.

The time-dependent threshold value itself also meets O(epsilon^2) using
only its amplitude bound. A spectral-difference criterion is derived as
an alternative sufficient condition. A translated-source family shows
that the criterion can grow logarithmically even with uniform source
norms. The spatially uniform kernel argument controls that family without
requiring the stronger criterion or a zero threshold coefficient.

## Finite diagnostics

A resonant free comparison source has fixed Gaussian spatial profile and
time coefficient epsilon^3 cos(m t), with m^2=16/3. The exact spectral
oscillator response is integrated numerically against the compact
observation. The Gaussian comparison and the smooth compact cutoff used
in the proof are explicitly distinguished in the replay instructions.

| epsilon | Critical time | Squared trace divided by epsilon^2 |
|---|---|---|
| 1/8 | 64 | 0.22220637 |
| 1/16 | 256 | 0.27419491 |
| 1/32 | 1024 | 0.30377360 |

At epsilon=1/32, halving the time step changes the squared trace by about
5.3e-7 relatively. Increasing frequency quadrature from 1536 to 2048 nodes
changes it by about 4.3e-9 relatively. Symmetric source translations of
16, 64 and 256 also satisfy the finite comparison bound. The weighted
spectral-difference comparison has a measured final logarithmic slope
1.5165, consistent with its analytically derived limiting coefficient 1.5.
All 36 recorded checks passed, including source provenance.

## Scope and continuation

The result is a self-audited conditional PDE proof using the predecessor
global small-energy bounds, quadratic-symbol bounds and exact source
decomposition. Independent review and formal verification are outstanding.
The numerical diagnostics provide finite checks of representative kernels
and source histories.

The next obligation from VD is to construct actual modal/continuum
correctors in a specified chart and estimate the integrated action-work
error as o(epsilon^2), retaining all three radiative channels and mixed
terms. This is required before admitting an actual critical damping law.
