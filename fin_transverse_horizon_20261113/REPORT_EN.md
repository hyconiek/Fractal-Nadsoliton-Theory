# Extended transverse-amplification horizon

This continues the strongest tested preparation: w=1.25, epsilon=0.1,
background psi phase pi/2, shape phase zero, transverse ring mode j=6.
The full nonlinear RP background and infinitesimal variational columns are
integrated together, with no perturbation back reaction. The equations,
static reference-energy norm and two-dimensional seed-space Gram analysis
are unchanged from the phase/gap campaign.

Three runs use T=200 and periodic [-230,230), dx=0.1: a static control at
dt=0.02, driven dt=0.02 and driven dt=0.01. Increasing the domain preserves
the propagation margin. The coarse driven trace on [0,100] is compared with
the preceding L=130 run at exactly the same dx and dt. This is a domain/
overlap consistency test, not a rigorous boundary-error estimate. Spatial
refinement was performed at T=100 previously; it is not repeated at T=200.
Fourier collocation without explicit dealiasing and RK4 remain shared methods.

Both individual initial seeds retain fixed identities throughout the run:
position seed (f/sqrt(gap),0) and velocity seed (0,f), with
f=sqrt(3/4)sech^2(x), gap=1.25-1.1519772119538363. Their gains are reported
alongside the largest generalized Gram eigenvalue. A maximizing combination
can change over time, so fixed-seed gains help assess whether growth is only
an artifact of re-optimizing initial phase. No supremum over all spatial
perturbations is claimed. Gains measure energy; norm gains are square roots.

RESULTS.md and comparison.json report fixed-time samples, the sampled peak,
static conservation, temporal refinement and domain overlap. Every half-unit
observation is retained in CSV. A peak or endpoint increase on this finite
interval is not a Floquet exponent, asymptotic growth rate or nonlinear
instability theorem. A subsequent decline would not prove global stability.

The model still assumes the conditional A11 coefficients, wave kinetics and
state map. No strict FIN source, AW extension, selector discharge, physical
units or experimental verification is supplied.

## Replay

From the repository root with Python and NumPy:

    OPENBLAS_NUM_THREADS=1 python fin_transverse_horizon_20261113/run.py
    python fin_transverse_horizon_20261113/analyze.py

The analysis checks finite positive Gram spectra, trace lengths, static
conservation, temporal refinement and agreement on the historical overlap.
Historical packages are read only. Input and output hashes are recorded.
