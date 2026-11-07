# The full actual cubic source has a logarithmic critical-time bound

The endpoint estimate now covers all cubic normal-form cascade terms.
Together with TB's direct contribution, it gives the actual squared trace
bound O(epsilon^2 [1+log(1/epsilon)]). The existing CS reduction transfers
this improvement to the integrated squared local-continuum norm and the
whole compact trace. Constants depend on the fixed model coefficients
and critical-window length.

The new estimate first proves that the quadratic multiplier maps H1 x H1
into L1. A bounded distorted test function then controls the low-frequency
momentum conversion. This allows the full physical phase projection to be
retained without requiring an extra derivative on the cubic source.

The complete threshold functional has two phase coordinates, gamma_1 and
gamma_2. The momentum coordinate contains the direct source pairing and
the cascade paired with Omega_psi r0. The position coordinate also carries
a cascade contribution. Both coordinates are included in the proof.

## Numerical comparison

The full signed bilinear sums were evaluated on projected Gaussian
comparison fields with nonzero phase velocities. Results for g4=1 are
shown below. These fields are specified snapshots rather than a simulated
prepared trajectory.

| Grid points | Half-domain | gamma_1 | gamma_2 |
|---|---|---|---|
| 256 | 20 | -0.00232475500596 | 0.00253383029164 |
| 512 | 20 | -0.00232475500596 | 0.00253383029160 |
| 768 | 30 | -0.00232475500596 | 0.00253383029160 |

The cascade momentum contribution is approximately 0.00461283463419.
The fixed-snapshot derivative of gamma_2 with respect to g4 is
-0.00107698077300, agreeing with the direct threshold-pairing formula.
A universal cancellation across arbitrary snapshots and coefficients
would therefore require an additional restriction.

An initial reality check detected the unpaired Nyquist frequency in the
finite Fourier representation. The final implementation excludes that
mode in both inputs and the retained bilinear output. It uses linear
convolution followed by truncation, retaining all sign combinations.
After correction all 40 checks passed, including grid/domain refinement,
phase reality, momentum duality, parameter dependence and source hashes.
The maximum recorded imaginary contamination is below 2e-17.

The analytical result is conditional on the predecessor global and
quadratic-symbol bounds. The logarithmic loss, actual threshold modulation,
modal-action errors and independent proof review remain open. The next
step is a time-dependent estimate of the combined threshold functional
along the actual prepared evolution and of its low-frequency remainder.
