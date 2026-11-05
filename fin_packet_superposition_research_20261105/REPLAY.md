# Replay

From the repository root:

```sh
python fin_packet_superposition_research_20261105/superposition_preview.py
python fin_packet_superposition_research_20261105/verify_research.py
```

The solver uses NumPy (recorded version in JSON), with the EF triangular
scheme adapted to complex Gaussian coefficients. R=4, S=1, m^2=16/3,
half-domain 64, dx=0.125; dt=0.02 and 0.01. The two atoms have relative
widths 1 and 1.5, phases 0 and pi/3, and mixture coefficients 0.6 and 0.4.
Each time step evolves the full quadratic source of the mixture, its shape
response and all three even channels. g4=1 is a fixed basis example.

The diagonal-only comparison is 0.6^3 Z[p1]+0.4^3 Z[p2]. Its difference
from Z[0.6 p1+0.4 p2] isolates the omitted mixed response. Doubling the
input checks cubic homogeneity. This check is supplementary: the energy
work diagnostic and time refinement test the numerical evolution.
This package reuses EF's method; EF has separate grid/domain diagnostics.
The present finite runs do not prove the infinite-series or uniform-in-R
claims, which rely on the analytical estimates. The standard-library
verifier checks bookkeeping, finite diagnostics and predecessor hashes.
