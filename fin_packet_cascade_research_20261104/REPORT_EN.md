# Retarded quadratic response with full packet tails retained

PC proves a finite-window local-energy bound for the odd continuous
response to the entire finite packet source -6Pc,o(K P_R^2). The packet
has width R and the time interval is S R^2. The integrated squared local
phase norm is bounded independently of R for fixed S.

The proof uses an exact quadratic normal form. Its positive-weight source
is O(R^-1), so the odd propagator's integrable decay controls its full
history. The initial correction is retained explicitly and controlled by
an all-time odd-continuum local-energy estimate from unweighted energy
data. No negative-weight local source error is inserted into a positive-
weight Duhamel estimate, and no spatial tail or outgoing wave is removed.

At packet amplitude epsilon the comparison response's squared integral
is O(epsilon^4), including when R=epsilon^-1. This is a packet-driven
linear comparison theorem, not a bound for the actual continuum generated
by the original nonlinear preparation. The actual O(epsilon^2) cubic trace
and modal action-law errors remain open.

Next: estimate the feedback of this odd response into the even packet
history, with the broad spatial product and initial transient retained.

## Numerical comparison

At m^2=16/3 and S=1, the finite-grid results were:

| R | dt | dx | Half-domain | Integrated local phase square | Relative energy/work defect |
|---|---|---|---|---:|---:|
| 4 | 0.02 | 0.0625 | 64 | 0.1202419965 | 0.0013060 |
| 4 | 0.01 | 0.0625 | 64 | 0.1202416189 | 0.0003267 |
| 4 | 0.01 | 0.03125 | 64 | 0.1202416189 | 0.0003267 |
| 4 | 0.01 | 0.0625 | 96 | 0.1202416189 | 0.0003267 |
| 8 | 0.02 | 0.0625 | 128 | 0.1313986707 | 0.0006113 |
| 8 | 0.01 | 0.0625 | 128 | 0.1313933086 | 0.0001527 |

The energy/work defect decreased by approximately four under time-step
halving. At R=4, spatial and domain refinement left the displayed local
integral unchanged at the shown precision. These tests cover two widths
and one mass/slow-time choice; they do not establish the analytical
uniform bound or simulate the full nonlinear trajectory.
