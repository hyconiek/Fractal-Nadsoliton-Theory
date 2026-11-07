# Threshold bands and the improved direct-cubic estimate

The actual direct cubic component now has integrated squared compact trace
bounded by O(epsilon^2 [1+log(1/epsilon)]) on the critical interval.
This improves its previous generic O(epsilon) bound. The proof uses the
actual global H1 bounds, the L1 size of cubic products and an explicit
L1 mapping estimate for the Darboux adjoint.

The full converted cubic forcing also has a precise sufficient spectral
budget. Its squared energy density is weighted by epsilon^-1 in the
threshold core and by |k|^-1 in intermediate bands. The square root of
this combined integral is then integrated in source time. A bound of
O(epsilon) for that quantity would close the cubic trace gate. Summing
bands inside the homogeneous Plancherel estimate avoids the additional
loss from a separate triangle inequality over every band.

A bounded spectral source density gives the logarithmic upper bound.
A positive power of vanishing at threshold would remove that loss.
For the direct source, a sufficient condition is zero pairing with the
threshold profile r0 together with a controlled first spatial moment.
These conditions require further analysis along the actual trajectory.

The complete normal-form cascade still needs endpoint control. The next
step is to establish its low-output Fourier bounds and evaluate the
combined direct-plus-cascade threshold coefficient, retaining all
projection and commutator terms.

## Finite verification

All 44 checks passed. Twelve model budgets confirm the factor-two
comparison between the dyadic square sum and the continuous spectral
weight. A comparison snapshot with projected even and odd fields checks
the Darboux commutation and threshold-pairing identities at 2048 and 4096
grid points. Relative commutation errors are below 4e-16. Its threshold
coefficient is approximately -0.00207900434 at both resolutions, showing
that parity and bound-mode removal permit a nonzero coefficient.

The snapshot supplies an algebraic control. The source estimate TB9 is
analytical and applies conditionally along the actual trajectory through
the predecessor global bounds. The remaining logarithmic loss, cascade
estimates and modal-action errors are open. Formal verification and
independent proof review are outstanding.
