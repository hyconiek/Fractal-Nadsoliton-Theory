# Three-channel work transfer for actual quadratic modal histories

The zero-data quadratic continuum responses driven by the actual canonical
products Z^2, C^2 and ZC now have a quantitative action-work approximation.
Their work equals CD's three-channel dissipative polynomial integrated
along the actual modal actions, with uniform error
O(epsilon^3 sqrt(1+log(1/epsilon))) through S epsilon^-2. This error is
o(epsilon^2).

The argument retains the nonzero initial envelopes, their switch-on
transients and their entire retarded history. MA's L2 bounds on envelope
derivatives control the history correction. All nonzero quartic carrier
frequencies are uniformly separated from zero, giving an additional
O(epsilon^4) averaging error. The mixed sum channel contributes equally
to the two action-work equations. The closed difference response and
cross-frequency terms are included in the calculation.

The actual action equations are now expressed as the dissipative
polynomial, the above small error, and explicit remainder work D_I,D_J.
Those last terms pair the modes with the difference between the whole
actual continuum and the constructed quadratic responses. Their current
absolute bound is O(epsilon^2). Establishing a little-o bound for that
specific work pairing is the next obligation before admitting damping.
The exact residual source equations and the two unpaid work integrals
are stated in RW8--RW10.

## Verification

All 42 checks passed with SymPy 1.14.0. Exact Laurent-polynomial averaging
reproduces the self-channel coefficients -9 d_a a^4 and -9 d_b b^4,
and the mixed coefficient -18 d_x a^2 b^2 in each action. It also checks
the reference energy-work normalization, zero means of closed and real
susceptibility contributions, and the complete nonzero frequency inventory.
A rational lower bound verifies the uniform carrier gap.

Finite tail-budget controls keep both the initial transient and envelope
history. Their normalized bounds decrease with epsilon consistently with
the analytical little-o estimate. Source hashes are verified. These are
algebraic and bound checks. The theorem is a self-audited conditional
analysis with independent review and formal verification outstanding.
