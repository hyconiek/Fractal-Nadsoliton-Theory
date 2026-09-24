# NL-05 — transport baseline after the NL-03/NL-04 results

## Result

`NEGATIVE_BASELINE_QUADRATIC_MODEL_TRANSPORTS_BUT_DOES_NOT_LOCALIZE`.

Because NL-04 did not source a unique kinetic law, NL-05 is explicitly
conditional.  Use the minimal formal Hamiltonian lift together with the
certified qualitative fact that transverse modes are gapped and the NL-02
Dirichlet coupling.  For one transverse normal mode this gives the discrete
Klein-Gordon baseline

`mu h rddot_i + rho h r_i + (kappa/h)(2r_i-r_(i-1)-r_(i+1))=0`, rho>0.

The static quadratic operator is strictly positive (generalized minimum
`rho=1` in the replay), so the only finite-energy decaying static solution on
the infinite line is zero.  Therefore the already certified quadratic content
cannot by itself produce a localized matter-like branch.

It does support propagation.  A deterministic N=512 wave-packet replay shows
spreading on the finite circle: IPR falls from
`0.045222` at t=0 to
`0.017520` at t=8 and peak amplitude from
`1.000` to
`0.432`.  This is transport/dispersion,
not a soliton certificate.

## Consequence

A localized branch now requires a genuinely nonlinear transverse FIN potential
or nonlinear intercell coupling with a separately sourced kinetic law.  It is
not enough to add momentum to the existing quadratic reduction.  The next test
must use the *actual* FIN transverse potential in the complete four-invariant
basis and check whether it has inequivalent minima/topological sectors before
launching continuation or stability numerics.
