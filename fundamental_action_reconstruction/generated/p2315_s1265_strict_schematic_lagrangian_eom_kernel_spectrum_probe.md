# P2315/S1265 — schematic Lagrangian EOM and strict-kernel spectrum

- Status: `OPEN_PARTIAL_SCHEMATIC_EOM_DERIVED_KERNEL_SPECTRUM_NONSELECTOR`
- Schematic EOM derived: `True`
- Full EOM derived: `False`
- Full Task-3 Lagrangian derived: `False`
- Strict kernel pair-plane degeneracy verified: `True`
- `4 ln 2` scalar scaling preserves degeneracy: `True`
- G1/G3 update: `OPEN_UNCHANGED`
- Theorem fingerprint: `a2e28afa1ade855d51139659b5457b01db7ffbf8bdf74c0d83f4ed9bf760a0f2`

## Concrete result
P2315 derives the normalized schematic Euler-Lagrange row from the Release-7 candidate Lagrangian and computes the frozen strict Z12 kernel spectrum.  The operator is circulant, so Fourier pair planes remain degenerate.  A scalar `4 ln 2` multiplier rescales eigenvalues but does not pick an orientation.

## Guardrail statement
P2315 does not export a full tensor EOM, a full theorem-grade Lagrangian, a selector source, QW-2191 discharge, or ToE closure.

## Next honest step
To get a real closure-relevant advance, add one missing nonlinear/tensor object: an explicit potential V with stability theorem, a full tensor-resolved EOM variation, or an internal selector/orientation source that breaks the Fourier pair-plane degeneracy before replaying P2302/P2281.
