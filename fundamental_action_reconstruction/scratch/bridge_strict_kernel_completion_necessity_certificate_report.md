# Strict Kernel Completion Necessity Certificate

Status: `all-three-explicit-factors-are-necessary-for-exact-pointwise-completion; phase-and-damping-are-shape-critical`

## Result

This finite-domain audit exhausts all subsets of the explicit completion factors
`A(d)`, `P(d)`, and `D(d)` that complete the historical legacy carrier into
the current strict working kernel on `d=0..11`.

## Factor definitions

- `alpha_normalization`: A(d)=1/alpha_geo; global amplitude normalization only
- `phase_frequency_transport`: P(d)=cos(omega_S*d+phi_S)/cos(omega_L*d+phi_L); sign/phase transport
- `damping_compression`: D(d)=(1+beta_tors*d)/(1+beta*d^eta); denominator compression beta_tors*d -> beta*d^eta
- `identity`: K_strict_gate(d)=K_legacy_ont(d)*A(d)*P(d)*D(d) on the audited finite domain

## Necessity summary

- Exact subsets without extra scalar: `['alpha_normalization+phase_frequency_transport+damping_compression']`
- Exact subsets up to one global scalar: `['phase_frequency_transport+damping_compression', 'alpha_normalization+phase_frequency_transport+damping_compression']`
- Missing-phase sign mismatch union: `[2, 3, 4, 5, 8, 9]`
- Minimum best-scalar L2 residual when phase is missing: `5.231383063045e-01`
- Minimum best-scalar L2 residual when damping is missing: `7.523311877932e-01`
- Max quotient identity residual: `1.110223024625e-16`
- Damping factor positive: `True`
- Damping factor strictly decreasing after d=0: `True`

## Subset enumeration

| subset | missing | max residual | best scalar | max residual after scalar | sign mismatches |
|---|---|---:|---:|---:|---|
| none | alpha_normalization,phase_frequency_transport,damping_compression | 2.691540315724e+00 | 5.121500127355e-02 | 8.638519110418e-01 | [2, 3, 4, 5, 8, 9] |
| alpha_normalization | phase_frequency_transport,damping_compression | 1.029220678636e+00 | 1.419981349405e-01 | 8.638519110418e-01 | [2, 3, 4, 5, 8, 9] |
| phase_frequency_transport | alpha_normalization,damping_compression | 2.147744179995e+00 | 1.444810269970e-01 | 5.915168020448e-01 | [] |
| damping_compression | alpha_normalization,phase_frequency_transport | 1.414306364516e+00 | 3.971815651799e-01 | 3.274771722084e-01 | [2, 3, 4, 5, 8, 9] |
| alpha_normalization+phase_frequency_transport | damping_compression | 6.518564876187e-01 | 4.005864660296e-01 | 5.915168020448e-01 | [] |
| alpha_normalization+damping_compression | phase_frequency_transport | 3.405761500938e-01 | 1.101221128299e+00 | 3.274771722084e-01 | [2, 3, 4, 5, 8, 9] |
| phase_frequency_transport+damping_compression | alpha_normalization | 1.749236466809e+00 | 3.606737602222e-01 | 1.110223024625e-16 | [] |
| alpha_normalization+phase_frequency_transport+damping_compression | none | 2.220446049250e-16 | 1.000000000000e+00 | 5.551115123126e-17 | [] |

## Guarded interpretation

The audit shows which explicit completion factors are necessary for the finite
Z12 completion identity.  It does not prove that strict nadsoliton dynamics must
generate those factors, does not prove `beta_tors -> chi_11`, and does not
transfer legacy physical roles onto `K_strict_gate`.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), or D(d) from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
