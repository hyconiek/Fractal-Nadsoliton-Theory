# Replay

From the repository root:

```bash
python3 fin_rate_coupling_derivation_20261016/verify_derivation.py
sha256sum -c fin_rate_coupling_derivation_20261016/INPUTS.sha256
```

From this directory:

```bash
sha256sum -c MANIFEST.sha256
```

Only the Python standard library is required. The script writes
CERTIFICATE.json and RESULTS.json in this new directory. Historical files
are read only. No simulations, broad scans, installations or network calls
are needed.

The certificate uses exact Fraction intervals throughout. Every fifth-root
endpoint is checked by integer powers, and every cosine uses two consecutive
alternating sums with decreasing terms. The discrete Fourier coefficients
use exact algebraic trigonometry with a bounded sqrt(3). All public decimal
endpoints are rounded outward using integer division. The freshly constructed
active spectral boxes are also checked against the inherited rational boxes.
No uncertainty is inferred from floating-point precision.

The 576 rational algebra fixtures check direct Gibbs finite differences
against the claimed coupling, including nonconstant diagonals and fixed
slot-dependent fields. The eight rate fixtures use ordinary binary64 and
are labelled separately. They do not prove stochastic process convergence,
physical sourcehood or measured error bounds.

The analytic proof, finite exact algebra and numerical checks are different
evidence layers. A reproducible script is not peer review or mechanized proof.
The SHA manifests bind actual local inputs and this package, not missing
historical uploads or independent physical custody.

The observed repository HEAD at entry was
761e0731cb4e45669605611ef0c3b38b58a0eff8. Prior kinetic-composition research
changes were already present and were preserved. This task changes only
`fin_rate_coupling_derivation_20261016/`; it makes no commit or AGENTS amendment.
