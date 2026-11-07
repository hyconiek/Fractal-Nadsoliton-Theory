# Replay

Run from the repository root with Python and NumPy.

```sh
python fin_forced_packet_interface_research_20261107/check_interface.py
```

The recorded NumPy version appears in verification.json. The script checks
the source spectrum against spatial Gauss-Legendre quadrature, refines the
annular integral and compares the causal time integral with a closed-form
forced oscillator solution. The time-integrated coefficient mass is also
recorded. The analytic bound on atomic annular mass is evaluated at four
widths. Source hashes are checked against INPUTS.json.

The diagnostics verify finite identities and quadrature convergence.
The uniform bounds and causal limit arguments are in the proof document.
MANIFEST.json binds package files other than itself.
