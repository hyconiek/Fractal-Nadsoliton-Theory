# Replay

From the repository root, first regenerate the numerical comparison:

```sh
python3 fin_packet_cascade_research_20261104/cascade_preview.py
python3 fin_packet_cascade_research_20261104/verify_research.py
```

The comparison solver requires NumPy (tested with 2.3.5); the verifier uses
only the Python standard library. Outputs are cascade_preview.json,
odd_spectral_preview.json and verification.json. Historical files are
unchanged. The checks distinguish exact finite controls, sampled spectral
quadratures, numerical comparison diagnostics and input hashes.

The solver uses a periodic Fourier grid and centered time differences.
Its experiment covers m^2=16/3, slow time S=1, R=4 and R=8. Time steps
are halved for both widths; grid and domain refinements are performed at
R=4. Outer-band diagnostics supplement, but do not replace, the domain
check. These finite runs cannot establish an all-R/all-band theorem, and
do not simulate the full nonlinear RP solution. The analytical proof is
separate and remains self-audited rather than formally verified.
