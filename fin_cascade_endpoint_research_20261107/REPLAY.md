# Replay

Run from the repository root with Python and NumPy.

```sh
python fin_cascade_endpoint_research_20261107/check_cascade.py
```

The script evaluates every output/input sign of CN's cross-channel
bilinear symbol. The two H1 inputs are transformed on a periodic grid.
Their Fourier coefficients are convolved linearly, truncated to the
retained frequencies and transformed back. The unpaired Nyquist mode is
excluded to preserve real-phase conjugacy. No wraparound convolution is
used. The fixed mass squared is 16/3.

Comparison snapshots use projected Gaussian positions and velocities,
with configurations (N,L)=(256,20),(512,20),(768,30). Verification records
reality, projection, momentum duality, fixed-snapshot g4 dependence,
grid/domain refinement and predecessor hashes. NumPy's version is saved
with the output. The analytical proof supplies the uniform bounds.
MANIFEST.json hashes package files except itself.
