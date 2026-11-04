# Replay

From the repository root:

```sh
python3 fin_threshold_packet_research_20261104/verify_research.py
```

Python standard library only. The script writes verification.json and
packet_preview.csv. It checks Gaussian moments, finite-time Fourier and
distorted-wave quadratures, convergence previews, real-carrier averaging,
scaling controls and pinned source hashes. The sampled error bound used
in the preview is not claimed as the analytical uniform theorem constant.
The analytical proof covers the asymptotics; these are finite numerical
controls, not a nonlinear PDE simulation or formal verification. Historical
inputs and outputs remain unchanged. No packet is added to the RP data.
