# Replay

From the repository root:

```sh
python fin_even_packet_feedback_research_20261105/feedback_preview.py
python fin_even_packet_feedback_research_20261105/verify_research.py
```

The solver uses NumPy (recorded run: 2.3.5). The verifier uses the Python
standard library. The solver runs six explicitly listed periodic Fourier
space / centered-time comparisons. The packet is evaluated from its exact
Darboux/free Fourier formula at every time. All driven responses start at
zero, with the second-order initial time step. Bound projections and parity
are applied to the sources. The three even channels retain cross terms in
the trace Gram matrix: use v=(1,1,g4), and v^T G v for any fixed signed g4.
The coefficient-one run is a basis recombination, not parameter fitting.

Energy/source-work defects, time, grid and domain refinement are numerical
diagnostics. They are not interval certificates or evidence of a uniform
large-R theorem. Read the proof for the analytical result and its assumptions.
The manifest hashes package files except itself; INPUTS pins predecessor
files. Verification does not regenerate or validate an analytical proof.
