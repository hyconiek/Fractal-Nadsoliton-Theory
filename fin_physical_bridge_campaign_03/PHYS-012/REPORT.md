# PHYS-012 — PUBLIC-SYSTEM-FIDELITY-STRESS-TEST

Date: 2026-09-29

Status: **PUBLIC-DATA NUMERICAL AUDIT**

## Question

Can the existing published coupled Potts-machine output be used as evidence
that a new FIN q=12 experiment will automatically satisfy the frozen
TV<=0.003 hardware envelope?

## Result

No.

The five public 16-channel Potts probability outputs give final 256-state
TV errors relative to their own theoretical targets of approximately:

- `probabilities_29V_illum14.csv`: 0.047338
- `probabilities_29V_illum2.csv`: 0.070508
- `probabilities_31V_illum0.csv`: 0.871177
- `probabilities_31V_illum1.csv`: 0.036161
- `probabilities_31V_illum8.csv`: 0.033901

The best is about **0.0339 TV**, more
than **11.3x** the frozen
0.003 gate.

The notebook accumulates about 4,999,500 samples in the final distribution.
A generic sampling-only upper estimate
`0.5*sqrt(256/n)` is only about 0.0036, so the observed
best-case ~0.034 discrepancy cannot be treated as finite-sampling noise alone.

## Consequence

The full published coupled-Boltzmann implementation must **not** be assumed
accurate enough for FIN.

This does **not** kill the VRSPAD platform. The FIN N=2 test is much smaller
(two q=12 nodes, 144 joint states, one 12-bin symmetry observable), and can be
calibrated directly. But its fidelity has to be measured rather than inherited
from the paper's qualitative statement that ideal distributions are replicated.
