# FIN PHYSICAL BRIDGE — HANDOFF AFTER PHYS-011–014

Date: 2026-09-29

## Main update

The engineered bridge has been sharpened from a one-node categorical shortcut
to a two-node q=12 physical interaction test.

For N=2, the exact FIN joint law is P(i,j) ∝ exp[(g/2)A7[i,j]); the pair
difference d=(j-i) mod 12 has exactly the frozen PHYS-007 histogram. Thus 24
channels and one 12-bin readout implement the existing test without an anchor
or postselection.

The VRSPAD-144 RTL/host design is structurally compatible with arbitrary
contiguous Potts node sizes and independent 144x144 weight rows.

## Negative evidence retained

Published q=4 coupled-machine output distributions show 3–7% TV error in their
better operating records, far above the FIN 0.003 gate. Therefore existing
system-level fidelity cannot be inherited.

## New numerical design result

Using only the public VRSPAD-144 calibration table at 31.5 V / 100 uA and the
public calibration rule, 26 channels meet the predeclared box
|T-27.25|<=1.5 and |baseline residual|<=1%. Freezing 24 of them, compiling
per-destination signed weights w_ab=round(T_a E_ab), and solving the complete
144-state two-node CTMC gives FIN pair-difference TV≈4.13e-4.

All four primary PHYS-007 countermodels remain within the 0.003 numerical
compilation envelope under the same rule, and the closest FIN-vs-FLAT contrast
remains ≈0.01161 instead of 0.01179 ideal.

## Current status

**READY FOR REAL Q12 PRECALIBRATION ONLY.**

Real validation remains closed until a real two-q12 configuration independently
passes its calibration gates. No external spend is authorized here.

Unchanged: task333 frozen, task335 unopened, task337 WIP, A7 source lane
NO_NEW_SOURCE.
