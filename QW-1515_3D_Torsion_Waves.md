# QW-1515: 3D Torsion Wave Simulation

**Date:** 2025-12-17 02:21:08.094299
**Status:** PARTIAL

## Parameters
- Grid: 64³ = 262144 nodes
- Wave speed: c_tors = 1.5109
- Source frequency: f_gw = 0.2

## Results

| Radius | Amplitude |
|--------|-----------|
| 3 | 0.000542 |
| 5 | 0.000275 |
| 7 | 0.000254 |
| 9 | 0.000261 |
| 11 | 0.000237 |

## Amplitude Scaling
$$A(r) = \frac{0.0009}{r^{0.59}}$$

**Expected (GR):** n = 1.0
**Measured:** n = 0.59
**Error:** 0.41

## Verdict
🟡 PARTIAL: n = 0.59 (close to 1.0)
