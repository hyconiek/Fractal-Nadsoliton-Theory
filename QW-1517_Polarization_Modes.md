# QW-1517: Gravitational Wave Polarization Modes

**Date:** 2025-12-17 02:26:01.565256

## Theory

GR predicts two polarization modes for gravitational waves:
- **Plus (+)**: Quadrupole pattern aligned with axes
- **Cross (×)**: Quadrupole pattern rotated 45°

## Method

- 2D grid simulation with rotating quadrupole source
- Source injects h+ ∝ cos(2ωt) and h× ∝ sin(2ωt)
- Detectors at cardinal directions and 45°

## Results

| Detector | θ+ amplitude | θ× amplitude |
|----------|--------------|--------------|
| North | 0.003585 | 0.004346 |
| East | 0.003585 | 0.004346 |
| South | 0.003530 | 0.004333 |
| West | 0.003530 | 0.004333 |
| NE (45°) | 0.003357 | 0.004154 |
| SE (315°) | 0.003265 | 0.004107 |

## Quadrupole Pattern

- North-East correlation (θ+): 1.0000
- 45° amplitude (θ×): 0.004154

## Verdict

✅ + polarization shows quadrupole pattern
✅ × polarization detected at 45°
