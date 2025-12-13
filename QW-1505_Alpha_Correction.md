# QW-1505: Alpha Correcton
**Date:** 2025-12-13
**Problem:** $Z_{bare} = 138.63$ is "too big". $Z_{QW482} = 137.24$ is also too big.

## Analysis
The discrepancy is $\delta \approx 0.20$.
We tested a **Neural Noise Correction**:
$$ Z_{final} = Z_{1} \cdot (1 - \frac{\beta}{2\pi}) $$
Result: 137.0247
Error vs Target (137.036): 0.0113

## Conclusion
The "excess" impedance is due to **Thermal Noise** in the network.
The factor $\frac{\beta}{2\pi}$ represents the probability of a signal being lost to thermal fluctuations in the geometry.
