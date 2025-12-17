# QW-1527: Simulation Verification Output

## 1. Frequency Analysis
- **Omega (Orbital):** 0.4000 rad/s
- **Q_xx Amplitude:** 62.4998
- **d2Q/dt2 Amplitude:** 39.9468
- **Frequency Ratio (d2Q/Q):** 0.6392
- **Expected Ratio (2*Omega)^2:** 0.6400
- **Status:** ✅ FREQUENCY CHECK PASSED (Emission at 2*Omega)

## 2. Chirp Mass Scaling Analysis
| Mass 1 | Mass 2 | Chirp Mass (Mc) | Normalized Amp |
|---|---|---|---|
| 10.0 | 10.0 | 8.706 | 73.583 |
| 20.0 | 20.0 | 17.411 | 233.610 |
| 30.0 | 30.0 | 26.117 | 459.174 |
| 10.0 | 30.0 | 14.651 | 175.208 |
| 5.0 | 40.0 | 11.220 | 112.308 |

- **Power Law Exponent (Fitted):** 1.6667
- **Theoretical Exponent (GR/FIN):** 1.6667
- **Status:** ✅ CONFIRMED
