# Corrected kinetic and memory notation

Old C03 fields `M0_coupled = 0.00995201522562` and `M1_derivative = -0.018039981780484` correspond to the scalar equilibrium quadratic channel for the UNNORMALIZED indicator of visible Y00: **K(0)** and **K'(0)**, respectively, not integrated M0/M1.

Here we use the normalized resolved basis vector `e00` (orthonormal in the reversible coordinates), so `K(0)=0.160181620467715` and `K'(0)=-0.290360640462747`. Multiplying by `piY00=0.062129570150189` gives the previous unnormalized channel values (up to numerical rounding).

The actual integrals of the scalar memory channel computed by a **numerical CG solve** are `M0=∫ K(t)dt ≈ 0.0916819627335244` and `M1=∫ t K(t)dt ≈ 0.0547291733152583`. These are not interval-certified global integrals; no independent rigorous Q-sector spectral-gap/tail certificate is supplied.

For all claims S is the parent **M=3** reversible generator. The M=4→M=3 static marginal is a different experiment and is NEVER substituted for piY.
