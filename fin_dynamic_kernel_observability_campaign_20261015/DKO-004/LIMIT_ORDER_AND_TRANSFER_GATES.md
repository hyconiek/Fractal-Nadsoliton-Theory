# DKO-004 — Mandatory limit order and resource interpretation

1. Fixed α∈(0,1), fixed finite fast horizon T: N→∞ from the microscopic labelled FIN model to the Gaussian OU observation process. Theorems DKO-001–003 live here.
2. Inside the **already obtained Gaussian family** only: fix h>0 and analyze α↓0, M↑∞. KL=Θ(Mα²) and sign-test achievability hold for this analytic family.
3. The formula M≈(N/k)² is a purely **formal Gaussian-family analogy** if α=k/N. It is **not** a valid finite-N FIN requirement because taking α_N↓0 and M_N→∞ simultaneously could break the snapshot CLT, tightness, soft freezing and MT post-relaxation residence on expanding horizons.
4. To promote such a statement one would need a uniform triangular observation/process theorem controlling α_N, T_N=M_Nh, path exits, soft-coordinate drift, non-Gaussian conditional block terms and calibration uncertainty. This package does not supply that theorem.
5. Nothing here estimates physical measurements, sensor count, seconds, iid replicates, finite sample N0, physical kernel source or externally calibrated observation errors.
