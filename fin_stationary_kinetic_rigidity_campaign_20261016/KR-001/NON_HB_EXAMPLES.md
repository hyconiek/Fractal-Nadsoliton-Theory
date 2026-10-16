# KR-001 — A genuine q12 non-HB member

Let `q^0_N(n,n−e_i+e_j)` denote each native reflected HB count rate. For `i≠j`, apply the **same-parity thinning**

`q^eta_N(n,n−e_i+e_j) = eta(g_N) q^0_N(n,n−e_i+e_j)` if `sigma_i=sigma_j`; otherwise `q^eta_N=q^0_N`, with `eta(g)=1/(1+g²)`.

Each permitted directed edge and its reverse receive the same strictly positive factor, so detailed balance and pi_D,N are unchanged. Thinning cannot violate the native count jump rate cap or the original per-label rate-one attempt bound. It is the identity at g=0, and for positive g differs on available same-parity edges. Every altered move has `Delta m=0`, whence `L^eta_N m=L^0_N m` **at every reflected state** and `Gamma^eta_N(m,m)=Gamma^0_N(m,m)` pointwise. Native smoothness and MT concentration give NC for this example. Thus the class is nonempty and not restricted to HB.

This changes same-parity kinetics but does not show universality of arbitrary opposite-parity activity variations. It neither sources nor selects a physical rate rule.

## Independent rational negative control

`verify_kr.py` builds a q4/N2 finite reversible count network with ten states and a six-edge conductance stress. It keeps the stationary measure and the **exact pointwise parity drift** while changing parity quadratic variation in two states (`max |Delta Gamma|=1/8`), with positive rates and count attempt capacity ≤N. This is an exact-arithmetic algebraic fixture, **not a q12/A7 stationary critical process**, so it cannot prove failure of NC in the FIN sequence.
