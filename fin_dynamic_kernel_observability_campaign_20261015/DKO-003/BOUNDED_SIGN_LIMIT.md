# DKO-003 — Bounded two-time sign theorem

Freeze α=1/2, h=1, fixed r∈{1,3,4}, model IDs A7 and A7+0.01(P3−P4). Set sign(0)=0 at finite N and

H_{r,N}=sign(V_{r,N}(t_N)) sign(V_{r,N}(t_N+1)).

By DKO-002's fixed two-time weak convergence, the limiting Gaussian pair (V_r(0),V_r(1)) is centered, nonsingular, with variance C_r(0)=γ_r^{-1}+1 and covariance C_r(1)=γ_r^{-1}e^{-γ_r}+e^{-1}. Nondegenerate Gaussian marginals assign zero probability to either sign boundary. By bounded almost-everywhere-continuous mapping, E H_{r,N}→s_r=(2/π) arcsin(ρ_r), where ρ_r=C_r(1)/C_r(0).

**Proof of the arcsin identity (self-contained):** let (X,Y) be a standard bivariate normal with correlation ρ. Write X=G, Y=ρG+√(1−ρ²)H, where G,H are independent isotropic standard normals. The vector (G,H) has a uniform angular coordinate. The two sign rays G=0 and ρG+√(1−ρ²)H=0 cut two disagreement wedges of total angular fraction arccos(ρ)/π; hence P(XY>0)=1−arccos(ρ)/π=1/2+arcsin(ρ)/π. Since sign(X)sign(Y)=2 1_{XY>0}−1 almost surely, its expectation is 2 arcsin(ρ)/π.

For the frozen mode3 ε=+.01, Λ3 increases, γ3 decreases, and at h=1 ρ3 increases strictly (see gauge theorem). Thus E H_{3,N} and P(V3(0)V3(1)>0) have strictly larger limiting values for the control. This is an ideal ensemble law difference for a bounded two-time event, not an unbounded covariance estimate nor a single-record identification guarantee. The zero finite-N outcomes are retained and vanish only in the limiting Gaussian marginal, not discarded by definition.
