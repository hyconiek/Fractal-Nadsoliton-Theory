# SRN-001 — Exact common-scale gauge and the independent-gain obstruction

For each observation time i, define raw w_i=(W_1(i),W_3(i)) and H(w)=1{|w_3|>|w_1|}. For any finite real a_i != 0,
H(a_i w_i)=1{|a_i||w_{3i}|>|a_i||w_{1i}|}=H(w_i).
This equality is **pointwise** for every finite q12 path and preserves H=0 on equality, including (0,0). Every a_i may vary arbitrarily across i and need not be stochastic, known, positive or independent. No new measurement apparatus is postulated.

This gauge is restricted to a *common* channel gain at each sample. It fails for separate gains a_{1i},a_{3i} (unless their absolute magnitudes coincide), offsets b_r, an erasing a_i=0, independent noise, nonlinear response or unknown assignment of the two Fourier features. Erasure sends all w to (0,0) and can change a strict event; it is not included. Finite multiplication must not overflow a computer representation if physically implemented, but this is a mathematical invariance only.

Let gamma0=1-Lambda3/Lambda6 and gamma+=1-(Lambda3+.01)/Lambda6 with 0<gamma+<gamma0<1. Under the **analytical stationary** Gaussian reference, W1~N(0,1), W3~N(0,v_gamma(alpha)), independent,
 v_gamma(alpha)=1+alpha(1/gamma-1).
If G1,G3 are independent standard Gaussians, W1=G1 and W3=sqrt(v)G3. Rotational symmetry in (G1,G3) gives
 p_gamma(alpha)=P(|sqrt(v)G3|>|G1|)=(2/pi) atan sqrt(v).
The angle of a standard planar Gaussian is uniform. Equality lines have zero Gaussian mass.

Because v_gamma+>v_gamma0, p_gamma+>p_gamma0 for every alpha>0. Differentiating at alpha=0:
 p'_gamma(0)=(1/gamma-1)/(2 pi),
 d_ang=p'_gamma+(0)-p'_gamma0(0)=[1/gamma+ -1/gamma0]/(2 pi)>0.
On a sufficiently small explicit existence interval (0,alpha_*] chosen by continuity of the positive derivative, p_gamma+(alpha)-p_gamma0(alpha)>=d_ang alpha/2. The interval is **existential**; no unjustified practical N0 is claimed. A useful pointwise derivative formula is
 p'_gamma(alpha)=(1/gamma-1)/[pi sqrt(v_gamma(alpha)) (1+v_gamma(alpha))].

Distinct channel gains erase identification in the single-time static Gaussian class: for fixed alpha, choose a1=1 and a3=sqrt(v_gamma0/v_gamma+) in the control model. Then its transformed W3 has the same variance as A7's W3 and transformed W1 has the same law. Thus the whole static Gaussian two-channel observation agrees. This is a nuisance overlap for *static limiting Gaussian laws*, not a full-history native impossibility theorem. A common gain cannot make the two events coincide, as the event is invariant and their probabilities are different.
