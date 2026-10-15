# SRN-004 — Frozen actual-parent angular classifier: conditional consistency

Fix N in 12N, k_N=floor(N^(63/64)), alpha_N=k_N/N, M_N=floor(N^(1/4)), eta_N=N^(-1/32), h=1. The fixed block is predetermined, with H_i=1{|W3(t_N+i)|>|W1(t_N+i)|}; ties are zero. Let A_N=M_N^-1 sum_{i=0}^{M_N-1}H_i and theta_N=(p_gamma0(alpha_N)+p_gamma+(alpha_N))/2. Predict the comparator if A_N>theta_N and A7 otherwise; this tie convention is frozen and correct for p_control>p_A7. The sole classifier input is (H_0,...,H_{M-1}).

The ideal gap Delta_N=p_control(alpha_N)-p_A7(alpha_N)>= (d_ang/2) alpha_N for all sufficiently large N by the analytical angular calculation. Under each reflected stationary model, a wrong or midpoint-tie decision implies |A_N-p_model(alpha_N)|>=Delta_N/2 >=(d_ang/4)alpha_N. Apply SRN-003 with that fixed c, including all ties. Using the exact integer floors (asymptotic equivalents only eventually), each error tends to zero by:

 eta_N^-2/(M_N alpha_N²)=O(N^-5/32),
 eta_N^-3 N^-1/4/alpha_N=O(N^-9/64),
 eta_N^-3 k_N^-1/2/alpha_N=O(N^-49/128),
 eta_N/alpha_N=O(N^-1/64),
 k_N^-1/2/alpha_N=O(N^-61/128),
 eta_N^-3 (N/sqrt(k_N)) pi_D(B_D)/alpha_N = polynomial(N)*pi_D(B_D)=o(1)
by the inherited MT positive departure energy and polynomial horizon. The slowest explicit polynomial error is N^-1/64. No useful numeric N0 or native constant is established.

At the common deterministic post-MT origin t_N, couple each original iid-uniform unreflected labelled parent to the corresponding stationary reflected labelled chain until the departure time. The coupling error is at most initial outside probability + preparation TV + exit probability through the *last recorded time*, bounded by the previous full M_N horizon. It is o(1) by the inherited MT assumptions for **each** operator. Transfer the single final bounded *wrong-classification event* by adding this coupling probability, once. Do not transfer scaled unbounded moments by TV or postselect surviving runs. Thus each original parent's model-conditional misclassification probability ->0.

Raw H is exactly invariant to independently time-varying common nonzero sample gains a_i. This conclusion does *not* transfer to independent channel gain, offsets, measurement noise, erasures, sign-only temporal records, minimax/necessary resource bounds or real detectors. It is conditional on independent admission of the inherited MT comparator local windows and four-coordinate native L1 proof; the finite diagnostics do not mechanically prove those analytical inputs.

The eventually-valid condition `k_N<=N/2` requires roughly `N>=2^64` (up to the exact floor and the N in 12N restriction). This window is mathematically asymptotic and not an experimentally feasible sample size.
