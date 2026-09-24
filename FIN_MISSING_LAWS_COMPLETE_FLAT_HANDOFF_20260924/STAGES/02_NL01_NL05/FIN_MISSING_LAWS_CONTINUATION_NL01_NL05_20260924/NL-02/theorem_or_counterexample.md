# Intrinsic 1D finite-volume theorem

Let points be cyclically ordered on S1 with positive gaps h_i and h=max h_i.
Use nodal control volumes v_i=(h_{i-1}+h_i)/2 and conductances c_i=kappa/h_i.
Then the weighted graph operator is conservative and self-adjoint in the
M-inner product. For f in C^4,

L f(theta_i) = (kappa/mu) [ f''(theta_i)
  + (h_i-h_{i-1}) f^(3)(theta_i)/3 + O(h^2) ].

Hence every sequence with h->0 is pointwise consistent. Shape regularity is
not required for consistency, though it affects rates/conditioning.
