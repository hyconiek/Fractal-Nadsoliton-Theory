# NL-02 — intrinsic dual measure on the FIN S1 carrier

## Result

`PASS_INTRINSIC_1D_DUAL_MEASURE_AND_IRREGULAR_MESH_CONVERGENCE`.

For ordered intrinsic phases `theta_i` on the sourced circle, define circular
gaps `h_i = theta_(i+1)-theta_i`, Voronoi/control volume

`v_i=(h_(i-1)+h_i)/2`,

mass `m_i=mu v_i`, and edge conductance `c_i=kappa/h_i`.  No ambient Euclidean
coordinate, Delaunay triangulation or external density is needed in one dimension.

The resulting operator is

`(Lf)_i = (kappa/(mu v_i))*[(f_(i+1)-f_i)/h_i + (f_(i-1)-f_i)/h_(i-1)]`.

It obeys exactly:

- weighted conservation `sum_i m_i (Lf)_i = 0`;
- `M L` is symmetric;
- `-f^T M L f = kappa sum_i (f_(i+1)-f_i)^2/h_i >=0`.

Taylor expansion gives, for smooth f,

`Lf_i = (kappa/mu)[ f''(theta_i) + (h_i-h_(i-1)) f^(3)(theta_i)/3 + O(h_max^2) ]`.

Thus `h_max -> 0` is sufficient for consistency, without a lower bound on the
smallest gap.  Smooth distorted meshes show second-order convergence because
adjacent gap differences are O(h^2).  IID/Haar ordered samples are more irregular;
the deterministic replay still converges and the observed median max-error slope
is `-0.866223` over N=128..4096.

Numerical exact-identity checks:

- `||ML-(ML)^T||_inf = 5.684e-14`;
- Dirichlet-form defect = `1.137e-13`.

This directly strengthens LAW-02 for the actually established 1D FIN geometry.
It does **not** solve higher-dimensional dual-cell construction or source kappa/mu.
