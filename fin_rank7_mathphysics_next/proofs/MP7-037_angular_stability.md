# MP7-037 — positive odd block at interior aligned stationary points

Assume all three paired mode amplitudes `a3,a4,a5` are strictly positive.  Let `phi=(phi3,phi4,phi5)` be their phase coordinates, with the alternating amplitude fixed (either sign is reducible by translation).

From MP7-007, after sign correction if needed,

`Z(phi)=sum_m c_m e^(i m.phi)`

with `c_m>=0`, `c_{-m}=c_m`, and exact support `L6` for positive alternating amplitude or `L12` when the alternating amplitude vanishes.  Absolute convergence with polynomial frequency weights follows from the factorial Bessel-series bound, so the series may be differentiated twice termwise.

At `phi=0`,

`grad_phi Z=0`,

and

`- Hess_phi log Z = [sum_m c_m m m^T]/Z(0) =: Q`.

The quadratic dual cost is phase-independent, hence

`Hess_phi Phi = Q`.

For any nonzero `x in R^3`,

`x^T Q x = [sum_m c_m (m.x)^2]/Z(0)`.

The support lattices `L6` and `L12` have rank three (explicit bases were proved in MP7-009), and every support coefficient is strictly positive in the interior. Thus some supported `m` has `m.x !=0`, and `Q` is **strictly positive definite**.

## Cartesian odd block

Let `(x_k,y_k)` be the Cartesian cosine/sine amplitude pair and use the phase convention

`x_k=s_k cos(phi_k)`, `y_k=-s_k sin(phi_k)`.

At `phi=0`, `dx_k/dphi_k=0`, `dy_k/dphi_k=-s_k`, and `d^2x_k/dphi_k^2=-s_k`. Therefore the chain rule gives

`Q = D H_odd D - diag(s_k partial_{s_k} Phi)`,

where `D=diag(s3,s4,s5)` and `H_odd` is the Cartesian sine/sine Hessian block.  This identity is a nonlinear-coordinate Hessian formula; away from stationarity the gradient term is essential.

At a full stationary point `grad Phi=0`.  Since all three `s_k>0`,

`H_odd = D^(-1) Q D^(-1) > 0`.

Hence every aligned stationary point with all three paired amplitudes nonzero has a strictly positive three-dimensional odd/sine Hessian block.

This does not cover boundary points with a zero paired amplitude, and it does not assert that arbitrary nonaligned stationary points can be rotated to aligned form.
