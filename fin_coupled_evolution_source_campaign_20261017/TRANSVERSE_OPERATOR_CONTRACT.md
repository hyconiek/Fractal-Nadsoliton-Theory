# UE-002: transverse operator contract

Start with the **user-supplied** real transverse interface on the infinite line

`L_perp = -d_x^2 I + A + 2Y Phi_K(x)^2`, `A=M2+K_s`.

For the following arguments `A,Y` must be *real symmetric* constant n×n matrices, the kinetic metric Euclidean positive definite, and `Phi_K(x)=v tanh(kx)` with `v,k>0`. This full statement is a **testable conditional interface**, NOT an assertion that actual 12-complex-component FIN satisfies it. Set `B=2v² Y` and `V_inf=A+B`. Then

`L_perp = -d_x² I + V_inf - B sech²(kx)`.

The natural self-adjoint realization is on `H²(R;R^n)` in `L²(R;R^n)`; its closed bounded-below quadratic form is on `H¹`. Bounded symmetric multiplication is a relatively bounded perturbation of the free Laplacian, so self-adjointness is immediate. Since `B sech²(kx)` vanishes at spatial infinity, essential spectrum starts at `lambda_min(V_inf)` and contains the continuous half-line(s) determined by eigenvalues of `V_inf` (spectral union `[lambda_min(V_inf),infinity)` for the present 1D equal-gradient operator).

**Mandatory full-sector caveat:** the actual action reportedly has 12 *complex* `Psi` fields, i.e., 24 real transverse coordinates, not just the 12 real coordinates of the submitted operator. A 12×12 slice operator cannot certify full transverse stability unless conjugation invariance and block Hessian comparison have been paid. Mixed real-imaginary blocks or phase-sensitive terms can create additional directions.

No wall translation-mode claim in `Phi` follows from examining only `Psi`, and no conversion of any eigenfrequency to physical mass is licensed.
