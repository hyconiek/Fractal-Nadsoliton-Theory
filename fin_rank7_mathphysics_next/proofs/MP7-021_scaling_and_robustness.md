# MP7-021 — scaling equivalence versus genuine robustness

Scientific state: **PROVED_ANALYTIC_SCOPE**.

Let `A=X X^T` and

`V_{g,A}(p)=D(p||u0)-(g/2)(p-u0)^T A(p-u0)`.

For every `c>0`,

`V_{g,cA}=V_{cg,A}`.

Equivalently, rescaling `A` by `c` can be absorbed into the reciprocal choice of gain
normalization. In the dual factorization, `X -> sqrt(c) X` together with
`theta -> theta/sqrt(c)` preserves the label field `X theta`; this is a normalization
equivalence, not a derivation of either `c` or `g`.

For a nonuniform positive diagonal rescaling of retained feature columns, write
`X' = X S`. At fixed unscaled label field, the probability law is unchanged while

`M' = S^T M S`.

Ordinary covariance eigenvalues are not invariant under this congruence. The
coordinate-correct curvature object is the generalized pair `(M',G')`, with the
quadratic metric transformed together with the coordinates, as in MP7-019.

If a stationary solution is instead allowed to move under a perturbation parameter
`alpha`, write `F(theta,alpha)=0`. At a nondegenerate stationary point,

`d theta/d alpha = -H^{-1} partial_alpha F`.

Therefore a quantitative robustness radius requires both a perturbation bound and a
lower singular/eigenvalue bound for `H`. Near a fold no uniform ordinary implicit-
function bound survives. This separates normalization equivalence from genuine model
robustness and does not source a physical kernel or gain.

## Quantitative global neighborhood from MP7-020

MP7-020 now supplies the full-domain C4 bound

`lambda2(M) <= tau0-m`, with `m >= 6500000000000000000000000000/121203210760018863485060548125326941`.

For a fixed probability law and a diagonal feature rescaling `S`, the min-max
principle gives, for every ordered PSD eigenvalue,

`lambda_k(S M S) <= ||S||_2^2 lambda_k(M)`.

Thus Target P is preserved whenever

`||S||_2^2 (tau0-m) <= tau0`.

Equivalently, if the positive spectral weights change by
`lambda'_k/lambda_k <= rho` for all retained columns, it is sufficient that

`rho <= tau0/(tau0-m) = 8120615120921263853499056724396905047/8120613495921263853499056724396905047`,

so `rho-1 ~= 2.0010803381e-07`. In column-amplitude scaling this is
`||S||_2-1 <= 1.00054011964e-07`.

This radius is deliberately conservative: it uses only the single worst global
margin and the operator-norm comparison. It is a certified neighborhood, not an
estimate of the actual failure point. When stationary points move with the spectral
weights, a separate implicit-function bound is still required locally.

## Quantitative moving-root sensitivity at coexistence

For the unscaled aligned field coordinates `J`, write

`F(J,d)=J-g D m(J)=0`,

where `D=diag(d3,d4,d5,d6)` with `d3=lambda3/6`, `d4=lambda4/6`,
`d5=lambda5/6`, `d6=lambda6/12`.  If `alpha_j=log d_j`, then at a stationary
root

`partial F_i/partial alpha_j = -delta_ij J_i`.

Thus, with `A_J=partial_J F`,

`dJ/d alpha_j = A_J^{-1} e_j J_j`.

Direct interval inversion on the certified coexistence root box yields the full
matrix recorded in `results/MP7-021_local_spectral_sensitivity.json`.  In
relative form its infinity-row-sum bound is

`||d log J||_inf <= 2.622056124 ||d log lambda||_inf`

for infinitesimal perturbations at this root.  All displayed matrix entries are
positive on the certified box.

This is a local derivative condition number, not by itself a finite-size
parameter rectangle.  A finite perturbation theorem would additionally keep
`A_J` nonsingular throughout the entire connecting parameter box.
