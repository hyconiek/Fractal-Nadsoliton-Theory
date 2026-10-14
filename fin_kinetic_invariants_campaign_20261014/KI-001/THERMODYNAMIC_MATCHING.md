# KI-001: primal/dual fold matching (mathematical proof)

## Premises and chart
The already supplied finite-copy model has q=12, the exact parameter-defined real matrix A=XX^T, Gibbs free energy F_g(p)=sum_i p_i log(12p_i)-(g/2)p^T A p on the open simplex, and theta-potential Phi_g(theta)=||theta||^2/(2g)-log[(1/12)sum_i exp(X_i theta)]. The fold (p_f,theta_f,g_f,v) is an isolated nondegenerate *source-certified* root, where p_f=softmax(X theta_f), theta_f=g_f X^T p_f, v^T v=1, and (I/g_f-X^T Sigma_f X)v=0. Positive remaining H7 eigenvalues and p_f>0 are inherited conditions, not established afresh by this campaign.

Define T={u:1^Tu=0}, Sigma_f=diag(p_f)-p_f p_f^T, r=Sigma_f Xv. Then 1^T r=0 and g_f X^T r=v. Thus the normalized thermal left functional ell_HB=g_f Xv satisfies ell_HB^T r=1. Note ell_HB is a covector in occupation space, not a new variable or physical clock.

## Hessian congruence and inertia
For tangent z, H_F=diag(1/p_f)-g_f XX^T represents the second variation of F, with constant covectors irrelevant. Introduce y=diag(p_f)^(-1/2) z, with y orthogonal to sqrt(p_f). Then

z^T H_F z = y^T[I-g_f B B^T]y,
B=(I-sqrt(p_f)sqrt(p_f)^T) diag(sqrt(p_f)) X.

The nonzero eigenvalues of B B^T are the nonzero eigenvalues of B^T B=X^T Sigma_f X=M_f. The seven eigenvalues of I-g_f M_f are g_f times those of H7, one zero and six positive. Four additional tangent dimensions carry eigenvalue 1 in the congruent form. Thus the constrained free energy has exactly one null and ten positive directions; this is a *congruence/inertia statement*, not equality of all Hessian eigenvalues in Euclidean coordinates.

Indeed H_F r = [diag(1/p_f)Sigma_f Xv]-g_f X X^T Sigma_f Xv = Xv-[p_f^T Xv]1-Xv = -[p_f^T Xv]1. Hence T^T H_F r=0. The full occupation null direction is r.

## Coefficients in precisely the same normalization
For g-derivative of the constrained gradient, partial_g grad F=-A p, so

 a_F=r^T partial_g grad F(p_f,g_f) = -r^T A p_f = -(v^T theta_f)/g_f^2 = a.

The energy term is quadratic and the entropy third derivative is -sum_i z_i^3/p_i^2. Writing y_i=(Xv)_i, m=sum p_i y_i, r_i=p_i(y_i-m),

 b_F=D^3 F(p_f,g_f)[r,r,r] = -sum_i p_i(y_i-m)^3 = b.

These are exact identities, with a,b inherited from the established dual fold orientation; no recalibration or numerical fitting is introduced.

## Transverse elimination does not alter cubic term
For local coordinates p=p_f+xr+w, w in any complementary stable tangent subspace, the constrained stationary w solving the transverse gradient equation obeys w=O(x^2+|g-g_f|) by the implicit-function theorem (invertible fast Hessian). The quadratic term H_F[r,w] vanishes, so the x^3 and (g-g_f)x coefficients in the reduced energy remain b/6 and a respectively. Terms of degree four and beyond can depend on the elimination chart; this is not a claim that the full reduced energy is exactly cubic.

## Numerical but nonproof checks
`kinetic_analysis.py` and `independent_verify.py` calculate this equality in binary64, with ~1e-16 differences. Their absolute numeric centers are not an exact-fold enclosure. Exact-zero and inertia statements above rely on the source-certified fold, including its original interval-backend conditions.
