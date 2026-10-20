# RP004-B — logarithmic-window strong remainder and endpoint trace estimate

Date: 2026-10-19. Disposition: **LOG_WINDOW_COMPARISON_PAID_CONDITIONALLY_ON_FROZEN_SMOOTH_CAUCHY_MODEL**. This is an analytic PDE estimate, not numerical validation, scattering or positive transported flux.

## 1. Fixed coupled system, data, symmetry and coefficients

In the inherited invariant real e1 subspace only (other complex species remain part of the full action but are exactly zero along this frozen datum), use

  phi_tt-phi_xx+2phi(phi²-1)+6phi psi²=0,
  psi_tt-psi_xx+(a+6phi²)psi+g4 psi³+g6 psi⁵=0,
  a∈[-3/4,-1/4], |g4|<=M4, |g6|<=M6,
  phi(0)=K=tanh x, phi_t(0)=0,
  psi(0)=epsilon f, psi_t(0)=0,
  f=sqrt(3/4) sech² x.

Here M4,M6 are declared finite coefficient bounds on the **unchanged supplied coefficients**; they are not fitted, nor an extra potential. By parity and uniqueness, phi-K is odd, psi is even, all other species stay zero and the kink translation mode is not excited for these prepared data. Do not generalize this special parity reduction to arbitrary 25-real-field data.

Set L=-∂²+4-6 sech² x and La=L+a+2. Let omega=sqrt(a+2), v=f cos(omega t), and let u be the zero-data retarded solution from RP004-A. Let eta=phi-K=epsilon²u+w, psi=epsilon v+r; the exact zero-data remainders satisfy

 w_tt+Lw=-6K(psi²-epsilon²v²)-6eta psi²-6K eta²-2eta³,
 r_tt+La r=-12K eta psi-6eta²psi-g4 psi³-g6 psi⁵.  (B.1)

Neither RHS was appended as an independent radiation law.

## 2. Positive graph energies and their uniform wave estimates

L on the odd subspace has point-spectrum floor 3 (odd shape mode) and continuum floor 4; the even translation zero mode is absent by *solution symmetry*. La on the even subspace has floor a+2>=5/4; the first eigenfunction f is retained. These are the inherited exact Pöschl–Teller spectral facts, and every bound is uniform in a on the frozen closed interval.

For A=L_odd or La_even, define

    G_A[y,y_t]²=||A y||²_2+||A^(1/2)y_t||²_2.

Since A=-∂²+V(x), V and its derivatives are uniformly bounded, and A>=cI>0 on the chosen parity subspace, elliptic graph-norm estimates and the positive quadratic form imply

    c_* (||y||_H2+||y_t||_H1) <= G_A[y,y_t]
                <= C_* (||y||_H2+||y_t||_H1).  (B.2)

A scalar spectral oscillator or Duhamel calculation establishes for y_tt+Ay=F with F∈L¹_t H¹_x and parity preserved,

    ||(y(t),y_t(t))||_(H2×H1)
      <= C [||(y(0),y_t(0))||_(H2×H1)+∫_0^t ||F(s)||_H1 ds].  (B.3)

Indeed one measures A y and A^(1/2)y_t; spectral sine/cosine operators have norm <=1, and ||A^(1/2)F||_2<=C||F||_H1. This route avoids commutator loss and uses only the fixed smooth potential. In 1D, H¹ embeds into L∞, and H¹ is an algebra. Thus H² positions/H¹ velocities control pointwise first-position-derivative and velocity traces on every endpoint uniformly in its location.

The reference equations and B.3 give uniformly:

    sup_t ||(v(t),v_t(t))||_(H2×H1)<=C_v,
    ||(u(t),u_t(t))||_(H2×H1)<=C_u (1+t).      (B.4)

The latter is linear in t because forcing S(1+cosΩt) is uniformly H¹-bounded; it contains DC/turn-on and the odd shape mode.

## 3. Uniform nonlinear bootstrap

Let tau=1+t, X(t)=||(r,r_t)||_(H2×H1), Y(t)=||(w,w_t)||_(H2×H1). Fix a<=-1/4 within its exact interval; no parameter sweep. Constants below depend on M4,M6 and fixed spectral/Sobolev bounds, not epsilon or t.

Assume X<=A0 epsilon³ tau², Y<=B0 epsilon⁴ tau³ over a maximal bootstrap interval and epsilon² tau³<=c0 sufficiently small. Then B.4 and 1D products imply

  ||psi||_H2<=C epsilon, ||eta||_H2<=C epsilon² tau,
  ||F_r||_H1<=C [epsilon³ tau + epsilon Y + epsilon² tau X]
               <= C epsilon³ tau[1+(A0+B0)epsilon² tau²],
  ||F_w||_H1<=C [epsilon X+X²+epsilon⁴ tau²+epsilon² Y]
               <= C epsilon⁴ tau²[1+(A0+B0)epsilon² tau²],       (B.5)

where omitted higher powers are absorbed by the stated smallness, not discarded or sign-selected. Fixed finite signed g4,g6 enter only the C constants. Apply B.3 to the zero-data remainders and integrate tau and tau². Choose A0,B0 >=4C with the relevant integration factors and epsilon0 so that (A0+B0) epsilon² (1+T)^3<=c1. The improved strict bounds close the continuity bootstrap:

    X(t)<=C_r epsilon³(1+t)²,
    Y(t)<=C_w epsilon⁴(1+t)³     for 0<=t<=T.    (B.6)

An explicit admissibility statement is: there are fixed K0>0 and epsilon0>0, depending only on the frozen family bounds, such that (B.6) holds whenever epsilon<=epsilon0 and epsilon²(1+T)^3<=K0. Since epsilon²[1+12 log(1/epsilon)]³->0, the **original unmodified window**

    T_epsilon=12 log(1/epsilon)

is admissible for all sufficiently small epsilon. No global-in-time H² estimate or nonlinear scattering is claimed.

## 4. Boundary current error on the exact frozen spatial ray

By the 1D trace embeddings, the same bounds hold for |r_t|,|r_x| <=CX and |w_t|,|w_x|<=CY at x=±R, with constants independent of R. The reference u has |u_t|+|u_x|<=C tau. The bound eigenfunction v has |v_t(t,±R)|+|v_x(t,±R)|<=C exp(-2R), while K'(±R)<=C exp(-2R).

Expand the **original** current J=-phi_t phi_x-psi_t psi_x; no subtraction of the adverse original terms:

 J_exact = -epsilon² K'u_t -epsilon² v_t v_x -epsilon⁴ u_tu_x
           -w_t(K'+epsilon²u_x)-epsilon²u_t w_x-w_t w_x
           -epsilon v_t r_x-epsilon r_t v_x-r_t r_x.          (B.7)

Consequently, denoting the signed two-endpoint integral F_exact and F_u=∫[ -u_tu_x(R)+u_tu_x(-R) ]dt, one gets

 |F_exact(R,T)-epsilon⁴ F_u(R,T)|
    <= C [epsilon² e^(-2R)(1+T)²
          +epsilon² e^(-4R)(1+T)
          +epsilon⁴ e^(-2R)(1+T)^4
          +epsilon⁶ (1+T)^5
          +epsilon⁸ (1+T)^7].                                  (B.8)

At *exactly* the unmodified R_epsilon=2 log(1/epsilon), T_epsilon=12 log(1/epsilon), the RHS is bounded by

    C [epsilon⁶ (1+T)^5 +epsilon⁸ (1+T)^7+epsilon¹⁰ (1+T)],

and is o(epsilon⁴ T). (Terms epsilon⁸ tau⁴, epsilon⁶ tau² have been dominated by epsilon⁶ tau⁵.) This is the required trace-admissible **same-trajectory** comparison at the frozen log window.

In particular if an independently *proved* RP004-A theorem supplied F_u(T/6,T)>=c_ref T for all sufficiently large T with c_ref>0 uniform in a, B.8 would imply F_exact(R_epsilon,T_epsilon)>0 for all sufficiently small epsilon. This is a **conditional implication**, not a statement that the missing moving-ray reference theorem has been proved.

## 5. Audit and exclusions

Re-derive B.1 from the original potential V=(phi²-1)²/2+(a+6phi²)psi²/2+g4 psi⁴/4+g6 psi⁶/6. The independent symbolic verifier checks potential derivatives, both remainder sign conventions and parity. It confirms 6K f²=2S. That is a bounded algebraic control, not a proof of the infinite-dimensional energy estimate B.3.

The parity gap, graph-energy bound, 1D products and bootstrap are analytic arguments with explicit dependence on the frozen coefficients; no unspecified dispersive decay has been assumed. The full 25-real-field orbital stability remains the inherited NS conditional result for the parent trajectory. No assertion is made about the natural selection of coefficients or about a photon, physical units, measured flux, general kernel-sensitive radiation or ToE.
