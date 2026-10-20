# Audit and triangular repair of the submitted strong remainder proof

Date: 2026-10-19. One directly motivated repair within the submitted graph-
energy method. Original B bytes/status are retained; the result is admitted
with this completion, not by treating its algebraic flags as a PDE proof.

## 1. What needed correction

The submitted (B.5) bounds the term `epsilon X` in F_w by an expression
whose dependence on the bootstrap coefficient A0 is multiplied by an extra
small factor. But `X<=A0 epsilon^3 tau^2` gives
`epsilon X<=A0 epsilon^4 tau^2` with **no** such factor. Choosing A0 and B0
both above a supposedly independent C does not by itself close those
displayed inequalities. This is not a counterexample to their target orders.

## 2. Pay the strong linear norm

For `A=L_odd` or `L_a,even`, inherited spectral gaps are 3 and at least 5/4.
The smooth potential is bounded uniformly in a. On these subspaces

    G_A^2=||A y||_2^2+||A^(1/2)y_t||_2^2

is uniformly equivalent to `||y||_H2^2+||y_t||_H1^2`: the gap first controls
L2, `y''=V y-Ay` controls H2, and the quadratic-form identity plus the gap
controls H1 velocities. The reverse bounds follow from bounded V. The
spectral sine/cosine group consequently gives, for zero data,

    ||(y,y_t)(t)||_(H2 x H1)<=C integral_0^t ||F(s)||_H1 ds.   (12)

Prove (12) first for smooth data or directly by spectral Duhamel, then
extend by density. There is no unbounded commutator or loss of a derivative.
In one dimension H1 is an algebra; H2 positions have H1 first derivatives.
Thus both velocities and position derivatives have bounded point traces,
uniform in the location of the trace.

The exact reference has `||(u,u_t)||<=C_u tau`, while `||(v,v_t)||<=C_v`,
where `tau=1+t`. These bounds do not assume the moving-cone theorem.

## 3. Fix constants in the correct order

Use the exact equations and zero remainder data from the shared contract:

    eta=epsilon^2 u+w, psi=epsilon v+r,
    F_r=-12K eta psi-6eta^2 psi-g4 psi^3-g6 psi^5,
    F_w=-6K(2epsilon v r+r^2)-6eta psi^2-6K eta^2-2eta^3.

Let `X=||(r,r_t)||`, `Y=||(w,w_t)||` in H2 x H1, and initially bootstrap
`X<=A epsilon^3 tau^2`, `Y<=B epsilon^4 tau^3`.
Require `A epsilon^2 tau^2<=1`, `B epsilon^2 tau^2<=1` and epsilon<=1.
Then, with constants **independent of A,B**,

    ||psi||_H2<=(C_v+1)epsilon,
    ||eta||_H2<=(C_u+1)epsilon^2 tau.                         (13)

The product estimates give `||F_r||_H1<=C epsilon^3 tau`, provided
`epsilon^2 tau<=1`. Applying (12) first improves

    X<=C_r epsilon^3 tau^2,                                (14)

with C_r fixed without A or B. Next the F_w equation gives

    ||F_w||_H1<=C[epsilon X+X^2+epsilon^4 tau^2].             (15)

Terms `eta^3` and `eta psi^2` are included in the last term using (13) and
the same smallness, not discarded. Substituting the **improved** (14), the
square term is `C_r^2 epsilon^6 tau^4`; its ratio to epsilon^4 tau^2 is
`C_r^2 epsilon^2 tau^2`. Under a fixed smallness bound this gives

    Y<=C_w epsilon^4 tau^3,                                (16)

where C_w can depend on the already fixed C_r. Now choose A=4C_r, B=4C_w,
and then K0>0 small enough that `epsilon^2(1+T)^3<=K0` pays all the
smallness conditions above. Equations (14),(16) strictly improve the
bootstrap and continuity closes it up to T. This triangular ordering avoids
the circular constant choice in the submission. Constants depend on the
fixed Sobolev/spectral/source bounds and finite M4,M6 only, not epsilon,t,a.

Local H2 x H1 wave theory follows from the same Duhamel estimate and local
polynomial Lipschitz products. Bounded norms exclude a finite continuation
failure before this T. For these smooth data higher regularity is propagated
on finite intervals. The inherited NS theorem supplies the independent
full-parent orbital statement; it is not used as a dispersive estimate.

Since `epsilon^2(1+12log(1/epsilon))^3 ->0`, the original logarithmic
window is admissible. No potential, preparation, radius or clock is changed.

## 4. Audit the complete signed current

At either fixed endpoint x=+/-R, exact expansion gives

    J=-epsilon^2 K' u_t-epsilon^2 v_t v_x-epsilon^4 u_t u_x
      -w_t(K'+epsilon^2 u_x)-epsilon^2 u_t w_x-w_t w_x
      -epsilon v_t r_x-epsilon r_t v_x-r_t r_x.               (17)

Here `J=-phi_t phi_x-psi_t psi_x`; no current is post-subtracted. K', v_x
and v_t at the endpoints are O(exp(-2R)), while u derivatives are O(tau).
Use (14),(16) and the trace embedding from section 2 on **each** endpoint,
then integrate. All terms, including the localized bound-mode and kink
interference, are paid by

    |F_exact-epsilon^4 F_u| <= C[
       epsilon^2 exp(-2R)tau_T^2 + epsilon^2 exp(-4R)tau_T
       +epsilon^4 exp(-2R)tau_T^4
       +epsilon^6 tau_T^5+epsilon^8 tau_T^7].                 (18)

For example `epsilon v_t r_x` has bound epsilon^4 exp(-2R)tau_T^3,
which is dominated by the displayed third term; `w_t w_x` gives the last
term. Thus (18) includes all signed cross terms even when their favorable
sign is unknown. At `R=2log(1/epsilon)`, `T=12log(1/epsilon)`, it is at most

    C[epsilon^6 tau_T^5+epsilon^8 tau_T^7+epsilon^10 tau_T]
      =o(epsilon^4 T).                                     (19)

This proves the submitted B target with a corrected constant argument.
No fixed-T remainder is silently used at a growing horizon, and no H1 x L2
estimate is mistaken for a derivative trace bound.
