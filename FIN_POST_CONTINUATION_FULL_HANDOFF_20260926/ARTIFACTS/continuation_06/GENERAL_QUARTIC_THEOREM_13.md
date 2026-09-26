# GENERAL-QUARTIC-THEOREM-13 — arbitrary retained quartic closure defect

Status: **PROVED WITHIN DECLARED FINITE-N HEAT-BATH + ME7 CLOSURE**.

This theorem is conditional on the same declared heat-bath generator and the
same local 7D maximum-entropy closure as KURTOSIS-THEOREM-11.  It is not a
FIN-derived physical law or measured observable.

## Statement

Let `Q=12`, `u=(1/12,...,1/12)`.  Let `V` be the retained real Fourier space
`k=3,4,5,6`, let `H` be the discarded `k=1,2` space, and denote the orthogonal
projectors by `P_V,P_H`.  Let `A7` be the strict retained self-adjoint operator.
For any real retained probe `phi in V`, define

    h = P_H(phi^2),
    k = P_H(phi * (A7 phi)).

For

    F_N(p) = (sqrt(N) phi^T p)^4

and `N` divisible by 12,

    lim N (L_full^3 - L_ME7^3) F_N(u)
      = -12 <h,h>_u + 2 g <h,k>_u.                 (1)

The `g^2` coefficient is identically zero.

For an eigenprobe `A7 phi=lambda phi`, (1) becomes

    (2 g lambda - 12) ||P_H(phi^2)||_u^2.

The pure strict `k=5` cosine/sine gives

    lambda5^2 (g lambda5 - 6) / 144,

so KURTOSIS-THEOREM-11 is a corollary.

## Generator expansion

Write `epsilon=1/N`, `p=u+epsilon d`, `v=P_V d`, `z=P_H d`.
The ME closure has

    p_ME = u + epsilon v + epsilon^2 w + O(epsilon^3),
    w = 6 P_H(v^2).

The heat-bath target has

    q = u + epsilon (g/12) A7 d
          + epsilon^2 (g^2/24)[(A7 d)^2-<(A7 d)^2>_u 1]
          + O(epsilon^3).

Expand the rate product into `G0 + epsilon G1 + epsilon^2 G2`.
Since `F_N=epsilon^2 f(d)`, `f=(phi.d)^4`, the coefficient of `epsilon`
in the third generator derivative is the sum of the six ordered compositions
whose h-orders add to two:

    G2 G0 G0, G0 G2 G0, G0 G0 G2,
    G1 G1 G0, G1 G0 G1, G0 G1 G1.

At `d=0`:
- `G2 G0 G0` vanishes because `G0^2 f` is constant;
- `G1 G1 G0` and `G1 G0 G1` vanish because the outer `G1` has zero
  coefficient at `d=0`;
- `G0 G2 G0` vanishes by the V/H orthogonality and the constant diagonal of
  the translation-invariant projector `P_V`;
- only `G0 G0 G2` and `G0 G1 G1` remain.

The independent 210-monomial producer confirms exactly this decomposition:
`002=-48||h||^2`, `011=+36||h||^2+2g<h,k>`, all other sequences zero to
roundoff, and the entire `g^2` tensor exactly zero in the producer arithmetic.

## The `G0 G0 G2` term

The full-minus-ME second-order departure term is

    Delta G2 = sum_ij [ z_i q1_j - u w_i ] Delta_ij.

The part linear in `g` has no degree-four component: its only possible degree-4
factor is proportional to

    sum_ij z_i (A7 d)_j (phi_j-phi_i)^2,

which is zero because both coefficient vectors have zero sum and `z.phi=0`.
Thus only the ME curvature `w=6 P_H(v^2)` contributes.

For `a=phi.d`, `h=P_H(phi^2)`, the degree-four part is

    -36 a^2 sum_n h_n v_n^2.

For one uniform jump `S=e_J-e_I`,

    E[(phi.S)(P_V S)_n] = (2/12) phi_n.

Using two independent jumps gives

    G0^2 [ a^2 sum_n h_n v_n^2 ](0)
      = (16/12) ||h||_u^2.

Therefore

    (G0 G0 Delta G2 f)(0) = -48 ||h||_u^2.       (2)

## The `G0 G1 G1` term: departure part

Factor `G1` into a departure operator and the target operator.  For a projector
`P` on mean-zero coordinates define

    B_P F = 12 E[(P d)_I (F(d+S)-F(d))].

The full departure uses the whole tangent projector `T`; ME7 uses `P_V`.
For mean-zero linear forms `L_x=x.d`, the degree-preserving identities needed
below are

    [B_P(L_x L_y)]_2
      = -L_x L_{P y} - L_y L_{P x},

    [B_P(L_x L_y L_z)]_2
      = L_x L_{P(yz)} + L_y L_{P(xz)} + L_z L_{P(xy)}.

A direct expansion of `B_P^2 a^4` gives its homogeneous quadratic part

    24 a L_{P(phi^3)}
    +72 <phi^2>_u a^2
    +12 a L_{P(phi P(phi^2))}
    +6 L_{P(phi^2)}^2.

Subtracting `P=P_V` from `P=T` and applying the outer `G0`, all V/H cross
terms vanish.  The two surviving contributions are `24||h||_u^2` and
`12||h||_u^2`, hence

    G0[(B_T^2-B_V^2)f](0) = 36 ||h||_u^2.         (3)

## The `G0 G1 G1` term: target part

Factor one power of `g` and write

    T_A F = E[(A7 d)_J (F(d+S)-F(d))].

Let `D=B_H=B_T-B_V`, `psi=A7 phi`.  The coefficient linear in `g` is

    G0 (D T_A + T_A D) f (0).

The required quadratic/cubic product identities are

    [T_A(L_x L_y)]_2
      = (1/12)[L_x L_{A7 y}+L_y L_{A7 x}],

    [T_A(L_x L_y L_z)]_2
      = (1/12)[L_x L_{A7(yz)} + cyclic].

After substitution, every quadratic term is a V/H inner product and vanishes
except

    (12/12) a L_{A7(phi h)}.

The outer `G0` gives

    (24/12) <h, P_H(phi psi)>_u
      = 2 <h,k>_u,                                (4)

where self-adjointness of `A7` was used:
`phi.A7(phi h)=(A7 phi).(phi h)`.

Equations (2)--(4) prove (1).

## Complete Z12 selection rule

Write the real retained probe using complex Fourier coefficients
`b3,b4,b5` and real `b6`.  The hidden Fourier coefficients of `phi^2` are

    H1 = 2( b4 conj(b3) + b5 conj(b4) + b6 conj(b5) ),

    H2 = 2 b5 conj(b3) + 2 b6 conj(b4) + conj(b5)^2.

Thus hidden `k=1` can be produced only by sector pairs

    (3,4), (4,5), (5,6),

and hidden `k=2` only by

    (3,5), (4,6), (5,5).

The corresponding weighted coefficients in `phi(A7 phi)` are

    K1 = (lambda3+lambda4)b4 conj(b3)
       + (lambda4+lambda5)b5 conj(b4)
       + (lambda5+lambda6)b6 conj(b5),

    K2 = (lambda3+lambda5)b5 conj(b3)
       + (lambda4+lambda6)b6 conj(b4)
       + lambda5 conj(b5)^2.

Parseval rewrites the theorem as

    C(phi) = -24(|H1|^2+|H2|^2)
             +4 g Re[conj(H1)K1+conj(H2)K2].

Hence a probe is a `g`-independent null control iff `H1=H2=0`.
Pure `k=3,4,6` satisfy this; pure `k=5` does not.

## Independent tensor certificate

`quartic_tensor_adjoint.cpp` computes the six ordered compositions once and
contracts their signed three-jump measures against all 210 homogeneous degree-4
monomials of the seven retained coordinates.  No probe is fitted.

`quartic_tensor_compare.py` independently constructs the RHS tensor from
`P_H(phi^2)` and `P_H(phi A7 phi)`.

Strict-kernel residuals:
- constant tensor: max absolute residual `2.50e-14`;
- coefficient linear in `g`: `5.11e-15`;
- `g^2` tensor: identically zero in the producer output.

The sequence-resolved producer gives:
- `002`: `-48 ||h||_u^2`;
- `011`: `+36 ||h||_u^2 + 2 g <h,k>_u`;
- `200,020,110,101`: zero to numerical roundoff.

## Methodological boundary

Promoted: exact algebraic theorem inside the declared finite-N heat-bath + ME7
model, with an independent full-tensor numerical certificate.

Not promoted: universal stochastic closure law, physical time, laboratory
observable, sourced dynamics, strict FIN physical prediction, SM/GR/ToE result.
