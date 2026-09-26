# FIN — RELATIONAL CLOCK / ESCAPE-PROFILE CONTINUATION

Date: 2026-09-26
Repository baseline inspected: hyconiek/Fractal-Nadsoliton-Theory, current visible HEAD fe14a6f4e436815635df54429102f23c22862296.

Status:
- analytic theorem inside the declared strict finite FIN operator model;
- independent regeneration/strengthening of the OP-002 route;
- numerical finite-window checks for the frozen strict kernel;
- NOT a derivation of SI time, a physical Born rule, or a selected microscopic dynamics.

## 1. Object

Let A be a real symmetric Laplacian-like operator with A 1 = 0. Fix a source i.
For j != i define

    a_j = A_{ji},    b_j = (A^2)_{ji},    c_j = (A^3)_{ji},
    V = sum_{j!=i} a_j^2.

Given an off-diagonal probability/readout record p_{ji}, define

    S = sum_{j!=i} p_{ji},
    r_j = p_{ji}/S.

Thus S says "how much escaped" and r says "where it went", conditioned on escape.
The unparameterized curve (S,r) is unchanged by any monotone reparameterization
of the external control variable.

## 2. Unitary-like channel

For U(t)=exp(-i t A),

    U_{ji} = -i a_j t - (b_j/2)t^2 + i(c_j/6)t^3 + O(t^4),

hence

    |U_{ji}|^2
      = a_j^2 t^2
        + [b_j^2/4 - a_j c_j/3] t^4
        + O(t^6).

Let

    u4_j = b_j^2/4 - a_j c_j/3,
    Q_U = sum_j u4_j.

Then

    S_U = V t^2 + Q_U t^4 + O(t^6),

and

    r_U = r0 + K_U S + O(S^2),

where

    r0_j = a_j^2/V,

    K_U,j =
      (1/V) [u4_j/V - a_j^2 Q_U/V^2].

Therefore, for any norm with K_U != 0,

    D_U(S) := ||r_U-r0|| = kappa_U S + O(S^2),

so the clock-free relational exponent is

    gamma_U = d log D_U / d log S -> 1.

## 3. Zero-velocity wave-like channel

For C(t)=cos(t sqrt(A)),

    C_{ji} = -(a_j/2)t^2 + (b_j/24)t^4 + O(t^6),

so

    |C_{ji}|^2
      = (a_j^2/4)t^4 - (a_j b_j/24)t^6 + O(t^8).

Let

    T = sum_j a_j b_j.

Then

    S_W = (V/4)t^4 - (T/24)t^6 + O(t^8),

and

    r_W = r0 + B sqrt(S) + O(S),

with

    B_j =
      a_j^2/(3 V^(3/2)) * [T/V - b_j/a_j]

for nonzero a_j.

Hence, when B != 0,

    D_W(S) := ||r_W-r0|| = kappa_W sqrt(S) + O(S),

and therefore

    gamma_W -> 1/2.

This separates the two declared dynamic categories without knowing the
control-to-time map.

## 4. Exact exceptional class

Assume a connected real symmetric circulant Laplacian and all off-diagonal
a_j are nonzero.

B=0 for every destination iff

    (A^2)_{ji}/A_{ji} = constant

for every j != i.

Circulant symmetry makes the diagonal residual constant, hence

    A^2 = c A + d I.

Applying both sides to 1 and using A 1 = 0 gives d=0, so

    A^2 = c A.

Thus every eigenvalue satisfies lambda(lambda-c)=0. Connectedness leaves
one zero eigenvalue and a single repeated positive eigenvalue. Therefore

    A = c (I - J/N),

which is precisely the equal-weight complete-graph Laplacian up to scale.

So the wave sqrt(S) coefficient vanishes identically only in this exceptional
one-positive-eigenvalue class. The strict FIN operator is not in that class.

## 5. Strict FIN numerical constants

Using the frozen strict kernel

    W_ij = cos(0.18575 d + 0.1625)/(1+d^1.8),  i != j,
    d = min(|i-j|,12-|i-j|),

and A=diag(W1)-W, the spectrum is

    0,
    0.754121154207079 (x2),
    1.577049514427609 (x2),
    1.961406861976445 (x2),
    2.199568849333209 (x2),
    2.298606272079096 (x2),
    2.342182041146300.

For one source node:

    V = 0.5379635797962747
    ||B||_2 = 0.08542176796653574
    ||K_U||_2 = 0.10683540605435281
    range_j [(A^2)_{ji}/A_{ji}] = 8.388027687996473

so the exceptional condition is very far from being satisfied numerically.

The leading heat/diffusion conditional profile is proportional to W_ji,
whereas the unitary/wave common leading profile is proportional to W_ji^2.
Their strict-kernel L2 separation at S->0 is

    ||r_heat,0 - r0||_2 = 0.20412204975578369.

Thus the three declared categories have:
- heat: different limiting destination profile already at order S^0;
- unitary: same squared-weight r0 but D ~ S;
- wave: same squared-weight r0 but D ~ sqrt(S).

## 6. Finite-window check

Exact matrix-function evaluation for the frozen strict operator gives
representative local slopes gamma=d log D/d log S:

    S=1e-6:   gamma_U = 1.00000055, gamma_W = 0.50046437
    S=1e-5:   gamma_U = 1.00000546, gamma_W = 0.50147268
    S=1e-4:   gamma_U = 1.00005456, gamma_W = 0.50470005
    S=4e-4:   gamma_U = 1.00021830, gamma_W = 0.50952929
    S=1e-3:   gamma_U = 1.00054606, gamma_W = 0.51531234
    S=1e-2:   gamma_U = 1.00550632, gamma_W = 0.55361578

A 1% asymptotic-slope tolerance survives approximately to:
- unitary: S ~ 1.8e-2,
- wave:    S ~ 4.3e-4.

These are numerical finite-window diagnostics, not apparatus power guarantees.

At equal escape level, the full conditional profiles are also distinct. For
example the L2 separation behaves as

    ||r_W(S)-r_U(S)|| / sqrt(S) -> 0.0854218...

and numerically is still about 0.08364 at S=0.01.

## 7. Exact scale-gauge no-go

For every c>0,

    exp[-i t (cA)] = exp[-i (ct) A],
    exp[-t (cA)]   = exp[-(ct) A],
    cos[t sqrt(cA)] = cos[(sqrt(c)t) sqrt(A)].

Therefore the UNPARAMETERIZED curve (S,r) is exactly invariant under a global
operator rescaling A -> cA for all three channel classes (with the appropriate
time reparameterization).

Consequences:
1. the escape-profile curve can classify relational dynamic category;
2. it can support a dimensionless internal progress coordinate once a category
   and normalized A are declared;
3. it cannot determine the absolute physical rate, seconds, hbar-like action
   scale, or a unique dimensional normalization.

This is the precise mathematical sense in which a time-like relational
structure can emerge before a physical clock scale.

## 8. Interpretation for FIN

This result strengthens a specific emergence picture:

    relation -> transformation -> changed relation

can carry an observable temporal signature before an external time coordinate
is introduced. The order/type of change can be encoded in relations among
observables themselves.

But the result does NOT show that FIN has derived physical time. The exact
A -> cA gauge proves that an absolute clock scale remains additional
information. It also does not select unitary, wave or heat dynamics from
statics alone; it only shows how to distinguish those declared categories
from records if the readout assumptions hold.

## 9. Next physical research atom

The next useful test is RELATIONAL-CLOCK-MEMORY-02:

Combine the clock-free escape-profile exponent with the two established memory
mechanisms:
1. dynamic tree storage (linear poles/residues),
2. finite-N hidden heat-bath memory.

Question:
Does gamma=1 versus gamma=1/2 survive after eliminating latent memory, or can
a legal hidden-memory channel change the leading exponent and alias the
dynamic category?

Acceptance:
- a theorem giving conditions under which gamma is invariant under hidden
  memory / calibrated memoryless observation; or
- an explicit counterexample showing that memory can change the exponent.

Only after that should the result be treated as a robust candidate physical
diagnostic.
