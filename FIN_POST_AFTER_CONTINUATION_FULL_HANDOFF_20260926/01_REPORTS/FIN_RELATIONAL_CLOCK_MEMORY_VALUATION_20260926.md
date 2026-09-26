# FIN — RELATIONAL CLOCK MEMORY / VALUATION THEOREM

Date: 2026-09-26
Scope: mathematical continuation of the clock-free escape-profile analysis.
Repository context: current visible HEAD fe14a6f4e436815635df54429102f23c22862296.

Status:
- exact asymptotic theorem for analytic positive/subprobability records;
- exact robustness/alias criteria for memoryless linear observation and causal
  linear memory filters;
- conditional connection to FIN tree Schur memory;
- NOT a physical clock derivation and NOT a claim that tree memory is the
  physical detector/readout law.

## 1. General order-of-vanishing theorem

Let the off-diagonal/readout vector have an expansion

    p(t) = t^m [a + t^q b + O(t^(q+delta))],

where m>0, q>0, delta>0,
A = 1^T a > 0, and p is nonnegative on a sufficiently small one-sided window.

Define

    S(t) = 1^T p(t),
    r(t) = p(t)/S(t),
    r0 = a/A.

Write B = 1^T b and

    eta = b/A - a B/A^2
        = (I-r0 1^T)b/A.

Then

    S(t) = A t^m [1 + (B/A)t^q + O(t^(q+delta))],

and

    r(t) = r0 + eta t^q + O(t^(q+delta)).

Since

    t^q = A^(-q/m) S^(q/m) [1+o(1)],

if eta != 0,

    ||r-r0||
      = ||eta|| A^(-q/m) S^(q/m) [1+o(1)].

Therefore the unparameterized relational exponent

    gamma := lim_{S->0+} d log ||r-r0|| / d log S

is

    gamma = q/m.

This number depends only on the orders of vanishing, not on a monotone
reparameterization of the external control/time variable.

## 2. Interpretation of the three declared FIN channels

For the records already used in the FIN channel programme:

### Heat/diffusion

    p_j(t) = W_ji t + O(t^2).

Generically m=1, q=1, hence gamma=1, but its limiting profile is

    r0,H ∝ W_ji.

### Unitary probability

    p_j(t) = A_ji^2 t^2 + O(t^4).

Thus m=2, q=2 and

    gamma_U=1,

with

    r0,U ∝ A_ji^2.

### Zero-velocity wave squared response

    p_j(t) = (A_ji^2/4)t^4 + O(t^6).

Thus m=4, q=2 and

    gamma_W=1/2,

with the same leading profile as the unitary record.

Consequently a more informative clock-free fingerprint is the pair

    F = (r0, gamma),

not gamma alone.

For the frozen strict FIN kernel:
- heat has a different r0 from unitary/wave;
- unitary and wave share r0 but have gamma 1 versus 1/2.

This separates the three DECLARED channel classes in the nonexceptional
small-S regime. It is not a universal classifier of all dynamics.

## 3. Why gamma is not a 'quantum signature'

The theorem immediately shows that gamma is only a ratio of Taylor/valuation
orders.

Any unrelated channel with the same pair (m,q) has the same gamma.
For example a classical hidden process engineered so that escape starts at
order t^4 and the first transverse profile correction starts at order t^6
also has gamma=1/2.

Thus gamma=1/2 must never be promoted to 'wave', 'quantum', or Born-rule
evidence without the declared channel class and readout assumptions.

## 4. Perturbation theorem: which hidden memory preserves gamma?

Let

    p_eps(t) = t^m [a + sum_s t^s c_s + t^q b + ...].

For each coefficient c define its transverse/projective component

    Pi_a(c) = c/A - a (1^T c)/A^2.

Let

    s_* = min { s>0 : Pi_a(c_s) != 0 },

including the intrinsic term b at s=q.

Then, provided the leading vector a itself is unchanged,

    gamma = s_*/m.

Therefore:

1. A lower-order correction s<q that is exactly proportional to a is
   shape-preserving and cancels from r=p/S.
2. A lower-order correction with any transverse component changes gamma to s/m,
   no matter how small its coefficient is, if the asymptotic S->0 limit is taken.
3. A perturbation that changes the leading order m or leading ray a changes the
   entire fingerprint and must be treated as a different observation law.
4. Cancellation of the first transverse coefficient raises gamma to the next
   surviving order; this is an identifiability singularity, not evidence for a
   new dynamical category.

This is an asymptotic theorem. On finite windows a very small lower-order
contamination may become visible only below a crossover scale.

## 5. Memoryless detector / calibration map

Let a fixed linear observation map M act as

    y(t)=M p(t).

Then

    y(t)=t^m [Ma+t^q Mb+...].

If 1^T Ma>0 and the transformed transverse correction is nonzero,

    Pi_{Ma}(Mb) != 0,

the exponent is still

    gamma_y=q/m.

So a MEMORYLESS linear confusion/efficiency map does not change the valuation
orders. It can, however, destroy identifiability if it projects the first shape
correction onto the leading ray or annihilates it.

An invertible M does not automatically guarantee preservation of the normalized
transverse component under an arbitrary choice of output summation functional,
but generic nonsingular calibrated maps preserve the order unless this explicit
collinearity condition is met.

## 6. Causal linear memory filter

Consider an analytic causal observation law

    y(t) = D p(t) + integral_0^t K(t-s) p(s) ds,

with

    K(u)=K0 + K1 u + O(u^2).

### 6.1 Direct feedthrough D != 0

If Da is nonzero, the direct term keeps leading escape order m.

The memory term starts as

    K0 a * t^(m+1)/(m+1).

Therefore the first possible memory-induced profile correction is one order
after the leading ray.

If

    Pi_{Da}(K0 a) != 0,

then

    q_eff = min(q,1)=1

and

    gamma_eff = 1/m.

For the declared channels this means, GENERICALLY under such a stateful
observation filter:
- unitary: gamma 1 -> 1/2,
- wave:    gamma 1/2 -> 1/4.

If K0 a is projectively parallel to Da, the t^(m+1) memory term cancels from
the normalized profile. Then the next term must be inspected; the original
gamma may survive.

### 6.2 Pure-memory readout D=0

Now

    y(t) =
      K0 a t^(m+1)/(m+1)
      + K1 a t^(m+2)/[(m+1)(m+2)]
      + K0 b t^(m+q+1)/(m+q+1)
      + ...

Generically the leading order becomes m+1 and the first transverse correction
appears one order later, so

    gamma_eff = 1/(m+1).

For a wave-like base record (m=4,q=2), a generic pure-memory readout therefore
has gamma=1/5, not 1/2.

Again this is a statement about a declared filter law, not an assertion that a
FIN tree kernel is literally a detector probability filter.

## 7. Connection to the admitted FIN tree-memory object

The repository's proof-grade conditional tree result has

    Lambda(z)=L_BB-L_BI(L_II+z C_I)^(-1)L_IB,

and time-domain kernel

    K_H(t)=L_BI exp[-C_I^(-1)L_II t] C_I^(-1)L_IB.

Hence

    K_H(0)=L_BI C_I^(-1)L_IB

is generically nonzero.

If a future physical bridge makes this K_H enter an escape-profile observation
through a convolution of the form in section 6, then the decisive robustness
test is NOT merely 'are the residues small?'. It is

    Is Pi_leading(K_H(0) a) = 0 ?

If yes, the leading memory correction is shape-preserving and gamma can survive.
If no, memory changes the asymptotic exponent irrespective of how small the
residue amplitude is.

This provides a concrete bridge condition between the existing tree-memory
mathematics and the clock-free observable programme.

## 8. Finite-N hidden heat-bath memory

The admitted stationary heat-bath result says the leading Gaussian hidden
residual decouples after the accepted linear transform and genuine stationary
reduced feedback begins at O(1/N).

That is an expansion in 1/N, not automatically an expansion in short time t.
It therefore does NOT by itself prove preservation of gamma.

Away from stationarity the repository explicitly leaves N^(-1/2) hidden
dependence open. The missing NONSTATIONARY EDGEWORTH task is therefore now
directly relevant to the emergent-time programme:

    determine the first t-order at which nonstationary hidden feedback changes
    the escape ray or its first transverse correction.

This is a sharper physical motivation for that P0 than 'complete the
asymptotic calculation'.

## 9. New no-go / positive result

Positive:
- a dimensionless temporal ordering exponent gamma=q/m can be defined from
  relations among observables without an external clock parameter;
- it is invariant under monotone reparameterization and global rate scale;
- memoryless calibrated linear observation generically preserves its order.

No-go:
- arbitrary hidden/stateful memory can change gamma;
- arbitrarily weak lower-order transverse leakage wins asymptotically;
- therefore gamma is not an intrinsic property of the static FIN operator A;
  it belongs to the typed triple

      (dynamics, preparation, observation/memory law).

This sharply narrows what 'time emerges from FIN' can currently mean.

## 10. Recommended next atom

NONSTATIONARY-RELATIONAL-CLOCK-01

Use the declared finite-N heat-bath generator and the already accepted
stationary Edgeworth coordinates. For a specified nonstationary initial family:

1. expand the visible transition/readout record simultaneously in t and N^-1/2;
2. project every correction with Pi_a;
3. identify the first nonzero transverse bidegree (t^s, N^-k/2);
4. decide whether the clock-free gamma is stable at fixed large N, in the
   N->infinity-first limit, and in coupled t=t(N) windows.

Acceptance:
- a uniform wedge in (t,N) where gamma is unchanged, with explicit remainder;
  or
- an explicit lower-order hidden correction proving asymptotic aliasing.

That would connect the current relational-clock result to FIN's own finite-N
memory mechanism rather than to an externally declared filter.
