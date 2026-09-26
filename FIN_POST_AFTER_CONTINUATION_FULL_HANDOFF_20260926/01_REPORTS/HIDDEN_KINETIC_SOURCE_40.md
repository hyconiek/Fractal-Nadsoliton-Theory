# HIDDEN-KINETIC-SOURCE-40
## OU inheritance, exact kinetic nonuniqueness, and the moving-phase closure obstruction

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact conditional projection theorem from the admitted heat-bath OU residual;
- exact same-statics/different-kinetics no-go for the new scalar u;
- exact moving-phase closure obstruction;
- no physical clock, SI rate, or unique kinetic law is derived.

This continues `RECIPROCAL-HIDDEN-GEOMETRY-ACTION-39`.

---

## 1. Existing ingredients

The continuation has produced the conditional chain

    H4 hidden field z
      -> localized intrinsic phase chi
      -> scalar u
      -> reciprocal geometry interaction E_int.

The admitted finite-N heat-bath lane separately gives, at leading Gaussian
order,

    dz = -z dt + dW_z

in the declared dimensionless heat-bath clock.

The repository also already proves that the same static energy can support
different kinetic categories (`PHY-021`) and that eliminating local storage
degrees of freedom generates positive-pole memory.

The question is whether any of these selects the kinetic law of u.

---

## 2. Frozen-phase OU inheritance

Write the two complex hidden Fourier amplitudes as

    z1  (k=1),
    z2  (k=2),

and the localized retained phase as

    chi = exp(i theta).

For fixed theta define the two reflection-even scalar channels

    u1 = Re[z1 exp(-i theta)],
    u2 = Re[z2 exp(-2 i theta)].

Under the isotropic OU residual

    dz_k = -z_k dt + dW_k,

with theta fixed,

    du1 = -u1 dt + dW_u1,
    du2 = -u2 dt + dW_u2.

Therefore any fixed linear combination

    u = alpha u1 + beta u2

also obeys a first-order OU/relaxational law

    du = -u dt + dW_u

with projected noise covariance determined by alpha,beta and the hidden
covariance.

Thus the finite-N heat-bath object supplies a **conditional first-order kinetic
law for u** without adding a new functional form.

But its time unit and even the choice of heat-bath microscopic update remain
declared premises, not FIN-derived physics.

---

## 3. Point-evaluation scalar

Under the zero-order locality selector of report 38, the canonical scalar is,
up to normalization,

    u = u1 + u2.

For a frozen intrinsic phase it therefore inherits the same OU decay rate 1 in
the declared heat-bath clock.

At uniform equilibrium, where the four real hidden components have covariance
I4/12, the point-evaluation scalar has equilibrium variance

    Var(u) = 1/6

under the orthonormal Fourier convention used here.

A global rescaling of the microscopic attempt rate changes the decay rate and
noise strength together while preserving the static equilibrium law, exactly as
in `DYN-002`.

---

## 4. Moving phase: exact connection terms

Now allow

    chi(t)=exp(i theta(t)).

Define comoving complex amplitudes

    w_k = z_k exp(-i k theta)
        = u_k + i v_k.

If the deterministic hidden OU part is

    dot z_k = -z_k,

then

    dot w_k
      = -w_k - i k dot(theta) w_k.

Hence exactly

    dot u_k = -u_k + k dot(theta) v_k,
    dot v_k = -v_k - k dot(theta) u_k.          (1)

Noise rotates covariantly in the same moving frame.

For k=1,2:

    dot u1 = -u1 + omega v1,
    dot v1 = -v1 - omega u1,

    dot u2 = -u2 + 2 omega v2,
    dot v2 = -v2 - 2 omega u2,

where

    omega = dot(theta).

Thus the scalar

    u=u1+u2

obeys

    dot u
      = -u + omega(v1+2v2) + noise.             (2)

Equation (2) is **not closed in u alone** whenever the intrinsic phase moves.

---

## 5. Exact closure obstruction

Take two hidden states with identical scalar u but different quadratures
(v1,v2).

At the same theta and the same nonzero omega, equation (2) gives different
instantaneous dot u unless

    v1+2v2

also agrees.

Therefore no autonomous first-order law

    dot u = F(u,theta,omega)

can reproduce all admissible hidden states.

At minimum an additional state variable carrying the weighted quadrature is
required, or one must integrate it out and accept memory.

This obstruction is kinematic: it follows from representation transport in a
moving relational frame, not from a chosen potential.

---

## 6. Constant angular velocity spectrum

For constant omega, each k-sector has deterministic eigenvalues

    -1 +/- i k omega.

So a purely relaxational OU mode in the fixed Fourier frame becomes a damped
rotating pair in the phase-locked frame.

The decay envelope remains exp(-t), but the relationally measured scalar
contains oscillatory components at omega and 2 omega.

This is an explicit mechanism by which a moving order parameter changes the
observed temporal category without changing the underlying OU generator.

---

## 7. Adiabatic elimination

If |omega| << 1 and the quadratures relax rapidly, setting

    dot v_k approximately 0

gives

    v_k approximately -k omega u_k

and therefore

    dot u_k
      approximately -(1+k^2 omega^2) u_k.

So the leading correction is second order in angular velocity; there is no
linear-in-omega correction to the decay rate after adiabatic elimination.

However, when omega changes on the relaxation timescale, this local
approximation fails and the eliminated quadratures generate history dependence.

---

## 8. Same static energy admits at least three inequivalent kinetic classes

Let the reciprocal interaction/self-energy define, after linearization, a
positive scalar stiffness K for u.

The same static energy

    E(u)=K u^2/2

is compatible with all of the following admitted/conditional structures.

### A. Local relaxational / heat-bath

    eta dot u + K u = f.

Transfer function:

    chi_R(s)=1/(K+eta s).

This is the category inherited from the declared OU/gradient lane.

### B. Inertial / wave-like

    M ddot u + K u = f.

Transfer function:

    chi_I(s)=1/(K+M s^2).

This is the category already used conditionally in the refinement-compatible
wave construction and in `PHY-021`.

### C. Passive finite-pole memory

Dynamic elimination of positive internal storage gives a stiffness of the
Stieltjes form

    Lambda(s)
      = K + sum_r b_r s/(s+lambda_r),

with

    b_r >= 0,
    lambda_r > 0,

possibly plus a local eta s term.

Thus

    chi_M(s)
      = 1/[K+eta s+sum_r b_r s/(s+lambda_r)].

This is the scalar version of the admitted tree-storage memory mechanism.

All three have the same static susceptibility

    chi(0)=1/K.

They differ dynamically.

---

## 9. Exact nonuniqueness theorem

Static Dirichlet/refinement data plus the reciprocal interaction determine K
and the coupling shape, but they do not determine whether the dynamic
denominator contains

    s,
    s^2,
    or positive-pole rational terms.

Therefore current FIN data do not select a unique temporal category for u.

This is stronger than saying that one coefficient is unknown: the **order and
state-space dimension of the evolution law remain nonunique**.

The admitted objects instantiate at least two of these categories directly:
- heat-bath OU -> local first-order relaxation;
- tree storage -> finite-pole memory.

The repository's conditional inertial refinement model provides a third.

---

## 10. Distinguishing observable

A small impulse/step experiment separates the classes without changing the
static energy.

For an impulse:

- relaxational response has a finite jump followed by monotone exponential
  decay;
- inertial response has zero displacement jump but finite velocity jump and
  oscillatory/second-order behavior;
- finite-pole memory produces a sum of additional relaxation scales and
  generally a higher Hankel rank.

Equivalently the high-frequency asymptotics differ:

    chi_R(s) ~ 1/(eta s),
    chi_I(s) ~ 1/(M s^2),

while finite-pole memory contains identifiable rational corrections.

Thus kinetic category is operationally distinguishable but not statically
sourced.

---

## 11. Consequence for emergent time

There are now two separate ways memory enters:

1. **microscopic hidden storage**, already present in tree/finite-N reductions;
2. **frame-induced hidden state**, created when the intrinsic relational phase
   itself moves and the phase-locked scalar is projected from H4.

The second mechanism is especially important for the emergence picture:
even a Markovian hidden law can become non-Markovian after expressing it only
through a scalar tied to a changing relational frame.

So "memory" need not be fundamental storage; it can also be generated by
coarse-graining a moving relational representation.

---

## 12. Boundary

Not concluded:
- OU is the physical law of u;
- the dimensionless OU rate 1 is one inverse second;
- inertia is selected;
- tree memory is the physical bath;
- theta is physical spatial position or physical time;
- GR, SM, QW-2191 or role-bearing L_total are closed.

---

## 13. Next research atom

### MOVING-PHASE-MEMORY-KERNEL-41

Eliminate the quadratures v1,v2 exactly for a prescribed theta(t).

For constant omega, derive the exact scalar memory/transfer kernel for
u=u1+u2 and determine its minimal realization order.

For slowly varying theta(t), derive the first nonlocal correction and identify
the geometric connection/holonomy term.

Acceptance:
- exact closed integro-differential equation or minimal state-space realization;
- proof of whether k=1 and k=2 can be distinguished from the scalar record;
- explicit separation between genuine storage poles and frame-induced poles.

This is the next direct test of whether FIN's relational phase can generate a
structured notion of temporal memory without inserting a new bath.
