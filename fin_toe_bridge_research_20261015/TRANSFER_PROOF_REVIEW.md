# Native temporal reconstruction and critical transfer convergence

Date: 2026-10-15.
Status: **CONDITIONAL_ANALYTICAL_PROOF_SUBMISSION_FOR_AUDIT**.
This document derives a mathematical consequence under the explicitly listed
inherited premises. It does not amend the historical admission ledger, claim
machine-checked or peer-reviewed proof, or establish physical quantum mechanics.

## 1. Scope and the genuinely different object

The object here is the already specified interacting, reversible FIN count
process and its stationary path law. It is not a newly supplied free field,
temporal lattice, mass term, wave equation, or microscopic dynamics.

Use the connected local count domain `D_N` and auxiliary reflected generator
`L_D,N` admitted in the MT review. Internal rates and the unit-attempt-rate
clock are unchanged; only attempted exits are suppressed in this **analysis
process**. Its invariant probability is the actual count Gibbs law conditioned
on `D_N`, including multinomial multiplicity. Write this law as `pi_N`.
Reflection is not asserted to be a physical boundary or controller.

Fix one lambda for which the inherited local theory applies, and use the exact
critical gain convention `g_N=g_c+lambda*N^(-1/2)`, with N in the admitted
multiple-of-12 sequence. If the required predecessor process theorem is only
available at lambda=0, all operator-limit conclusions below are restricted to
lambda=0 until that premise is paid. No uniformity over lambda is implied.

The source of A7, the heat-bath kinetic rule, critical tuning, initial law and
attempt clock remain supplied mathematical premises. The new construction
does not derive any of them from nadsoliton ontology.

## 2. Explicit inherited premises to bind before admission

The proof requires the following precise facts, not a report's PASS label.

**P1 — Reversible parent.** The finite count chain on `D_N` is connected and
has generator `L_D,N`, row-sum convention, positive invariant weights `pi_N`,
and detailed balance. The admitted MT construction provides the relevant
candidate. Count and labelled generators must not be interchanged silently.

**P2 — Static critical law.** With `z_N(n)=N^(1/4)*m(n)`, the pushforward of
`pi_N` converges weakly to

```text
rho_lambda(dz) = Z_lambda^(-1) exp[-V_lambda(z)] dz,
V_lambda(z) = z^4/12 - a*z^2/2,
a = chi*lambda.
```

**P3 — Stationary critical two-time law.** For every fixed `t>=0` and bounded
continuous f,g, the **stationary reflected chain**, accelerated by `sqrt(N)`,
satisfies

```text
E_pi_N[g(z_N(X_0))*f(z_N(X_(sqrt(N)*t)))]
    -> E_rho_lambda[g(Z_0)*f(Z_t)],
```

where the limiting stationary diffusion is

```text
dZ_t = (a*Z_t-Z_t^3/3) dt + sqrt(2) dB_t.
```

For complex tests use conjugation on g. Joint stationary finite-dimensional
convergence is needed for the later all-cylinder claim, but the operator proof
itself uses only the displayed two-time convergence, including time `2t`.

P3 must be traced through the actual UF/MT/OBS proofs: a deterministic-start
UF theorem alone is not P3. Its extension to the local equilibrium initial
law, inherited fast moments, localization and reflected/full-process boundary
payment must really be available. No new numerical trajectories pay this
logical dependency. A proof only for a fixed compact-stopped process does not
automatically identify the unbounded stationary diffusion used here.

**P4 — Actual-parent path transfer, used only in section 8.** At an allowed
deterministic late origin `t_N`, the unchanged iid-start parent's observed
path on each fixed scaled horizon is TV-close to the stationary local
comparison, with initial and exit errors retained. P4 is not necessary for
the exact reflected operator theorem, and does not make the iid process
stationary.

These facts are to be bound to the real MT/UF/OBS reviews and completions,
not to unavailable original upload wrappers. In particular consult:

- `fin_critical_metastability_review_20261015/INTAKE_20261015.md` and
  `PROOF_COMPLETIONS.md`;
- `fin_uniform_critical_review_20261014/INTAKE_20261014.md` and the actual
  critical-process proof it admits;
- `fin_intrinsic_observation_review_20261015/INTAKE_20261015.md` and
  `PROOF_COMPLETIONS_AND_CORRECTIONS.md`.

## 3. Exact native temporal reflection positivity

Let X be a two-sided stationary version of the finite reversible chain.
Let theta reflect time, `(theta X)_t=X_(-t)`. For a bounded complex-valued
functional F of the nonnegative-time path, set `h(X_0)=E[F|X_0]`.

The past and future are conditionally independent given the **full count
state** X_0. Detailed balance identifies the conditional reflected past law
with the conditional future law. Consequently

```text
E[conj(F(theta X))*F(X)]
    = E[conj(E[F(theta X)|X_0])*E[F(X)|X_0]]
    = E[|h(X_0)|^2] >= 0.
```

The same argument gives the sesquilinear reflection form
`(F,G)_OS=<E[F|X_0],E[G|X_0]>_L2(pi_N)`. Quotienting its null space and
completing therefore recovers `L2(pi_N)`, provided time-zero functions are
included. If the initial cylinder algebra uses strictly positive times only,
finite-chain stochastic continuity recovers time-zero functions in its
closure. This is an explicit finite-dimensional temporal reconstruction,
not an appeal to an unverified full relativistic OS theorem.

The result holds for every bounded future cylinder and its L2 closure. It
does not require the observed soft coordinate alone to be Markov. Restricting
to functions of that coordinate inherits positivity, but need not recover
the entire finite count-state Hilbert space.

This is distinct from Program24's explicitly added free massive temporal
completion in `FIN_Programs_21_30_Executed_Monograph.md`: here positivity is
proved from the admitted interacting native path measure. No free covariance
or extra time coordinate is postulated.

## 4. Transfer generator and exact edge factorization

On `H_N=L2(pi_N)`, define

```text
A_N = -sqrt(N)*L_D,N,
T_N(t) = exp(-t*A_N).
```

Detailed balance makes A_N self-adjoint and nonnegative. Its closed quadratic
form is exactly

```text
E_N(u,u) = sqrt(N)/2 * sum_(x,y)
           pi_N(x)*L_D,N(x,y)*|u(y)-u(x)|^2.
```

There is no unspecified mobility or noise normalization in this formula.
The `sqrt(N)` factor is the already admitted critical time rescaling.

Writing `D_pi=diag(pi_N)`, the unitary map to the counting-measure space gives

```text
H_N = -sqrt(N)*D_pi^(1/2)*L_D,N*D_pi^(-1/2),
H_N*sqrt(pi_N) = 0.
```

For unordered configuration edges x,y and
`c_xy=pi_N(x)*L_D,N(x,y)=c_yx`, this matrix is the sum of positive rank-one
terms

```text
sqrt(N)*c_xy * (|x>/sqrt(pi_x)-|y>/sqrt(pi_y))
              (<x|/sqrt(pi_x)-<y|/sqrt(pi_y)).
```

It is a stoquastic transfer representation on **configuration space**. Its
edges are not physical spatial adjacency, and its spectrum is measured in
the supplied/rescaled model clock. The word Hamiltonian here names the
mathematical transfer operator, not a derived conserved physical energy.

## 5. Direct three-term proof of strong semigroup convergence

Take `H=L2(rho_lambda)` and the dense linear test space `C=C_b(R)` (complex
or real consistently). Define the linear maps

```text
J_N f(n) = f(z_N(n)).
```

P2 gives `||J_N f||_N -> ||f||_H`, including all inner products by
polarization. These are the asymptotically isometric core maps required
for varying Hilbert spaces. They need not be injective at finite N, nor
does their range need to equal the full finite count-state space.

The confining scalar SDE has a conservative Feller semigroup T(t), and its
invariant reversible law is rho_lambda. Here the Feller property needed is
`T(t)C_b subset C_b`, not strong continuity in the supremum norm on all C_b.
Local smoothness and `z*b(z)=a*z^2-z^4/3` provide the usual nonexplosion and
continuous-dependence localization proof. On H, T(t) is a strongly continuous
self-adjoint contraction semigroup.

For `f in C_b` and a fixed t, put `q=T(t)f in C_b`. Expand the actual native
norm, without assuming a Markov generator for z_N:

```text
||T_N(t)J_N f-J_N q||_N^2
 = <J_N f,T_N(2t)J_N f>_N
   -2 Re <J_N q,T_N(t)J_N f>_N
   +||J_N q||_N^2.
```

Self-adjointness and the semigroup law give the first term. P3 at `2t`
and at t, and P2 for the last term, give the limits

```text
<f,T(2t)f> = ||T(t)f||^2,
<q,T(t)f> = ||T(t)f||^2,
||q||^2 = ||T(t)f||^2.
```

Therefore

```text
||T_N(t)J_N f-J_N T(t)f||_N -> 0.
```

This is stronger than merely converging matrix elements. It works because
reversibility identifies the **norm at t with a correlation at 2t**. It does
not invoke, or establish, a uniform fast-fibre spectral gap.

To extend from core sections, use the standard varying-L2 definition: a
sequence `u_N` strongly approaches `u in H` if it can be approximated in the
iterated-limit native norm by `J_N f_j`, with `f_j -> u` in H. Native and
limiting semigroup contractions, together with the core result, transfer
strong convergence to every such section. This pays the required operator
notion; it is not fixed-space operator-norm convergence.

## 6. Resolvents, Mosco scope and what does not follow

For any alpha>0,

```text
(alpha+A_N)^(-1) = integral_0^infinity exp(-alpha*t)*T_N(t) dt.
```

For a core f, the corresponding difference from
`J_N (alpha+A)^(-1)f` tends to zero by the preceding norm convergence,
Minkowski's inequality and domination by
`2*||f||_infinity*exp(-alpha*t)`. The limiting resolvent preserves C_b by
the Feller property and dominated convergence. Uniform resolvent bounds
extend the statement to arbitrary strongly convergent sections.

Kuwae-Shioya Theorem2.4 now identifies the resulting strong resolvent/
semigroup convergence with Mosco convergence of the closed forms, **after
the map/space contract above and P1-P3 are verified**. This is a legitimate
way to pay Mosco without separately repeating an all-sequence liminf proof.
The theorem is not being applied to the non-Markov projection z_N itself:
its native input is A_N on the full count space.

Do not infer compact convergence, convergence of every low eigenvalue,
`gap_D ~ constant*N^(-1/2)`, or convergence of multiplicities. Those require
additional asymptotic compactness/compact spectral convergence. A simple
abstract warning is `H_N=R^2`, `H=R`, `J_N x=(x,0)`, `A_N=diag(0,1)` and
`A=0`. The lifted strong semigroups converge exactly, while the extra
eigenvalue1 persists. Bounded-energy vectors in the extra direction violate
asymptotic compactness. The inherited ten stable linearized directions do
not automatically remove this separate operator-level obligation.

If desired, bounded-observable spectral measures can be investigated with
their own tightness proof. Analytic continuation of numerical correlations
is neither required here nor a stable inverse problem paid by this theorem.

## 7. Exact critical differential operator and its domain

Let A be the nonnegative operator of the closure of

```text
E(f,f) = integral |f'(z)|^2 rho_lambda(dz),  f in C_c^infinity(R).
```

Its Markov realization is the conservative diffusion generator `A=-L`,
where on this core

```text
L = d^2/dz^2 + (a*z-z^3/3)*d/dz.
```

Use the unitary `Wf=sqrt(rho_lambda)*f` from H to `L2(dz)`. The transformed
closed form is the closure of

```text
Q(psi,psi) = integral |(d/dz+V_lambda'(z)/2)*psi(z)|^2 dz,
psi in C_c^infinity(R).
```

Define the transfer operator by this closed form, equivalently its
Friedrichs realization. Its differential expression is

```text
H_lambda = -d^2/dz^2 + (V_lambda')^2/4 - V_lambda''/2
 = -d^2/dz^2 + z^6/36 - a*z^4/6
   +(a^2/4-1/2)*z^2+a/2.
```

The statement specifies a domain by form closure; it does not equate an
arbitrary formal differential expression with a self-adjoint operator.
Smooth positive multiplication maps the compactly supported smooth cores
onto each other, paying the unitary equivalence of the form closures.

The normalized ground state is `psi_0=Z_lambda^(-1/2)*exp(-V_lambda/2)`.
Cutoffs approximate the constant function in the weighted gradient form,
so psi_0 belongs to the transformed form domain and has form energy zero.
It consequently belongs to the operator kernel. Factorization proves
nonnegativity even though the displayed potential is negative in places.
At lambda0 the expression is

```text
H_0 = -d^2/dz^2+z^6/36-z^2/2,
psi_0 proportional to exp(-z^4/24).
```

This is a derived one-coordinate transfer theory. No physical Planck
constant, particle mass, physical position, spatial momentum or Born
measurement postulate has been obtained by writing it this way. Distinct
operators sharing the admitted quartic critical limit share this limiting
object; its existence does not select nature's A7.

## 8. What the unchanged late iid-start parent receives

P4 transfers probabilities and expectations of **bounded finite-window path
functionals**, with the actual initial and exit error. For a two-sided window
centered at t_N, its earlier endpoint must still lie after the paid local
relaxation time and its later endpoint before the admitted trapping horizon.
No conditioning on successful survival is silently inserted.

For example, if its path law is within delta_N of the stationary comparison
in the convention `TV=sup_event |P-Q|`, then for `|F|<=C` the real part of
its reflection form is at least `-2*C^2*delta_N`. Thus asymptotic RP for such
records can be inherited with a paid scope. This is not exact stationarity
or exact finite-N reflection positivity for the iid-start process.

Neither a bounded path-TV approximation nor the limit operator theorem
automatically transfers unbounded moments, relative entropy, an operator
norm, the whole global Gibbs spectrum, or the actual parent's global
stationary dynamics. H_N above belongs to the explicitly named reflected
stationary comparison. Its physical use must not erase that distinction.

## 9. Audit checklist, failure gates and research value

The submitted algebraic implications are complete under P1-P3. Before they
are promoted to an accepted FIN theorem, the execution/audit agent must:

1. Bind the actual local domain, count law, generator and critical clock,
   and verify that predecessor P3 really covers stationary reflected
   two-time laws at every fixed t and2t. If missing, export precisely that
   lemma, not a new fitted diffusion or an unqualified operator PASS.
2. Check the varying-Hilbert definitions and the three-term norm identity
   independently, including `J_N T(t)f`, complex conjugation and dense-core
   extension. Matrix dimensions alone are not this check.
3. Bind the limiting conservative diffusion/form realization and Friedrichs
   domain. Do not insert a kinetic prefactor or change normalization during
   the ground-state transformation.
4. Keep stationary operator statements separate from bounded actual-parent
   predictions. A useful finite-N error rate, practical burn-in, and a
   numerical N0 are not part of this conditional asymptotic proof.
5. Preserve the physical promotion gate: reconstructing a transfer Hilbert
   space does not source actual complex preparations, noncommuting physical
   instruments, composition, interference or a Born readout law.

The proposed result changes a real mathematical dependency: an admitted
native statistical mechanism would supply an entire temporal positive
hierarchy and its controlled critical transfer structure, without another
microscopic law. Its honest terminal status can be
`DERIVED_CLASSICAL_TEMPORAL_TRANSFER_THEORY`. That is more informative than
an empty-input audit, but not evidence that FIN has derived QM or ToE.

## 10. Primary reading and theorem-use limits

- [Henley, From classical to quantum dynamics at Rokhsar-Kivelson points](https://arxiv.org/pdf/cond-mat/0311345),
  sections3.4 and4, equations3.7 and4.1-4.5: detailed-balance symmetrization
  and the distinction between imaginary-time correlation correspondence
  and actual real-time probability evolution. No physical RK model is
  imported into FIN by using this algebra.
- [Kuwae and Shioya, Convergence of spectral structures](https://www.intlpress.com/api/bgcloud-front/resource/pdf/volume/1805786134762446850-1805786134762446850-e776617e92e2a310c0026f43fb47be64.pdf),
  Assumption2.1, Definition2.4, Definitions2.11-2.13 and Theorem2.4; contrast
  Theorem2.6. The varying-space maps, closed forms and self-adjoint
  contraction semigroups are the explicit prerequisites used above.
- [Grothaus and Wittmann, Mosco convergence of gradient forms with non-convex interaction potential](https://arxiv.org/pdf/2105.05140),
  version2024, section1.1 and Theorem3.4: an additional account of
  varying-L2/resolvent/Mosco equivalences and domain obligations. Its
  application-specific gradient-density criteria are not assumed to hold
  automatically for the FIN jump chain.
- [Jaffe and Ritter, Quantum field theory on curved backgrounds via the Euclidean functional integral](https://www.arthurjaffe.com/Assets/documents/CurvedSpace.pdf),
  Definition1.5, sections2.1,2.5 and2.6, Theorems2.7-2.9: useful benchmark
  separating positivity, quotient Hilbert space, heat semigroup and formal
  unitary continuation. Their Riemannian/static-background hypotheses are
  not present in FIN and their spatial QFT theorem is not applied here.

All cited sections were inspected in captured primary PDF text. The direct
finite-state conditioning and three-term arguments above are supplied so
that resemblance to the literature is not itself counted as proof. No
external experiment, spending, data acquisition, physical source promotion,
selector/legacy-role transfer, QM/GR/SM or ToE closure is licensed.
