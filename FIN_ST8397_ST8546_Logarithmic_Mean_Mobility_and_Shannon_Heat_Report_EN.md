# FIN ST8397–ST8546 — Logarithmic-Mean Mobility and Shannon Heat

## Executive theorem

Let \(p\) be a strictly positive probability vector on the twelve FIN vertices
and let

\[
H(p)=\sum_i p_i\log p_i.
\]

Define the logarithmic mean

\[
\Lambda(a,b)=\frac{a-b}{\log a-\log b},
\qquad
\Lambda(a,a)=a,

\]

and the state-dependent Onsager operator whose edge conductances are

\[
c_{ij}(p)=W_{ij}\Lambda(p_i,p_j).

\]

Then \(K_p\) is symmetric, positive semidefinite and mass-conserving, and

\[
\boxed{K_p\log p=Ap.}

\]

Consequently the Shannon gradient flow is exactly the strict heat equation:

\[
\boxed{
\dot p=-K_p\nabla H(p)=-Ap.}

\]

A nonuniform twelve-state replay gives residual

\[
1.1102230246251565\times10^{-16}.

\]

## Edge-local uniqueness

Consider symmetric edge-local mobilities

\[
c_{ij}=W_{ij}M(p_i,p_j).

\]

Exact agreement with heat flow requires, on every edge,

\[
M(a,b)(\log a-\log b)=a-b.

\]

Therefore

\[
\boxed{M(a,b)=\Lambda(a,b)}

\]

for \(a\neq b\), and continuity fixes \(M(a,a)=a\). The logarithmic mean is
thus unique in the declared Shannon, symmetric, edge-local, heat-matching
class.

## Lyapunov structure

Along the flow,

\[
\frac{dH}{dt}
=-(\log p)^TAp
=-\frac12\sum_{i,j}W_{ij}(p_i-p_j)
(\log p_i-\log p_j)\le0.

\]

The uniform distribution is the unique equilibrium on the connected strict
graph. The heat semigroup preserves mass, positivity and the simplex.

## Falsification: entropy–metric duality

For a general strictly convex entropy

\[
F_\phi(p)=\sum_i\phi(p_i),

\]

define

\[
M_\phi(a,b)=
\frac{a-b}{\phi'(a)-\phi'(b)}.

\]

Then the corresponding Onsager operator also satisfies

\[
K_{\phi,p}\nabla F_\phi(p)=Ap.

\]

Thus heat flow alone does not select an entropy–metric pair. Shannon selects
the logarithmic mean only after Shannon entropy is independently privileged.
The FIN informational interpretation motivates that choice but is not yet a
uniqueness theorem over all convex information functionals.

## Dual dynamics boundary

The result closes the finite heat shadow only. The unitary channel

\[
U_t=e^{-itA}

\]

is Hamiltonian/isometric rather than an entropy gradient flow. A wave channel
requires a phase-space or symplectic extension. Shared \(A\) organizes these
channels but does not turn them into one dynamics.

## Refinement

For any supplied positive product refinement

\[
A_{24}=A_{12}\otimes I_2+I_{12}\otimes B_q,

\]

the same logarithmic-mean construction gives

\[
K^{(24)}_p\log p=A_{24}p.

\]

It therefore transports functorially after \(A_{24}\) is chosen. It does not
select the fibre rate \(q\), its level sequence or the refinement clock.

## Configuration-space limitation

The theorem acts on the twelve-vertex probability simplex. It does not derive
the 4096-state configuration generator, ternary Gibbs coupling, circulation or
hidden-state memory. Those objects live in a different algebra and their
moment hierarchy does not close to the vertex heat equation.

## Physical boundary

Multiplying the mobility by \(c>0\) produces \(cA\) and restores the clock
torsor. Shannon entropy is dimensionless. Physical heat, temperature and time
still require a bath, \(k_B\), energy calibration and timestamps.

## Gate

| Requirement | Result |
|---|---|
| logarithmic-mean construction | PASS |
| exact Shannon heat-gradient identity | PASS |
| edge-local uniqueness | PASS |
| entropy Lyapunov theorem | PASS |
| refinement functoriality for supplied \(A_{24}\) | PASS |
| strict Shannon uniqueness among all entropies | FAIL |
| strict refinement-rate law | FAIL |
| configuration-kinetics closure | FAIL |
| physical clock and temperature | FAIL |
| OA evidence | FAIL |

\[
\boxed{5/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

## Deepest interpretation

FIN gains a canonical time-symmetric frenetic geometry for its finite heat
channel once Shannon entropy is fixed. This is the first explicit mobility
section derived from \(W\), state \(p\) and an information functional. It is
not a complete physical section because entropy choice, clock, refinement,
other dynamical channels and apparatus remain open.

The highest-value next theorem is to test whether composability, locality and
information additivity uniquely select Shannon entropy inside FIN. Without
that theorem the result must remain Shannon-conditional.

No physical temperature, clock, continuum, laboratory evidence, Standard
Model, gravity, \(L_{\rm total}\), or Theory of Everything closure follows.
