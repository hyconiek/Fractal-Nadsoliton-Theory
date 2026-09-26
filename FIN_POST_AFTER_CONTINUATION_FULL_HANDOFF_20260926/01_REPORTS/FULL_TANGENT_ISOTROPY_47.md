# FULL-TANGENT-ISOTROPY-47
## Heat-bath reset is unique under full centered-space isotropy, but FIN symmetry does not imply that premise

Date: 2026-09-26

Status:
- exact uniqueness theorem for reset dynamics;
- exact representation-theoretic no-go for deriving the premise from D12 alone;
- no physical clock or microscopic update is selected.

---

## 1. Uniqueness theorem

Let Q be a continuous-time Markov generator on a finite interior state p.

Assume:

1. Q 1 = 0;
2. p is stationary;
3. every p-centered function f obeys

       Q f = -gamma f

   for one gamma>0.

Every function decomposes uniquely as

    f=(p^T f)1 + f_0

with p^T f_0=0.

Therefore

    Qf = -gamma f_0
       = gamma[(p^T f)1-f].

Hence as a matrix on functions,

    boxed:
    Q = gamma(1 p^T-I).

So the complete-reset/heat-bath generator is unique up to the one global clock
rate gamma.

No detailed-balance assumption is even needed once full centered-space
isotropy is imposed.

---

## 2. Why D12 does not imply full isotropy

At the uniform state the eleven-dimensional centered label representation of
D12 decomposes into real Fourier sectors

    k=1,2,3,4,5,6.

The k=1,...,5 sectors are two-dimensional real irreducible planes and k=6 is
one-dimensional.

A D12-invariant self-adjoint generator may therefore assign independent rates

    gamma_1,...,gamma_6

to these inequivalent sectors.

D12 only forces cosine/sine degeneracy inside each k-plane.

Thus symmetry permits six kinetic rates, not one.

The heat-bath choice

    gamma_1=...=gamma_6

is a special additional isotropy.

---

## 3. Retained/hidden split

FIN's current strict split is

    hidden:   k=1,2
    retained: k=3,4,5,6.

So even at the maximally symmetric uniform state, D12 permits:
- two hidden-sector rates gamma_1,gamma_2;
- four retained-sector rates gamma_3,...,gamma_6.

The stationary heat-bath model happens to choose them all equal at the
microscopic label-refresh level.

That equality is not a consequence of the static operator spectrum.

---

## 4. Localized state

A localized state reduces the stabilizer from D12 to a smaller subgroup
(reflection for the aligned representative).

Therefore symmetry constraints become weaker, not stronger.

Fisher-44 explicitly showed phase/hidden mixing and hidden anisotropy there.

So localization cannot recover full tangent isotropy by symmetry alone.

---

## 5. Refinement does not help

Metric-edge refinement constrains spatial/edge scaling laws such as

    c(ell) proportional 1/ell

and additive storage

    m(ell) proportional ell.

It does not identify the six internal Fourier irreducible representations of a
single twelve-label cell.

Therefore the intercell refinement theorem has no typed route to enforce

    gamma_1=...=gamma_6

inside the strict intracell state space.

---

## 6. Exact source boundary

The reset generator can now be characterized very sharply:

    heat bath
    = reversible/local refresh
      + full centered-space kinetic isotropy.

The unresolved physical question is precisely the source of that kinetic
isotropy.

Possible future source principles include:
- a permutation symmetry larger than D12 at the microscopic update level;
- a maximum-entropy/minimum-information update theorem;
- a universal refresh mechanism acting on labels before the strict operator
  breaks the symmetry;
- a role-transfer theorem connecting a more symmetric legacy layer to the
  strict cell.

None is currently established.

---

## 7. Consequence for the one-FIN emergence programme

This result suggests a concrete possible hierarchy:

    deeper update law with larger kinetic symmetry
        ->
    equal microscopic relaxation
        ->
    strict operator / localization breaks static symmetry
        ->
    state-dependent Fisher and memory geometry emerge.

This would allow kinetic simplicity at the deepest layer and anisotropic
effective physics after symmetry breaking.

But this is a hypothesis until the larger-symmetry update law is actually
sourced.

---

## 8. Next atom

### PRE-STRICT-KINETIC-SYMMETRY-48

Search the legacy/pre-strict lane for an admitted permutation-equivariant
microscopic object that can act before D12/strict differentiation.

Test whether:
1. S12-equivariant one-label Markov dynamics uniquely gives reset/isotropic
   kinetics up to rate;
2. the strict operator can then emerge in statics while retaining or
   predictably breaking that kinetic isotropy;
3. this supplies any real progress on the open legacy->strict role-transfer
   gate.

Acceptance:
- a typed bridge from a more symmetric update layer into strict FIN; or
- a no-go showing that the legacy/strict objects are not related strongly
  enough to transfer kinetic symmetry.
