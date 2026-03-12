# N464 Current First Strict QW‑2191 Site‑Distribution Permutation‑Invariant Objective O(2)‑Cut Nonuniqueness Boundary Theorem

Status: `N464_DISCHARGED_CURRENT_FIRST_STRICT_QW2191_SITE_DISTRIBUTION_PERMUTATION_INVARIANT_OBJECTIVE_O2_CUT_NONUNIQUENESS_BOUNDARY_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Close, at theorem level (not probe level), a broader strict claim pattern that continues to recur around `T165` after `N463`:

```text
“Use an information-theoretic objective of squared site-amplitude distributions (entropy / KL / divergence / etc.)
to select a unique theta on the QW-2191 O(2) family.”
```

`N463` already closes the **Shannon entropy** special case.

`N464` packages the stronger, still honest, current-state boundary:

```text
any objective that depends only on the squared site-amplitude probability distributions up to site permutation
(in particular, is invariant under ring translations)
inherits a nontrivial theta periodicity and therefore cannot yield a unique O(2)-cut.
```

This theorem does **not** deny that a future strict Shannon‑typed selector could exist.
It closes only this large permutation‑invariant / translation‑invariant objective family as a strict uniqueness source.

## Strict-admissible evidence reused

1. `QW-2190`
   - fixed 12‑octave ring scaffold and real Fourier mode basis on the site index `x=0..11`.
2. `QW-2191`
   - continuous `O(2)` rotation family in each degenerate two‑mode subspace.
3. `N463`
   - the explicit “translation induces theta periodicity” computation (reused as a lemma pattern).
4. `P423` (probe-level only)
   - numerical non-uniqueness seen for several naive objectives (non-theorem evidence only).

## Setup (rotated Fourier modes and squared site distributions)

Let `n=12` and define the unit-norm real Fourier pair in `R^n`:

```text
c_m(x) := sqrt(2/n) * cos(2π m x / n),
s_m(x) := sqrt(2/n) * sin(2π m x / n),
```

for `m ∈ {1,2,...,(n/2-1)}` and `x ∈ {0,1,...,n-1}`.

For a rotation angle `θ ∈ R`, define:

```text
u_m,θ := cos(θ) c_m + sin(θ) s_m,
```

and the squared site-amplitude distribution:

```text
p_m,θ(x) := (u_m,θ(x))^2,   with  Σ_x p_m,θ(x) = 1.
```

Let `T_k` be the cyclic shift operator on site indices:

```text
(T_k f)(x) := f(x-k mod n).
```

## Theorem

### Claim 1. Ring translations induce nontrivial theta periodicity of the squared site distributions.

For each integer `k`:

```text
u_m,θ+2π m k/n  =  T_k u_m,θ.
```

Therefore:

```text
p_m,θ+2π m k/n  =  p_m,θ ∘ (x ↦ x-k),
```

so the distributions differ only by a permutation of site labels.

*Proof.* Same direct cosine computation as in `N463`:

```text
u_m,θ(x) = sqrt(2/n) * cos(2π m x/n - θ).
```
∎

### Claim 2. Any permutation‑invariant objective on `p_m,θ` cannot yield a unique minimizer in θ.

Let `F` be any functional on probability distributions on `{0,1,...,n-1}` such that it is invariant under site permutations.
In particular, it is invariant under cyclic shifts:

```text
F(p ∘ (x ↦ x-k)) = F(p)   for all integers k.
```

Define the objective:

```text
J_m(θ) := F(p_m,θ).
```

Then for all integers `k`:

```text
J_m(θ+2π m k/n) = J_m(θ).
```

So `J_m` is periodic with a nontrivial period `2π m/n`, hence cannot have a unique global minimizer modulo `2π`:
if a minimizer `θ*` exists then all values

```text
θ* + 2π m k/n,  for k=0..(n/gcd(n,m)-1),
```

are also minimizers and are distinct modulo `2π`.

Therefore `J_m` cannot supply a canonical strict-core `O(2)` cut on the `QW-2191` rotation family.
At best it can select a discrete orbit of size at least `n/gcd(n,m)` unless additional structure is introduced.

*Proof.* By Claim 1, `p_m,θ+2π m k/n` is a permutation of `p_m,θ`; by permutation invariance of `F`, the objective is equal. ∎

## Consequence (for T165 discipline)

This theorem closes, as a strict-core uniqueness route, *all* objectives of the form:

```text
θ ↦ F( p_m,θ )   where F is permutation-invariant on the site labels.
```

This includes (as special cases) Shannon entropy, Rényi entropies, and common divergence-to-reference objectives such as
KL divergences to the uniform distribution, provided they are built only from the unlabeled site distribution.

Therefore, discharging `T165` requires **additional strict structure** not erased by ring translations, e.g.:

1. a strict non-translation-invariant datum on the ring (a marked site, directed edge, or equivalent), or
2. an explicit quotient-safe mechanism proving that the residual discrete ambiguity is an admissible tracked residual symmetry,
   rather than a hidden selector slot.

`N458/N459` already show that `Omega_16` equipartition does not canonically induce such a nontrivial asymmetry on the phase quotient
orbit carrier without adding new symmetry-breaking/selection structure.

## What N464 does not prove

`N464` does not prove:

1. impossibility in principle of a future strict Shannon‑typed selector ingredient,
2. discharge of `T165`,
3. discharge of `QW-2191`,
4. strict-core theta export (`T159`) or strict-core selector closure,
5. ToE closure.

