# Exact short-lag two-time fingerprints on the same FIN strict kernel

Date: 2026-10-19. Purpose: exhibit **explicit, mathematically incompatible time laws** over the same published `A=sI-W`, with no new intrinsic source claim. This is a concrete corollary of FIN's known coherent/diffusive distinction; not a new claimed foundational result.

For each `i != j`, strict `W_ij>0`. Define a stochastic matrix `P(t)=exp(-tA)` and a unitary matrix `U(t)=exp(-itA)` from the *same* real symmetric operator A. Then Taylor expansion at t=0 gives

```
P_ij(t) = W_ij t + O(t^2)                 (i != j),
|U_ij(t)|^2 = W_ij^2 t^2 + O(t^4)         (i != j).
```

The O(t^4) remainder follows from reality/symmetry of A: `U_ij(t)` has an odd purely imaginary part and an even real part, hence its squared magnitude is even in t. For i=j, `P_ii(t)=1-s t+O(t^2)` while `|U_ii(t)|^2=1-t^2 sum_{j != i}W_ij^2+O(t^4)`.

Thus the off-diagonal transition probability has **linear** small-time scaling in the Markov realization and **quadratic** scaling in the coherent amplitude-realization followed by *an expressly additional Born-type endpoint readout*. This is not merely common-clock rescaling by a finite constant. A time reparametrization singular at zero might alter exponents, but would require its own clock/consistency premises.

`|U(t)|^2` is one-time stochastic when considered as entrywise moduli squared (it is doubly stochastic); it is **not** in general a Markov semigroup satisfying `B(t+s)=B(t)B(s)` and must not be promoted to a history law without an instrument specification. Exactly this additional measurement choice is why this diagnostic is not a derivation of QM or physical time from FIN.

This is a minimal, explicit **model discrimination** available from the strict kernel after the categories of evolution/readout are declared. The kernel by itself still does not decide which is nature's rule.
