# EXACT-FINITE-N-GIBBS-HEAT-BATH-61
## R7P-111 admits an exact single-label maximum-entropy Gibbs sampler

Date: 2026-09-26

Status:
- exact finite-N conditional-distribution theorem;
- exact detailed balance for the labelled-copy Gibbs measure;
- exact occupation-count Markov generator;
- converges to the declared self-consistent heat-bath target at large N.

## 1. Finite-copy measure

Take the R7P-111 labelled-copy law

    Pi_N(x_1,...,x_N)
      proportional to
      12^(-N)
      exp[
        (g/(2N))
        sum_{a,b}
        A7[x_a,x_b]
      ].

The rank-seven operator A7=X7 X7^T is symmetric and circulant, so its diagonal
A7[j,j] is independent of j.

## 2. Exact conditional target

Hold every label except label a fixed.

Let m_j be the counts of the other N-1 labels.

Terms depending on a candidate new value j are

    (g/N) sum_k A7[j,k] m_k
      + (g/(2N)) A7[j,j].

The diagonal term is independent of j and cancels from normalization.

Therefore

    boxed:
    q_j^{(-a)}
      =
      softmax_j[
        (g/N)(A7 m)_j
      ].

This is the exact single-label Gibbs conditional.

## 3. Maximum-entropy update

By report 54, for this fixed conditional target the maximum-conditional-entropy
one-event update is to forget the old label and sample the new one independently
from `q^{(-a)}`.

Thus an exact finite-N heat-bath/Gibbs sampler exists:

    choose label a,
    compute the field from the other N-1 labels,
    redraw x_a from q^{(-a)}.

No Metropolis acceptance step is required.

## 4. Detailed balance

If two microstates x and x' differ only in label a, the Gibbs conditional gives

    Pi_N(x) q_{x'_a}^{(-a)}
      =
    Pi_N(x') q_{x_a}^{(-a)}.

Therefore the single-label refresh chain is reversible with respect to Pi_N.

Since all conditional probabilities are strictly positive, repeated updates
make the finite state chain irreducible.

So Pi_N is its unique stationary law.

## 5. Occupation-count generator

At count state n with total N, a label currently of type i is chosen at rate 1.

After removing that label, the remaining count vector is

    m=n-e_i.

The exact off-diagonal count transition rate is

    boxed:
    n -> n-e_i+e_j

at rate

    n_i q_j^{(i)}(n),

where

    q^{(i)}(n)
      =
      softmax[
        (g/N) A7(n-e_i)
      ].

The count process is Markov because the conditional depends only on counts.

Its invariant count law is exactly the multinomial degeneracy times the
R7P-111 Gibbs tilt.

## 6. Relation to the infinite-N heat bath

Write p=n/N and

    q(p)=softmax[g A7 p].

Then

    q^{(i)}(n)
      =
      softmax[
        g A7 p
        -(g/N)A7 e_i
      ].

Hence

    q^{(i)}(n)
      = q(p)+O(1/N)

uniformly on compact interior state sets.

The count drift is

    dot p
      =
      sum_i p_i q^{(i)}(n)-p,

so

    dot p
      = q(p)-p+O(1/N).

Therefore the previously declared deterministic heat-bath equation is the
large-N limit of an exact finite-copy reversible Gibbs sampler.

## 7. Source boundary

This closes a consistency gap between:
- the conditional finite-copy Gibbs realization;
- the maximum-entropy reset rule;
- the deterministic heat-bath mean equation.

Still supplied:
- g;
- X7 / rank-seven mediator choice;
- the per-label update clock rate;
- any physical interpretation of the labels.
