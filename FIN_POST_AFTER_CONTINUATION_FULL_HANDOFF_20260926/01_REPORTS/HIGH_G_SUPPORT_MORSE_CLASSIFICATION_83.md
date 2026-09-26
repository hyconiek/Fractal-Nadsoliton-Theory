# HIGH-G-SUPPORT-MORSE-CLASSIFICATION-83
## Support size controls the asymptotic Morse index of concentrated FIN stationary branches

Date: 2026-09-26

Status:
- exact conditional asymptotic theorem;
- application to known branches is partly numerical because branch convergence
  to each proposed support must be established separately.

## 1. Stationary equation

For the retained finite-dimensional FIN model,

    p_g = softmax(g A7 p_g),

and

    theta_g = g X^T p_g.

Assume a stationary branch has a limit

    p_g -> p*

with support S.

## 2. Equal-field condition on the support

For any two labels i,j with positive limiting weights,

    log[p_i(g)/p_j(g)]
      =
      g[(A7 p_g)_i-(A7 p_g)_j].

The left side remains O(1).

Therefore

    (A7 p*)_i=(A7 p*)_j

for all i,j in S.

Thus p* is determined by:
- equal A7 field on S;
- normalization;
- positivity.

For a candidate support this is a finite linear system.

## 3. Outside-support inequality

For j not in S, a genuine exponentially concentrating branch requires

    (A7 p*)_j < m,

where m is the common field on S.

A strict gap gives exponential suppression of outside labels.

## 4. Hessian limit

The dual Hessian is

    H_g
      =
      I/g - X^T[
        diag(p_g)-p_g p_g^T
      ]X.

Let

    C_* =
      X^T[
        diag(p*)-p*p*^T
      ]X.

Then

    H_g -> -C_*

on directions where C_* is nonzero, while directions in ker(C_*) retain the
positive `1/g` regularization.

Hence, for sufficiently large g,

    boxed:
    Morse index(H_g)=rank(C_*).

## 5. Affinely independent support

If the feature vectors

    {X_i : i in S}

are affinely independent, then

    rank(C_*)=|S|-1.

Therefore

    boxed:
    asymptotic Morse index
      = |S|-1.

This gives a direct branch interpretation:

    one-label concentration
      -> index 0 minimum;

    two-label concentration
      -> index 1 saddle;

    three-label concentration
      -> index 2 saddle;

    four-label affinely independent concentration
      -> index 3 saddle;

and so on.

## 6. Existing numerical examples

The current continuation campaign already displays this pattern:

### Main localized branch
high-g concentration on one dominant label:
    index 0.

### Generic i=5 branch
three-label support:
    index 2.

### Reflection branch from atlas i=1
numerically approaches a two-label 1/2-1/2 mixture:
    index 1.

### Several higher-index reflection / pure-k4 branches
numerically approach four-label mixtures:
    asymptotic index 3.

These examples are consistent with the theorem.

## 7. Physical interpretation boundary

This is a theorem about the geometry of the FIN stationary exponential-family
landscape.

A support label is not thereby a physical particle, location, or field state.

The useful structural statement is narrower:

> concentration complexity and saddle index are linked.

The number of independently competing surviving labels determines the number
of negative covariance directions.

## 8. Research consequence

The high-g stationary atlas can now be classified by candidate supports rather
than by blind continuation alone.

For each branch:
1. infer the dominant support numerically;
2. solve the equal-field support equations;
3. check the outside-field gap;
4. predict the asymptotic Morse index;
5. then validate the branch continuation.

This should turn the large-g frontier into a finite combinatorial search over
self-consistent supports.
