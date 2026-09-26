# STRICT-A-CYCLIC-GENERATOR-185
## The strict C12 Laplacian is exactly a polynomial in one cyclic permutation

Date: 2026-09-26

Status:
exact algebraic identity in the declared strict vertex basis.

Let P be the one-step cyclic shift on the 12 internal FIN labels.

The strict shell weights are


    w_1
      = 0.469985672645020

    w_2
      = 0.192043551690103

    w_3
      = 0.091428614277925

    w_4
      = 0.047029168745650

    w_5
      = 0.024131223363630

    w_6
      = 0.011070817321442


For d=1,...,5 define the shell Laplacian

    L_d
      =
      2I-P^d-P^(-d).

For the antipodal shell d=6 there is only one distinct neighbor, so

    L_6
      =
      I-P^6.

Then the strict operator satisfies exactly

    boxed:
    A_strict
      =
      sum_(d=1)^5 w_d L_d
      +w_6 L_6.

The numerical reconstruction residual is at roundoff level.

Thus the complete strict internal Laplacian is generated algebraically by one
cyclic permutation P and its inverse.

## Interpretation boundary

This is an internal C12 carrier statement.

It does not prove that physical multicell units are ordered by the same P.

But it shows a mathematically striking compatibility:

    one information-preserving cyclic transformation
      -> cycle skeleton
      -> all cyclic distance shells
      -> strict-type radial Laplacian.

So if an external/multicell P were independently sourced, a strict-like
interaction family could be built as a function of transformation distance
without separately inserting coordinates.


Direct reconstruction residual:

    ||A_polynomial-A_direct||_F
      = 5.875e-16.
