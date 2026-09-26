# CRT-Z12-FACTORIZATION-204
## The strict internal carrier factorizes exactly as Z12 ≅ Z3 × Z4

Date: 2026-09-26

Status:
exact finite-group / matrix identity.

Because gcd(3,4)=1,

    boxed:
    Z12 ≅ Z3 × Z4

through the Chinese remainder map

    j
      ->
    (a,b)
      =
    (j mod 3, j mod 4).

The inverse map is

    j
      =
      4a+9b
      mod 12.

## 1. One-step cyclic transformation

Let P12 be the unit shift

    j -> j+1.

Under the CRT basis permutation:

    boxed:
    P12
      =
      P3 tensor P4.

The numerical residual is exactly zero.

## 2. Independent factor shifts already live inside P12

Since

    4 ≡ 1 mod 3,
    4 ≡ 0 mod 4,

we have

    boxed:
    P12^4
      =
      P3 tensor I4.

Similarly,

    9 ≡ 0 mod 3,
    9 ≡ 1 mod 4,

so

    boxed:
    P12^9
      =
      I3 tensor P4.

Thus the single cyclic carrier contains two commuting coordinate shifts.

## 3. Strict operator quotient subspaces

The subspace of functions depending only on

    j mod 3

is exactly invariant under strict A.

Its compressed operator is

    A3
      =
      c3 L_C3,

with

    c3≈0.733189616444

and spectrum

    {0,
     2.199568849333,
     2.199568849333}.

The invariance residual is at machine precision.

The subspace depending only on

    j mod 4

is also exactly invariant.

Its compressed operator is a weighted C4 Laplacian with:
- nearest-neighbor weight ≈0.58554551;
- opposite weight ≈0.39515792;

and spectrum

    {0,
     1.961406861976,
     1.961406861976,
     2.342182041146}.

## 4. Meaning

The strict carrier is not merely a twelve-point cycle.

Algebraically it admits an exact:

    base Z3
      ×
    fiber Z4

coordinate factorization.

No physical-space interpretation is made here.
