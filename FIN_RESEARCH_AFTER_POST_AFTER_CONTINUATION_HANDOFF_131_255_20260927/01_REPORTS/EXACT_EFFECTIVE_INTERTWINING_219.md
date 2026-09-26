# EXACT-EFFECTIVE-INTERTWINING-219
## The 12-state effective chain projects exactly to the Z3 chain at the level of the full semigroup

Date: 2026-09-26

Status:
exact theorem for any D12-circulant localized-state generator;
numerical residuals are roundoff checks for N=3..8.

Let Q12 be the effective localized-basin generator with shell rates q_d.

Define the lift

    (R f)(j)
      =
      f(j mod 3).

Define

    k
      =
      q1+q2+q4+q5

and

    Q3
      =
      k
      [ -2  1  1
         1 -2  1
         1  1 -2 ].

Then, because every state in one residue class has exactly the same total rates
to each other residue class,

    boxed:
    Q12 R
      =
    R Q3.

This is the strong intertwining condition that was missing from the earlier
naive refinement test 132.

## Semigroup consequence

By induction,

    Q12^m R
      =
    R Q3^m

for every m>=0.

Therefore

    boxed:
    exp(t Q12) R
      =
    R exp(t Q3)

for every t.

So the Z3 variable is an EXACT Markov lumping of the 12-state effective chain,
not merely a first-derivative closure.

## Numerical replay

For N=3,...,8 the matrix residual

    ||Q12 R-R Q3||

is at machine precision.

Powers through m=6 also agree at roundoff.

## Importance

This is the first clear positive realization of the strengthened 132 criterion
inside the current campaign.

The hierarchy is:

    exact microscopic chain
      --memory-aware approximate reduction-->
    12 localized-state Markov chain
      --exact intertwining-->
    3-state Z3 Markov chain.

Thus the only paid approximation occurs in the microscopic-to-basin reduction.

The basin-to-Z3 reduction itself is dynamically exact once the circulant
effective chain has been obtained.

## Boundary

This does NOT prove exact intertwining directly from the microscopic chain to
Z3.

Reports 141-144 show that microscopic elimination carries short memory.

The exact semigroup result begins only after that memory has been absorbed into
the effective 12-state generator.
