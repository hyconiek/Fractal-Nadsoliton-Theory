# UPDATE-COMMUTATOR-GRAPH-201
## Noncommutation of local refresh operators also gives a complete dependency graph

Date: 2026-09-26

Status:
exact symmetry statement plus sparse-matrix replay for N=2,3.

Let K_a be the one-step kernel that refreshes labeled copy a while leaving the
others fixed.

If two local updates represented independent/noninteracting events, one would
expect

    K_a K_b
      =
    K_b K_a.

In the leave-one-out Gibbs process this generally fails because refreshing a
changes the conditional field used when b is refreshed, and vice versa.

Thus

    [K_a,K_b]
      !=0

for every a!=b.

## Replay

For N=2:

    ||[K_0,K_1]||_F / sqrt(|Omega|)
      ≈ 0.684069438653.

For N=3 all three unordered copy pairs give the same value by exchangeability:

    ≈ 0.298860528902.

No pair is distinguished as non-neighboring.

## Consequence

A dependency graph based on noncommuting update operations is again

    K_N.

So both operational definitions tested so far:

    intervention influence

and

    update noncommutation

derive the SAME qualitative incidence:

    complete weak coupling.
