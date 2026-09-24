# NL-12 — is FIN q-refinement the Markov category needed by Čencov?

Status: **NO_GO_NAIVE_FULL_Q_REFINEMENT_IS_NOT_A_CONGRUENT_MARKOV_EMBEDDING**

There is an exact obstruction before any probability fitting is attempted.

Under the natural carrier subdivision `Z_q -> Z_(2q)`, old node `j` is fine
node `2j`.  A fixed low Fourier character restricts correctly:

`exp(2 pi i k (2j)/(2q)) = exp(2 pi i k j/q)`.

Hence the smooth `k=3,4,5` modes are refinement-compatible.

The parity coordinate is different.  The fine model chooses its own Nyquist
mode `k=q` on `Z_(2q)`.  On the old nodes,

`exp(2 pi i q (2j)/(2q)) = 1`.

But the coarse model's parity is

`exp(2 pi i (q/2)j/q)=(-1)^j`.

So the "Nyquist parity at each resolution" is **not the same character under
refinement**.  This is an exact representation-level mismatch.

This explains why a naive q-to-2q probability aggregation is not an exact
congruent embedding for the full localized branch.  The large-q Fisher metric
may converge extremely well, but convergence is weaker than the categorical
hypothesis required for a Čencov uniqueness argument.

Possible repairs are now explicit:
1. treat the Z2 fiber as an independent persistent variable instead of
   re-identifying it with the new Nyquist mode at each q;
2. drop parity before defining the refinement category;
3. construct a different explicit Markov/refinement morphism and prove that it
   preserves the relevant FIN family.

Until one of these passes, Fisher cannot be promoted from candidate metric to a
FIN-sourced kinetic law.
