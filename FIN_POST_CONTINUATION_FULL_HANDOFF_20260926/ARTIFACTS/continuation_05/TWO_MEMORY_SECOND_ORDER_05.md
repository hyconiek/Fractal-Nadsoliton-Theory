# TWO-MEMORY-05 — second-order obstruction and cubic+quartic repair

Status: **CONDITIONAL IDENTIFIABILITY / NO-GO PAIR**.

The first-order shifted quartic tomography reconstructs hidden k=1,2 only as the
combination entering the symmetric jump kernel.  Because (f_j-f_i)^4 is
symmetric under i<->j, first order depends on p'+q'.

A bounded exact counterexample shows that first-order quartic data alone do not
fix the O(epsilon^2) constant.  For a hidden k=2 direction h, compare:
A) p=softmax(epsilon h), q=u;
B) p=q=softmax(epsilon h/2).
The complete first-order shifted quartic responses are identical.  Yet the
12-shift averaged second-order pure-k5 coefficient is 0 in A and 3/64 in B for
a unit cosine or sine k=2 direction.  k=1 gives zero in both examples.

There is, however, a target-blind repair using an odd jump statistic.  For the
same mixed retained probe f_t=c3+t c4:

quartic departure-channel coefficients:
  k1: t(9 t^2+10)/4,   k2: 3 t^2/4;
cubic departure-channel coefficients:
  k1: -3 t^2/8,        k2: -3 t/4.

Quartic first order reconstructs p'+q'; cubic first order reconstructs p'-q'.
Parity in t separates k1 from k2 and the 12-shift DFT separates cosine/sine.
Therefore all four hidden components of p' and q' are individually recovered.

An exact fixture with departure hidden coefficients (2,3,5,7) and target
(-1,4,6,-2) recovers both tuples exactly.  Once p'_2 and q'_2 are known, the
shift-averaged O(epsilon^2) cross contamination of an unnormalised pure c5
quartic is exactly

  (3/16) (p'_{2c} q'_{2c} + p'_{2s} q'_{2s}).

For the fixture this equals 3 exactly.  Single-channel p'' and q'' terms drop
out of the uniform shift average because they have zero total mass and the
shift-averaged single-channel kernel is translation independent.

Boundary: this second-order repair assumes the coupling enters through smooth
normalised departure/target distributions on the declared jump kernel.  More
general changes of the jump law can introduce new second-order objects.  No
physical coupling strength or clock is sourced.
