# FH: finite admissible histories for the already authorized hypothesis

Date: 2026-10-20. One analytical route, performed by the same author without
delegation at the user's request. Not formal verification or peer review.

Model: **FIN-HYP-HIST-ME-001**, unchanged. This package constructs a finite
ensemble representation of that hypothesis; it does not adopt uniform
history counting as a new physical law or replace primary FIN heat bath.

The finite state set S has size m. K is symmetric, 0/1, connected, with
exactly one self continuation K_xx=1. Supply pi>0, sum pi=1. Reuse the
unique positive scaling and its selected transition law

    pi_x=z_x(Kz)_x, c*_xy=K_xy z_x z_y,
    P_xy=c*_xy/pi_x, s_x=P_xx=z_x^2/pi_x.

For positive integer visit counts n summing to L, a cycle is a **rooted
word** (x_0,...,x_(L-1)), with x_L=x_0, every K edge legal, and exactly
n_x departures from x. All distinct words count once, including different
root positions when their words differ. Do not divide by rotation-orbit
sizes: periodic words can have nontrivial stabilizers. U_n is uniform on
this finite set if it is nonempty; an empty set defines no probability law.

For irrational pi use explicitly stated rational occupation approximations;
do not claim exact equality L pi=n. At native size N the histories are of
labelled configurations in (Z/qZ)^N, not uniformly counted count paths.
The global Gibbs pi_N, A7, gain, reference, one-self convention and initial
local preparation remain supplied premises. The cycle length is an analysis
parameter, not a physical clock or cosmological periodicity.

Targets: exact completion-count sampler; a counting lower bound; whole-block
KL/TV, not just empirical flows; a paid native diagonal limit using existing
IS/AD estimates. Preserve the HB mismatch and all physical/source gates.
