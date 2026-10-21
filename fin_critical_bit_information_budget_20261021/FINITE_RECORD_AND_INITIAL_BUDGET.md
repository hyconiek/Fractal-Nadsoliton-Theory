# Actual finite-parent bits: the prefix resource cannot be omitted

Date: 2026-10-21. Conditional analytical source-interface inequality.

Fix critical lag t>0 and M transitions. Observe the unchanged iid parent at
model times tau_N+sqrt(N)jt, j=0,...,M, with the actual TR polynomial
origin. Let epsilon be its path-TV comparison bound to stationary local
reflection through model horizon sqrt(N)Mt. In particular TR supplies
epsilon->0 when sqrt(N)Mt<=exp(aN), every fixed0<a<b. The sqrt(N) clock
factor is not dropped.

The reflected binary block has entropy >=M c_N(t). A maximal coupling of
the two complete blocks has error <=epsilon. Fano on the binary-word
alphabet, of cardinality2^(M+1), gives for epsilon<=1/2

    H_actual(Y_0,...,Y_M)
      >=M c_N(t)-epsilon(M+1)log2-b2(epsilon).             (1)

Thus no unbounded-information transfer is taken from TV for free: the
finite alphabet and its actual record length pay the continuity term.

Now consider an ACTUALLY SUPPLIED competing source block Yhat and a joint
law with mean bit mismatch e. Suppose its observable entropy is paid by

    H(Yhat_0,...,Yhat_M)<=I_N+H M t.                      (2)

This is a premise about that source's finite block, not a consequence of
its KS entropy bound. It may be proved from a particular seed/history or
regularity contract, but none is invented here. Conditional entropy chain
rule and concavity give H(Y|Yhat)<=(M+1)b2(e). Combining with (1),

    I_N/M+Ht+(1+1/M)b2(e)
      >=c_N(t)-epsilon(1+1/M)log2-b2(epsilon)/M.           (3)

For M_N->infinity, I_N/M_N->0, a common finite H in CRITICAL time units,
and the stated horizon, (2)–(3) yield the corresponding positive error
floor at any t with c(t)>Ht. Full infinite stationary joining is not
assumed for this finite-record version; the block budget and coupling are
instead explicit. Error independence and a Markov readout are not used.

## BI-005. A concrete closed finite-carrier obstruction

Consider ANY competing deterministic source whose complete initial random
state has at most 12^N possibilities. Allow an arbitrary initial law and
arbitrary fixed time-dependent deterministic readout maps. With no fresh
random input or additional random initial parameter, data processing gives

    H(Yhat_0,...,Yhat_M)<=H(initial state)<=N log12.        (4)

This is an input class, not a claim that the primordial nadsoliton has
finite cardinality. It includes an attempted deterministic reinterpretation
on the same N labelled q12 slots, without additional initial information.

Take the single frozen sequence M_N=N^2 and fixed critical lag t>0.
Its model-time record horizon sqrt(N)N^2 t is polynomial and hence falls
within the paid TR window for all sufficiently large N. Substituting
I_N=N log12 and H=0 in (3), then using c_N(t)->c(t)>0, gives

    liminf_N b2(e_N)>=c(t),
    liminf_N e_N>=b2^(-1)(c(t))>0.                       (5)

This concerns actual unpostselected iid-parent records. No full stationary
source joining, finite-H KS premise, Markov source, independent error,
copy symmetry or restriction on computational complexity is required.
A fixed infinitely precise coefficient in a deterministic law does not
add ensemble entropy; a randomly chosen coefficient/history would be
additional initial information and must be counted.

More generally a deterministic source matching these blocks with e_N->0
must satisfy liminf log(cardinality(initial state))/N^2>=c(t), if its
initial entropy is paid solely by that finite cardinality. This is a
necessary capacity bound, not a construction or an assertion of physical
finite precision. A continuous or richer initial state may pay more
information; it is not excluded by the finite-cardinality corollary.

## A zero-KS control with a large initial resource

A fixed binary de Bruijn cycle of order L has period2^L. Uniform random
phase yields B fair independent bits for every B<=L, hence block entropy
B log2, although its infinite deterministic periodic extension has entropy
rate zero. Its phase entropy is L log2 and pays those finite blocks.
The verifier uses L=5, not a large-word search.

This illustrates the old FH limit-order warning, not a new counterexample
to ID. If a source randomly selects a whole long history, its initial
budget includes that choice, not merely the phase of one frozen cycle.
The accepted FH ensemble can carry such a large resource; it is NOT
refuted by omitting it from (2). Neither finite-window agreement nor
infinite-dimensionality is an independently derived natural source law.

## Scope

N first/fixed-lag convergence and positive coarse-lag rounding are paid.
No small-time entropy expansion uniform in arbitrary t_N->0 is asserted.
No practical native entropy/error at an example N, numerical N0, optimal
source or finite-H property of the nadsoliton is concluded. Every source
must bind its own carrier/history, readout, preparation, clock map and
initial-information budget. F01 remains OPEN_SOURCE.
