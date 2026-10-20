# SL-002/003: sampled histories, weak paths and an exact TV kill test

Date: 2026-10-20. Use the fixed finite source and A from
[SCALING_AND_GENERATOR_PROOF.md](SCALING_AND_GENERATOR_PROOF.md).
The main preparation is exactly pi for both source and limit. A separately
supplied initial mismatch epsilon_c costs an additional initial TV term
once; it is not multiplied by the number of observations.

## 1. Finite semigroups and the entire sampled record

Finite conservative Q_c,Q have stochastic semigroups, hence induced
infinity norm one. Their construction by uniformization at rate A is an
auxiliary proof of the already derived generator law, not an extra physical
Poisson bath added to the finite parent. The exponential series is a convex
mixture of stochastic powers.

Let h=1/c. Since P_c=I+hQ_c, Taylor's integral remainder and contraction give

    ||exp(hQ_c)-P_c||_infinity<=h^2||Q_c||^2/2<=2A^2/c^2.

Differentiating the matrix product of the two semigroups gives the direct
Duhamel bound, without commuting generators,

    ||exp(hQ_c)-exp(hQ)||_infinity
       <=h||Q_c-Q||_infinity<=4A^2/c^2.

Consequently, for R_h=exp(Q/c),

    ||P_c-R_h||_infinity<=6A^2/c^2,
    sup_x TV(P_c(x,.),R_h(x,.))<=3A^2/c^2.             (1)

Maximal row coupling at successive steps bounds **the whole sampled
record**, not just its endpoint: through n=floor(cT) transitions,

    TV(Law(Y_0,...,Y_n),Law(Z_0,Z_h,...,Z_(nh)))
      <=min(1,epsilon_c+3A^2 T/c).                      (2)

Here Y has the actual source kernel P_c; Z has the derived CTMC generator
Q. Their marginals remain correct after a mismatch by continuing with
their original transition rows. Thus no independence of observations is
assumed. With a common initial law the corresponding operator estimate is

    sup_(0<=t<=T)||P_c^floor(ct)-exp(tQ)||_infinity
       <=(6A^2 T+2A)/c.                                (3)

The power difference is a telescoping sum of stochastic left/right
factors, of norm at most n times the one-step bound; no commuting of P_c
and R_h is needed. The final term pays t-floor(ct)/c. For the primary pi preparation both
processes are stationary. Uniform row bounds also allow any separately
supplied common initial law of these **hypothesis kernels**; that is not
a transfer from the unchanged native iid-start FIN process.

## 2. A quantitative full-path time-warp coupling

Let X_c(t)=Y_floor(ct), a cadlag path on [0,T]. Give the finite state set
the discrete metric. Use the Skorokhod J1 time-warp distance obtained by
minimizing max(sup|lambda(t)-t|,sup d(X(t),Z(lambda(t)))) over increasing
continuous bijections fixing 0,T. This supplies the usual J1 topology;
no physical interpretation is assigned to a time warp.

Couple the source and CTMC grid skeletons as in (2). One can first generate
each CTMC segment, then sample the coupled source endpoint using the
finite conditional endpoint coupling; each marginal remains its original
law. No physical feedback or source intervention is involved.

The CTMC's outgoing rates are bounded by A. In an interval of length h,
the probability of two or more actual jumps is at most (Ah)^2/2, since
the actual jumps are dominated by uniformization attempts and the Poisson
second factorial moment is (Ah)^2. Summing over n cells costs <=A^2T/(2c).
Also exclude any CTMC jump after max(0,(n-1)h), an interval of length at
most 2h; its probability is <=2A/c. This pays both the last full cell and
any endpoint remainder, including the case T is an exact grid point.

On the remaining event the skeletons agree, every earlier cell has at
most one CTMC jump, and neither path changes near the endpoint. Map each
source jump time kh to its corresponding CTMC jump in ((k-1)h,kh], and
interpolate linearly, fixing 0,T. The map is strictly increasing and its
displacement is at most h. The state paths agree after this warp. Thus

    Pr(d_J1,T(X_c,Z)>1/c)
      <=min(1,epsilon_c+[(7/2)A^2T+2A]/c).              (4)

For n=0 the source is constant and the same endpoint estimate suffices.
For a bounded Lipschitz path functional, this coupling bounds its expectation
error by its Lipschitz constant/c plus twice its sup norm times (4).
Continuity/local boundedness at the fixed target law gives weak J1 path
convergence. There is no extrapolation from a trajectory CLT to spectral
gaps, filter convergence or infinite time.

In the actual source family c=L(q-1), (4) concerns the retained finite
coordinates of the supplied finite global entropy source. Its exact lumping
was paid before taking L->infinity. The visible carrier and pi stay fixed;
this is not the interacting full-N/g mean-field FIN limit.

## 3. Full continuous-path TV does not converge

Let q_x=sum_(y!=x)Q_xy. Source nonconstant paths jump only at the finite
grid times kh. A finite-rate CTMC has zero probability of a jump at any
specified deterministic grid time. Therefore the jumped parts of these
two path measures are mutually singular. Their constant paths are shared
atoms. With initial laws mu_c,mu, the exact full-path TV is

    TV(Law(X_c),Law(Z))
      =1-sum_x min(mu_c(x)P_c,xx^floor(cT),
                   mu(x)exp(-q_x T)).                 (5)

This is not a claim that the **entire** path measures are singular: their
no-jump atoms overlap. To prove (5), partition into constant paths, the
source's nonconstant grid-supported paths and the remaining CTMC paths;
only the first part can carry a common submeasure. On finite-state grid
paths this is also the elementary sum of minimum atom masses.

If mu_c->mu, Q_c->Q and c->infinity, then

    P_c,xx^floor(cT)->exp(-q_xT),
    TV(Law(X_c),Law(Z))->1-sum_x mu(x)exp(-q_xT).        (6)

For a nontrivial connected source and T>0 this limit is strictly positive.
Thus (2) and (4) must not be renamed complete continuous-path TV convergence.

An exact control is K=ones on two states, pi=(1/2,1/2), T=1. Then
P_c has diagonal (c+1)/(c+2), off-diagonal 1/(c+2), and Q has unit flip
rates. For integer c>=1,

    [(c+1)/(c+2)]^c>exp(-1),
    full-path TV=1-exp(-1),                             (7)

even though sampled-record TV tends to zero and the time-warp coupling
converges. The inequality uses log(1+x)<x, not floating quadrature. This
ideal timestamp-sensitive distinction is not an apparatus or sample-budget
claim. No posterior/information/sufficiency limit follows from weak J1
alone; only the separately paid grid experiment has the TV guarantee (2).
