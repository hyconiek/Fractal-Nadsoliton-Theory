# FINITE-WINDOW-PROCESS-ERROR-THEOREM-309
## Exact path-TV composition theorem, explicit rare-history accounting, and a held-out four-readout FIN validation

Date: 2026-09-27

Status:
- exact finite-alphabet theorem for arbitrary finite observation windows;
- direct exact microscopic four-readout training calculations for N=3,4,5;
- reversible rank-120 spectral calculations with explicit path-TV tail allowances for N=6 and held-out N=7;
- rank-140 convergence check for held-out N=7;
- the current t1=12 contract PASSES as a weighted/path-level finite-window theory but does NOT give a useful small strict history-uniform bound over every observed history;
- this is still an effective-theory result, not a derivation of a fundamental physical carrier, A7, g, space, QM or GR.

## 1. Research question

Reports 306-308 established that one frozen FIN effective law predicts one-, two-, and three-readout statistics on held-out N=7 after an operational preparation map and an initial-slip layer.

The remaining issue is logical rather than merely numerical:

    when do one-step conditional errors compose into a controlled error for
    an arbitrary finite observation window?

The projected microscopic process is not assumed to be Markov.

The effective prediction is the frozen 12-state hidden Markov model from 304-308, including the preparation/initial-slip prior and the declared readout channel.

## 2. Exact finite-window composition theorem

Let P and Q be arbitrary probability laws on a finite observed path

    (Y1,...,Ym).

They need not both be Markov.

For an observed history

    h_r=(y1,...,yr),

define

    p_r(.|h_r)=P(Y_{r+1} in . | H_r=h_r),
    q_r(.|h_r)=Q(Y_{r+1} in . | H_r=h_r),

and the history-dependent conditional defect

    delta_r(h_r)=TV[p_r(.|h_r),q_r(.|h_r)].

Then

    boxed:
    TV(P_{1:m},Q_{1:m})
      <=
    TV(P_1,Q_1)
      +
    sum_{r=1}^{m-1} E_{H_r~P}[delta_r(H_r)].

### Proof

For one extension step,

    P_{1:r+1}(h,y)=P_{1:r}(h)p_r(y|h),
    Q_{1:r+1}(h,y)=Q_{1:r}(h)q_r(y|h).

Add and subtract P_{1:r}(h)q_r(y|h):

    |P(h)p(y|h)-Q(h)q(y|h)|
      <=
    P(h)|p(y|h)-q(y|h)|
      +
    |P(h)-Q(h)|q(y|h).

Sum over h,y and divide by two:

    TV(P_{1:r+1},Q_{1:r+1})
      <=
    TV(P_{1:r},Q_{1:r})
      +
    E_{P_{1:r}} delta_r(H_r).

Iterating proves the theorem.

No Markov assumption is used.

## 3. Two corollaries: strict-uniform and high-probability

If for every history with positive microscopic probability

    delta_r(h)<=d_r,

then

    boxed:
    TV(P_{1:m},Q_{1:m})
      <=
    epsilon_1 + sum_r d_r,

where

    epsilon_1=TV(P_1,Q_1).

This is the strict history-uniform form.

However, a supremum may be dominated by a small exceptional set.

Let G_r be a declared good-history set such that

    P(H_r notin G_r)<=beta_r

and

    sup_{h in G_r} delta_r(h)<=d_r.

Since TV<=1,

    E delta_r
      <=
    (1-beta_r)d_r + beta_r.

Therefore

    boxed:
    TV(P_{1:m},Q_{1:m})
      <=
    epsilon_1
      +
    sum_r [(1-beta_r)d_r+beta_r].

The exceptional probability is explicit; no rare history is silently discarded.

The sharpest quantity available from a known microscopic path law remains the exact weighted defect

    bar_delta_r=E_P delta_r.

## 4. Four-readout validation protocol

Keep the frozen protocol of 308:

    t1 = 12,

and three subsequent intervals

    rho_pred Delta t = 0.5.

Preparation class:

    mu_(N,kappa)(n)
      proportional to
    pi_N(n) exp(kappa n0/N)

inside localized basin J=0.

For the four-readout validation grid use

    kappa = 0,1,...,12.

Readout noise remains

    eta=0.05.

No four-time parameter is fitted.

## 5. Exact/spectral microscopic path calculation

For N=3,4,5 the full four-readout microscopic law was computed directly:

    J_abcd
      =
    mu exp(Q t1)
      F_a exp(Q Delta t)
      F_b exp(Q Delta t)
      F_c exp(Q Delta t)
      F_d.

For N=6 and N=7 the reversible representation

    S=D_pi^(1/2) Q D_pi^(-1/2)

was used with the 120 largest eigenmodes.

For N=6:

    lambda_120 = -1.24892022644309,
    max eigenpair residual = 2.48e-14.

For N=7:

    lambda_120 = -1.14086606627006,
    max eigenpair residual = 1.30e-14.

A rank-140 N=7 rerun changes:

    max four-time path TV by 1.08e-10,
    weighted conditional defects by at most 4.88e-11,
    strict conditional suprema by at most 1.29e-9.

Thus the reported history-profile structure is numerically converged far beyond the scientific error scales below.

The conservative omitted-spectrum allowance is still retained for path-TV certification.

## 6. Training four-time envelope

Leave-one-N-out predictions on N=3..6 give the following four-readout path errors:

    N=3:  7.796616 % TV
    N=4:  3.914367 % TV
    N=5:  3.543169 % TV
    N=6:  3.262968 % TV spectral approximation
            <= 3.576512 % TV after spectral allowance.

Therefore the frozen training envelope is

    boxed:
    B_4time = 7.796616 % TV.

The maximum weighted chain bound over training is

    16.221994 %,

and the maximum strict-uniform chain bound is

    24.503958 %.

The theorem is deliberately conservative; every tested path error is below its weighted composition bound.

The training envelopes of the individual weighted conditional defects are

    step 1: 5.114993 %,
    step 2: 4.882069 %,
    step 3: 4.898567 %.

The corresponding strict suprema are

    step 1: 7.723128 %,
    step 2: 8.169698 %,
    step 3: 8.614458 %.

## 7. Held-out N=7 four-readout result

The frozen cross-N law gives

    rho_pred = 0.0125026797791533,
    Delta t = 39.9914265447069.

The rank-120 held-out four-time microscopic-versus-predicted error is

    2.759606 % TV.

The conservative omitted-spectrum allowance is

    1.103435 % TV.

Hence

    boxed:
    TV_true(N=7, four readouts)
      <=
    3.863040 %.

This is below the frozen training envelope

    7.796616 %

by

    boxed:
    3.933576 percentage points.

Thus the four-readout held-out validation PASSES.

## 8. History defects on held-out N=7

The rank-120 history-conditioned estimates are:

### Weighted defects

    step 1: 1.581224 %
    step 2: 1.214956 %
    step 3: 1.148425 %.

Together with the first marginal error they give a weighted chain bound of

    5.611384 %.

This is again above the actual rank-120 path error

    2.759606 %,

as required by the theorem.

### Strict history suprema

The corresponding maxima over all observed histories are much larger:

    step 1: 20.834410 %
    step 2:  7.003981 %
    step 3:  7.341084 %.

The resulting strict-uniform chain bound is

    36.519574 %.

This is mathematically valid as a bound on the truncated microscopic law but scientifically too loose to be called a useful uniform closure certificate.

The worst step-1 branch carries about

    1.50 %

of microscopic history mass.

Thus the poor strict supremum is not generated by a measure-zero branch, but neither is it representative of typical histories.

## 9. Explicit rare-history accounting

For N=7, optimizing the transparent split

    good-history conditional defect d
      +
    exceptional history mass beta

for each step gives representative worst preparation points:

### Step 1

    d = 5.176351 %,
    beta = 3.615717 %,
    d+beta = 8.792068 %.

### Step 2

    d = 2.417880 %,
    beta = 0.410231 %,
    d+beta = 2.828111 %.

### Step 3

    d = 2.391669 %,
    beta = 0.032180 %,
    d+beta = 2.423849 %.

The first transition after the preparation boundary remains the only materially nonuniform one.

The later history defects are both small on average and small outside explicitly tiny exceptional sets.

The strict N=7 conditional numbers are spectral estimates, not interval-certified conditional bounds; conditioning can amplify a small path-law perturbation.  Their rank-120/rank-140 numerical stability is nevertheless at the 1e-9 level.

## 10. Does an auxiliary memory coordinate remain necessary?

The answer is nuanced.

The current data do NOT support a persistent order-one hidden-memory defect:

- weighted next-step errors are about 1.1-1.6%;
- step-2 and step-3 strict suprema are about 7%;
- the four-readout path law passes held-out validation;
- rank-120/rank-140 history defects are numerically identical at the reported precision.

Therefore an additional persistent physical memory coordinate is NOT forced by 309.

However, the current t1=12 contract does NOT satisfy a strong small strict-uniform closure over every observed first history, because the step-1 supremum is about 20.8%.

Operationally, the frozen 12-state hidden-state posterior conditioned on the observed history remains the correct effective representation.  The seven-bin observation alone is not itself an autonomous Markov state.

Thus the strongest justified statement is:

    boxed:
    weighted/high-probability finite-window closure: PASS;
    strict small history-uniform closure at t1=12: NOT ESTABLISHED.

## 11. Exploratory burn-in scan

This scan was performed only AFTER the held-out t1=12 result and is therefore exploratory, not part of the 309 certification.

For N=7, the maximum first-step conditional defect behaves as:

    t1=8:   15.15 %
    t1=12:  20.83 %
    t1=16:  25.25 %
    t1=24:  10.29 %
    t1=32:   3.77 %.

The weighted defect remains around 1.3-1.7% throughout and falls to

    1.30 %

at t1=32.

Therefore the large strict step-1 defect at t1=12 does not look like persistent long memory.  A later burn-in may yield a genuinely useful uniform history bound.

This cannot be claimed yet because t1=32 was inspected after seeing the held-out N=7 data.

## 12. Comparator discrimination at four times

Give the comparator the SAME predicted state at t1 and change only the same-rho/same-exit transition generator.

The minimum predicted FIN-versus-comparator four-time separation is

    boxed:
    12.814146 % TV.

Subtracting the frozen training four-time envelope gives a pre-heldout margin

    boxed:
    5.017530 percentage points.

Using the actual held-out certified prediction error gives the stronger triangle lower bound

    boxed:
    TV(true microscopic FIN, comparator)
      >=
    8.951106 %.

At nominal 5% readout noise the minimum Chernoff information is

    C_min = 0.0149379161,

for which the standard equal-prior Chernoff bound falls below 5% at

    boxed:
    155 independent four-readout trajectories.

This is a model-versus-model count, not a complete laboratory sample-size guarantee.

## 13. Scientific verdict

309 separates three logically different statements.

### PASS: exact error-composition theorem

For arbitrary finite observed paths, path TV is bounded by the initial marginal error plus the microscopic-history-weighted one-step conditional defects.

This theorem does not assume that the microscopic projection is Markov.

### PASS: held-out finite-window validation

The frozen FIN effective law predicts the held-out N=7 four-readout path within

    <=3.8631 % TV

including the conservative spectral-tail allowance, comfortably inside the pre-existing training envelope.

### PARTIAL / NOT ESTABLISHED: strict uniform Markov closure

The present t1=12 protocol has a 20.8% worst first-history conditional defect.  Therefore a small strict history-uniform theorem has not been obtained.

The data instead support a weighted/high-probability process theorem with explicit rare-history accounting.

## 14. Epistemic boundary

This result strengthens the effective-theory lane only.

It does not establish:
- a fundamental simultaneous physical carrier;
- physical spatial incidence;
- a source for A7 or g;
- dimensional units;
- quantum mechanics;
- general relativity;
- a fundamental Markov law of nature.

It shows that one frozen effective FIN law, with an operational preparation map and an explicit observation model, controls increasingly long finite observation windows on a held-out finite-N microscopic process.

## 15. Next research target

The correct next P0 is not a fifth readout by brute force.

It is:

    310 — HISTORY-UNIFORM-LATE-BURN-CLOSURE

Training N=3..6 must be used to choose a burn-in criterion BEFORE opening the corresponding held-out N=7 test.

The goal is to determine whether there exists a predeclared late boundary layer after which

    sup_history delta_r(history)

is uniformly small across the operational preparation class.

If yes, the weighted theorem can be strengthened to a genuine uniform finite-window theorem.

If no, introduce one explicit auxiliary memory coordinate and test whether that produces the missing uniform closure.
