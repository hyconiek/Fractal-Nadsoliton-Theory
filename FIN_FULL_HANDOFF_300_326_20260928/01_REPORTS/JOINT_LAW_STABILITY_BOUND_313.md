# JOINT-LAW-STABILITY-BOUND-313
## Exact joint/path stability theorem and an auditable preparation + transition + reduction error budget

Date: 2026-09-27

Status:
- exact finite-alphabet total-variation theorems;
- exact two-time decompositions for grouped holdouts N=3..6;
- retrospective validation on already-opened N=7 and direct microscopic N=8;
- an empirical, pre-frozen grouped-N triangle envelope for the next untouched holdout;
- no distribution-free cross-N theorem is claimed.

## 1. Research question

Reports 311-312 showed a sharp asymmetry:

- strict posterior sets can become extremely wide after conditioning on a rare first readout;
- the full observed two-time joint law remains much more stable and generalizes well to N=7 and direct N=8.

The goal of 313 is to replace posterior-based certification by a theorem with three separately auditable sources of error:

1. microscopic -> same-N effective reduction;
2. same-N -> cross-N predicted transition;
3. actual -> predicted preparation state.

The main object is the full observed joint/path law, not a posterior conditioned on one rare outcome.

## 2. Two-time setup

Let `i` be a latent 12-state localized label after burn-in.

Let

    p_i

be the actual latent prior, and

    p_hat_i

its cross-N prediction.

Let `E(a|i)` be the fixed first-readout channel.

Let

    H_N(b|i)

be the same-N effective downstream observed channel and

    H_hat(b|i)

its cross-N prediction.

Define

    J(p,H)_(a,b)
      = sum_i p_i E(a|i) H(b|i).

Let `M(a,b)` be the exact microscopic two-time observed law.

Introduce the intermediate laws

    S = J(p,H_N),
    C = J(p,H_hat),
    J_hat = J(p_hat,H_hat).

Thus

    microscopic M
      -> same-N S
      -> cross-N transition C
      -> full prediction J_hat.

## 3. Theorem 313-A — exact auditable two-time decomposition

By the triangle inequality,

    TV(M,J_hat)
      <=
    TV(M,S)
      +
    TV(S,C)
      +
    TV(C,J_hat).

Define

    eps_red   = TV(M,S),
    eps_trans = TV(S,C),
    eps_prep  = TV(C,J_hat).

Then exactly

    boxed:
    TV(M,J_hat)
      <=
    eps_red + eps_trans + eps_prep.

No Markov assumption is made about the microscopic projected process.

The three terms have distinct meanings:

- `eps_red`: same-N microscopic/history/reduction error;
- `eps_trans`: error of the cross-N transition law with preparation held fixed;
- `eps_prep`: error caused by the preparation/prior map with the predicted transition held fixed.

This is the decomposition used below.

## 4. Theorem 313-B — channel bounds for transition and preparation

For fixed `p`,

    TV[J(p,H_N),J(p,H_hat)]
      <=
    sum_i p_i TV[H_N(.|i),H_hat(.|i)]
      <=
    sup_i TV[H_N(.|i),H_hat(.|i)].

For the preparation term define the pair-output channel

    K_hat_i(a,b)=E(a|i) H_hat(b|i).

Its Dobrushin coefficient is

    eta(K_hat)
      = max_(i,j) TV(K_hat_i,K_hat_j)
      <= 1.

Then

    TV[J(p,H_hat),J(p_hat,H_hat)]
      <=
    eta(K_hat) TV(p,p_hat)
      <=
    TV(p,p_hat).

Crucially, neither bound contains a factor

    1 / P(Y1=a).

Therefore rare first-readout outcomes cannot blow up the joint-law error in the way they blow up a conditional posterior radius.

## 5. Theorem 313-C — finite-window extension

Let a same-N effective latent Markov model have initial prior `p0` and transition kernels

    T_1,...,T_(m-1),

and let the cross-N prediction use

    p_hat0,
    T_hat_1,...,T_hat_(m-1).

After applying arbitrary observation channels at every time, data processing and the standard telescoping-product argument give

    boxed:
    TV(P_sameN^Y, P_pred^Y)
      <=
    TV(p0,p_hat0)
      +
    sum_r sup_i TV[T_r(i,.),T_hat_r(i,.)].

If `P_micro^Y` is the true microscopic observed path law, then

    boxed:
    TV(P_micro^Y,P_pred^Y)
      <=
    eps_red_path
      +
    TV(p0,p_hat0)
      +
    sum_r eps_trans,r,

where

    eps_red_path = TV(P_micro^Y,P_sameN^Y).

Report 309 supplies an independent exact composition theorem for this first term:

    eps_red_path
      <=
    initial marginal reduction error
      +
    sum_r E_micro[delta_sameN,r(history)].

Combining 309 and 313 therefore gives a full finite-window architecture:

    preparation error
      + cross-N transition errors
      + microscopic history/reduction defects
      -> observed finite-window path error.

This is an exact theorem once the individual terms are bounded.

## 6. Grouped-N calibration: N=3..6

The two-time decomposition was replayed for leave-one-N-out predictions on N=3,4,5,6.

| N | reduction | transition | preparation | triangle sum | actual joint TV |
|---:|---:|---:|---:|---:|---:|
| 3 | 4.18938% | 0.93502% | 0.45786% | 5.58227% | 5.10999% |
| 4 | 1.82547% | 0.38799% | 0.69035% | 2.90289% | 1.86779% |
| 5 | 2.08888% | 0.18046% | 0.20183% | 2.47118% | 1.92710% |
| 6 | 1.04354% | 0.69409% | 1.52963% | 3.26726% | 1.95136% |

Every actual joint-law error lies below the corresponding explicit triangle sum, as required.

### Two possible grouped envelopes

Taking the maximum of every component independently gives the very conservative

    B_red   = 4.18938%,
    B_trans = 0.93502%,
    B_prep  = 1.52963%,

hence

    B_union = 6.65403% TV.

This is valid if each component is only known through a separate grouped maximum.

A less wasteful but still auditable empirical envelope keeps the three terms paired within the same grouped holdout and takes the maximum only after summation:

    boxed:
    B_313 = 5.582267457% TV.

This is the primary empirical envelope frozen for the next untouched holdout.

For comparison, the direct final-joint training envelope is

    5.109993% TV,

but it does not expose which physical/modeling component generated the error.

## 7. Retrospective validation on already-opened N=7 and N=8

These are NOT new blind holdouts in task 313. They were opened in earlier tasks and are used only to test whether the decomposition behaves sensibly.

### N=7

    reduction   = 0.75908%,
    transition  = 0.51445%,
    preparation = 1.94203%,
    triangle    = 3.21556%,
    actual joint= 2.10890%.

### N=8

    reduction   = 0.38336%,
    transition  = 0.25271%,
    preparation = 1.42957%,
    triangle    = 2.06564%,
    actual joint= 1.74418%.

Both are well below the grouped N=3..6 triangle envelope

    5.58227%.

The theorem is therefore numerically non-pathological across N=3..8.

## 8. Error-source crossover

The composition of the error budget changes with N.

At N=3 the reduction/history term contributes about

    75.0%

of the triangle sum.

At N=7 the preparation term contributes

    60.4%,

and at N=8

    69.2%.

Thus the present high-N bottleneck is not microscopic memory or transition dynamics.

It is the cross-N preparation/prior map.

This is consistent with the independent diagnoses of reports 311 and 312.

## 9. Comparator implication

Using the empirical grouped triangle envelope only as a retrospective diagnostic:

- N=7 predicted FIN-vs-comparator joint separation was 6.23376%;
- N=8 predicted separation was 6.69236%.

Subtracting `B_313` gives positive retrospective margins:

    N=7: 0.65149 percentage points,
    N=8: 1.11009 percentage points.

This should NOT be relabeled as a blind certification, because both N=7 and N=8 had already been opened before report 313.

The scientific use is prospective: `B_313` is now frozen before any direct microscopic N=9 opening.

## 10. What is now frozen for a future N=9 test

Primary auditable empirical threshold:

    boxed:
    B_313 = 0.05582267457195117 TV.

A future pre-open prediction can be called empirically pre-certified against a declared comparator only if its predicted observed joint-law separation exceeds this number.

After opening the holdout, the following must be reported separately:

1. microscopic-to-same-N reduction term;
2. same-N-to-cross-N transition term;
3. preparation/prior-map term;
4. their triangle sum;
5. the actual full joint-law error.

The threshold must not be widened after observing N=9 and still be called pre-certified.

## 11. Epistemic boundary

The inequalities in sections 3-5 are exact probability theorems.

The numerical envelope `B_313` is not a theorem in N. It is an empirical grouped-N calibration from the finite set N=3..6, retrospectively consistent with N=7 and N=8.

Nothing here derives:

- a fundamental physical carrier;
- the source of A7 or g;
- physical space or dimensional units;
- QM or GR;
- a distribution-free large-N guarantee.

## Verdict

P0-A 313 is PASS.

The observed process now has a clean and auditable stability architecture:

    microscopic reduction
      + transition transfer
      + preparation transfer
      -> finite-window observed-law error.

The main practical bottleneck at N=7,8 is preparation transfer.

The next P0 should therefore not add memory to the generator. It should model the strongly low-dimensional preparation residual and freeze that correction before the next direct microscopic holdout.
