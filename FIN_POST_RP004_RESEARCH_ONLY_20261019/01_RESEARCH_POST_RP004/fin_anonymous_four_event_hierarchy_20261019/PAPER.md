# FIN mathematical research: first discriminating anonymous four-event histories and a finite moment hierarchy on Z12

**Date:** 2026-10-19. **Type:** class-conditional mathematical theorem and independent exact replay. **Not:** sourced intrinsic nadsoliton dynamics, native q12 N2 interacting heat-bath law, a physical time law, experimental evidence, or an externally refereed paper.

## Abstract

The inherited comparison model provides two strictly positive, non-isomorphic, reversible and dihedrally symmetric continuous-time jump generators on the twelve-state carrier. These generators are isospectral and consequently yield identical stationary two-time return probabilities at every model-time separation and identical return probabilities after any specified number of actual jumps. We show that the **complete anonymous equality-pattern record of four consecutive actual jump states nevertheless distinguishes the two laws, and that no anonymous equality record of at most three consecutive jump states does so**. The event `ABAB` has exact difference `11/(24*113^3)`; the total variation distance between full four-index equality-partition distributions is `11/(12*113^3)`. This supplies a precise, model-specific hierarchy: two-time scalar return data fail; three-index equality patterns still fail; a four-index time-ordered relational correlation separates the frozen pair, with no clock calibration. A second theorem shows that all alternating equality records of lengths at most twelve recover the **unordered eleven-element multiset** of relative outgoing jump probabilities (through Newton's identities), while leaving circular shell assignment, overall activity and intrinsic source law unlicensed. We supply a separate, rational counterpair proving that alternating records of *every* length do not determine cyclic shell assignment. An optional continuous-time corollary gives a positive cubic small-lag signature, but assumes the common model clock. All model-generating equations remain premises.

## 1. Frozen sources, carrier and method

The actual local attached input is `fin_relational_event_history_20261019/PAPER.md` and `verify.py`, with its `AGENTS(1).md` and `GEN001-gen956(1).md` snapshots. The original source explicitly classifies the 12-state one-label radial `Q=-A` jump process as a **comparison model**, not the native many-copy N2 LOO heat bath. The exact pair below is copied as frozen mathematics and all identities are freshly replayed. The AGENTS record preserves `OPEN_SOURCE` for a natural kinetic provider and for geometry and units. Existing six-mode labelled reconstruction and its normalization are not rerun as a new discovery.

Let `Z12={0,...,11}` be an *assumed finite label set*. Fix six positive shell weights `w_d>0` with radial distance `d(i,j)=min(|i-j|,12-|i-j|)`. The 12×12 generator has `Q_{ij}=w_{d(i,j)}` for `i!=j` and `Q_{ii}=-s`, where `s=2 sum_{d=1}^5 w_d+w_6`. This is a reversible continuous-time Markov law, supplied rather than derived from FIN ontology. Its embedded chain of **actual, non-self jumps** is `J=I+Q/s`, with stationary uniform distribution `pi_i=1/12`, `J_{ii}=0`, and `J_{ij}=w_{d(i,j)}/s` for `i!=j`. We observe only equality relations among the *ordered* states `Y_0,Y_1,...` of this chain, not their labels, distances, Fourier directions, or duration between actual jumps.

The inherited exact counterpair is

```
w = (16/15+sqrt(3)/30, 1, 61/60, 1, 16/15-sqrt(3)/30, 1)
v = (43/40+sqrt(3)/40, 59/60, 1, 1, 43/40-sqrt(3)/40, 31/30)
s = 113/10.
```

Both `J_w` and `J_v` are positive away from their diagonal, and their continuous-time generators differ by a swap of the two doubly-degenerate eigenvalues `lambda_1,lambda_2`. Thus `tr J_w^n=tr J_v^n` for every positive integer n, and their stationary two-time returns coincide at all model times. The earlier paper proved they are not conjugate by any permutation: their entrywise off-diagonal cubic moments differ by `11/24000` before dividing by `s^3`.

## 2. Theorem E3 — the first separating anonymous event record

**Theorem E3.** In the frozen exact counterpair, the equality-pattern distribution of any block of at most three consecutive actual jump states is identical. For four consecutive states, the event

`ABAB := {Y_0=Y_2 and Y_1=Y_3}`

satisfies

```
Pr_w(ABAB)=431959/51944292,
Pr_v(ABAB)=7645/919368,
Pr_w(ABAB)-Pr_v(ABAB)=11/(24*113^3)=11/34629528 > 0.
```

The **complete** anonymous four-index partition distribution has

```
TV(P_w^eq, P_v^eq)=11/(12*113^3)=11/17314764 ≈ 6.3529598209e-7.
```

Thus four is the minimum number of consecutive event-indexed states needed to distinguish **this pair** using only state-equality patterns. This is a minimum within this stated observation contract, not a universal bound on all sources or observables.

**Proof.** Every pair of consecutive actual jump states is unequal, because `J_ii=0`. For three consecutive states the only nontrivial anonymous equality event is `Y_0=Y_2`, whose stationary probability is `(1/12)tr(J^2)`. The two chains are isospectral, hence their entire equality-pattern laws up to length three coincide.

For four consecutive states, define `a_r=(1/12)tr(J^r)` for r=2,3 and `b=Pr(ABAB)`. Translation invariance and reversibility yield the exact alternating probability

```
b = (1/12)*sum_{i,j} J_ij J_ji J_ij
  = sum_{j != 0} J_0j^3
  = [2*sum_{d=1}^5 w_d^3+w_6^3]/s^3.
```

At four event times the only five possible equality partition types, and their probabilities, are:

| Pattern | Probability |
|---|---|
| ABAB | b |
| ABAC | a_2-b |
| ABCB | a_2-b |
| ABCA | a_3 |
| ABCD | 1-2 a_2-a_3+b |

Indeed `ABCA` corresponds to `Y_0=Y_3`; overlapping relations `Y_0=Y_2` and `Y_1=Y_3` coincide only in `ABAB`; any other pair of nonadjacent coincidences would imply a forbidden adjacent self-jump. The two chains have identical a2 and a3 by isospectrality, whereas the fixed cubic off-diagonal moment gives

```
b_w-b_v=[11/24000]/(113/10)^3 = 11/(24*113^3).
```

The four of five pattern probabilities that depend on b differ by `+Δb,-Δb,-Δb,+Δb`; `ABCA` agrees. The total variation distance, one half of the sum of absolute category differences, equals `2|Δb|`. This also proves minimality relative to ≤3 event states. □

**Operational interpretation.** If one independently prepared a sample of four states from either of these two *mathematical* stationary laws, and both hypotheses had equal prior weights, the maximum theoretical one-shot classification probability from the equality partition would be `1/2+|Δb| ≈ 0.500000317648`. This is tiny and is not an experimental claim. Neither access to an intrinsic observer nor independent repetitions is demonstrated by FIN.

## 3. Theorem E4 — alternating-history moment tomography and its boundary

Define the anonymous length-`r+1` alternating equality event `A_{r+1}` by

```
Y_0=Y_2=Y_4=... and Y_1=Y_3=Y_5=...,
```

with two different (unnamed) alternating states, automatically enforced by `J_ii=0`. For all r≥1 its exact probability is

```
p_r := Pr(A_{r+1})
     = sum_{j != 0} J_0j^r
     = 2*sum_{d=1}^5 (w_d/s)^r+(w_6/s)^r.
```

Here `p_1=1` follows from row stochasticity, and `p_2,...,p_11` are observable probabilities of alternating equality records of lengths 3 through 12.

**Theorem E4 (unordered jump-edge recovery).** Knowledge of `p_1,...,p_11` uniquely determines the multiset

```
{w1/s,w1/s,w2/s,w2/s,w3/s,w3/s,w4/s,w4/s,w5/s,w5/s,w6/s}
```

with multiplicities, even if some relative shell weights coincide. In particular, `w6/s` is the unique multiset value occurring with odd multiplicity. The remaining five relative shell weights are recovered as a multiset, **without their assignment to shell labels 1,...,5**. The common scale `s` is not recovered.

**Proof.** Write the eleven outgoing probabilities of any row of J as `x_1,...,x_11`. The anonymous alternating events deliver the power sums `p_k=sum_j x_j^k` for k=1,...,11. Newton's identities recursively calculate elementary symmetric polynomials:

```
e_0=1,
e_k=(1/k) sum_{j=1}^k (-1)^(j-1)*e_{k-j}*p_j, k=1,...,11.
```

Consequently all coefficients of the monic polynomial

```
P(z)=z^11-e1*z^10+e2*z^9-...-e11
    =(z-w6/s) prod_{d=1}^5 (z-w_d/s)^2
```

are determined. Its roots with multiplicities uniquely determine the unordered eleven-element set. Because five shell probabilities each contribute *two* copies and the antipodal shell contributes *one*, exactly one **distinct numerical value** has an odd total multiplicity; it must be w6/s. Subtract one occurrence of that value to recover the other five as an unordered multiset. Multiplying all unnormalised w by a common positive c does not change J or any event-indexed history. □

**Sharp limitation.** Consider two strictly positive rational shell vectors `q=(1,2,3,4,5,6)` and `q'=(2,1,3,4,5,6)`. Their total exit rate is 36 and they have **identical** `p_r` for every r≥1 (not merely r≤11), because they have the same outgoing-rate multiset. Nevertheless their anonymous three-step return probability is different:

```
(1/12)*tr(J_q^3)-(1/12)*tr(J_q'^3)=-11/7776 != 0.
```

Therefore the entire alternating-record moment hierarchy still does **not** identify the circular assignment of the weights. Other equality-history shapes contain additional incidence information. This independent control also prevents mistaking unordered jump strengths for derived spatial geometry.

## 4. Continuous-time four-point corollary (separate observation contract)

If a **common model time** is supplied, consider actual CTMC states at times `0,t,t+u,t+u+v` and require the first and third to match and the second and fourth to match. Writing `K_t(d)=(e^{tQ})_{0d}`, reversibility and translation invariance give

```
G_Q(t,u,v) = sum_{d=0}^{11} K_t(d)*K_u(d)*K_v(d).
```

For two isospectral circulant generators the d=0 term agrees exactly because `K_t(0)=tr(e^{tQ})/12` is fixed by the spectrum. For d≠0, `K_t(d)=t*w_{d(0,d)}+O(t^2)` at t↓0. Hence for the frozen witness pair and equal small lags τ:

```
G_w(τ,τ,τ)-G_v(τ,τ,τ)
  = (11/24000)*τ^3 + O(τ^4).
```

This offers a second, time-calibrated version of the discriminating four-event motif. It does **not** eliminate the need for an independent clock; the jump-indexed Theorem E3 avoids that additional premise.

## 5. FIN implications and class/source audit

| Obligation | What is actually paid | What is NOT paid |
|---|---|---|
| Anonymous histories | Four-event equality distinguishes the exact frozen previously isospectral counterpair; 3-event equality does not | All laws identifiable from four events |
| Finite moment hierarchy | Alternation words through twelve events reconstruct unordered normalised jump strengths exactly | Cyclic shell assignment; physical geometry |
| Alternative histories | A nonalternating return word can detect circular placement missing from all alternating moments | Unique generic intrinsic evolution |
| Relative time | Jump count is an ordered internal event index | Natural waiting-time law, physical seconds or time arrow |
| FIN F01 / F03 | Concrete observation restrictions inside an already chosen generator class | Any **independent** activity, kinetic or nonlinear coefficient source |
| FIN F06 | Explicit carrier-labelled limitations | Spacetime, interactions of multiple physical sectors, strict-to-field bridge |

The 12 labels are **states**, not physical positions, and the anonymous relation `same state` is an **assumed observer/readout**. Our Markov premises already supply the rate law. Even if all event-index histories identified its generator up to relabelling, that would reconstruct an assumed law from its outcomes, not establish why primordial information must obey that law. No `AGENTS.md` or existing FIN source is changed.

## 6. Reproducibility and next conditional target

`python verify.py` separately reconstructs both jump matrices, checks isospectrality, computes exact cubic moments and partition distributions, enumerates all 12^4 possible trajectories as a floating independent control, verifies the exact Newton-polynomial reconstruction for rational and radical data, and checks an alternative isomoment control. It does not purport to mechanically check all theorem prose.

**Recommended next step, only upon a real provider:** If FIN supplies independently (not from a preselected Q) any joint history observable, characterize its allowed *equality types* and composition. Test whether it separates the E2 witness using `ABAB` and whether it can determine the shell assignment or activity. If FIN only supplies static/trace/return information, keep `OPEN_SOURCE` and the demonstrated class obstruction. Do not add a bath, Poisson clock, detector, ring geometry, or postulate merely to force a positive result.

**Disposition:** `CONDITIONAL_FOUR_EVENT_RELATIONAL_SEPARATION_AND_UNORDERED_MOMENT_TOMOGRAPHY; F01_F03_F06_OPEN_SOURCE`.
