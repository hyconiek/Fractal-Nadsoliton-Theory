# 316 — SIMPLEX-NATIVE-PREPARATION-CORRECTION
## Probability-preserving preparation correction, exactified N9 labels, and repair of the raw N9 pseudo-distribution

Date: 2026-09-27

Status: finite-N grouped cross-validation on N=3,...,9; frozen before any N=10 microscopic opening.

## Question

Task 314 found real low-dimensional preparation residual structure, but its additive Fourier correction could leave the probability simplex. Task 316 asks whether the residual can be used in a probability-preserving correction without materially degrading any N-group, while retaining report 313 as the independent process-level safety envelope.

Predeclared Pareto rule:

- every corrected distribution must be a genuine probability vector;
- no grouped leave-one-N-out max preparation TV may worsen by more than 0.001 absolute relative to the canonical simplex-projected baseline;
- among admissible models choose the smallest rank within 0.0005 TV of the best global worst grouped-LOO error.

## 1. N9 boundary classification was tightened before fitting

The previous N9 core construction left 38,096 threshold-ambiguous states. D12 symmetry reduces them to 1,640 orbit representatives.

Task 316 ran the full deterministic V_g descent on all 1,640 representatives and propagated labels by the exact D12 action.

After this pass:

- unresolved/nonlocalized states: 1,484;
- unresolved stationary mass: 0.000172948739;
- localized mass: 0.999827051261.

Thus the previous 0.1927% core ambiguity was reduced to a 0.0173% genuinely nonlocalized remainder for the preparation analysis.

## 2. Important repair: the raw N9 canonical prior was not in the simplex

The raw task-315 N9 Fourier reconstruction had

    min_j p_j = -0.006659734156.

So it was a pseudo-distribution, not a valid probability vector.

Task 316 therefore makes simplex sanitation mandatory:

    p_base^Delta = Proj_Delta(p_base_raw),

where Proj_Delta is the standard Euclidean projection onto the probability simplex.

For N=3,...,8 the existing baselines were already in the simplex and projection changes them only at numerical roundoff. For N9 the maximum projection displacement was

    TV(p_raw,p_projected) = 0.0133194683125.

This is a formal correction to the probabilistic interpretation of task 315, not a new post-hoc fit.

## 3. Strict logit/exponential correction fails

A natural candidate was

    p_corr,j proportional to (p_base,j + epsilon) exp([U z(m)]_j).

Ranks 1-3, several low-order m-laws, and fixed epsilon values from 1e-8 to 1e-3 were tested by grouped leave-one-N-out.

Result:

    zero Pareto-admissible candidates.

The best variants had worst LOO TV above 7.5%.

Reason: N8/N9 lie close to the simplex boundary. Log-ratio coordinates magnify small/zero shell probabilities and are not a stable cross-N coordinate system here.

This route is rejected.

## 4. Full projected additive correction also fails Pareto

The next probability-preserving family was

    p_corr = Proj_Delta(p_base^Delta + U_r z(m)).

A full rank-1/2/3 correction lowers the global worst case substantially, from the baseline 1.9231% to about 1.49% TV, but degrades some small-N groups by more than the predeclared +0.1 percentage-point tolerance.

Therefore the unshrunk correction is not promoted.

## 5. A single global shrinkage passes

The final declared family was

    p_corr = Proj_Delta(p_base^Delta + s [mu + u z(m)]),

with:

- one residual direction u;
- z(m) quadratic in the operational preparation mean m=E[n0/N];
- a single global shrinkage s on the fixed grid 0,0.05,...,1.

Grouped leave-one-N-out selects

    rank = 1,
    z = c0 + c1 m + c2 m^2,
    s = 0.20.

This is the simplest nonzero Pareto-admissible model and the best worst-case model among the selected parsimony class.

### Grouped LOO max TV

| N | baseline | 316 correction | change |
|---|---:|---:|---:|
| 3 | 0.4367% | 0.4662% | +0.0295 pp |
| 4 | 0.7472% | 0.8079% | +0.0607 pp |
| 5 | 0.2970% | 0.3922% | +0.0952 pp |
| 6 | 0.7417% | 0.6311% | -0.1105 pp |
| 7 | 1.9231% | 1.8386% | -0.0845 pp |
| 8 | 1.5099% | 1.4247% | -0.0852 pp |
| 9 | 1.5337% | 1.4462% | -0.0875 pp |

Hence

    boxed:
    B_316,prep = 0.0183859524612 TV.

The grouped mean TV is

    0.007781302670.

No group violates the +0.001 absolute Pareto tolerance.

## 6. Geometry after N9

In probability-residual coordinates, the first PCA direction carries

    69.9934%

of total residual variance, while the first two carry

    97.1884%.

Thus the preparation residual remains strongly low-dimensional, but the safely transferable component is weaker than the raw two-PC geometry suggested in task 314.

The shrinkage s=0.20 is therefore scientifically meaningful: only a small common component of the empirical residual is stable enough to transport across N without harming some groups.

## 7. Retrospective N9 probability repair

Using exactified N9 labels and the already frozen N9 transition:

- raw task-315 pseudo-prior minimum: -0.00665973;
- simplex-projected valid prior max TV vs microscopic N9: 1.53366%;
- selected 316 correction max prior TV: 1.42510%;
- raw pseudo-joint max TV: 2.48365%;
- simplex-projected valid joint max TV: 1.80473%;
- selected 316 corrected joint max TV: 1.72112%.

Therefore the task-315 holdout conclusion survives the repair and becomes cleaner probabilistically. The report-313 safety envelope remains much larger:

    B_313 = 5.582267% TV.

The 316 correction must not replace B_313 as a process certificate.

## 8. Verdict

PASS, but only as a modest optional central correction.

What passed:

1. mandatory simplex projection repairs nonphysical raw preparation predictions;
2. a rank-1 quadratic-in-m residual direction with global shrinkage s=0.20 passes the predeclared grouped-N Pareto rule;
3. it improves N=6,...,9 consistently and never worsens N=3,...,5 by more than 0.1 percentage point;
4. the correction is a probability vector by construction.

What did not pass:

1. full logit/exponential correction;
2. full-strength low-rank additive correction;
3. any claim that the residual law is an asymptotic N->infinity theorem.

Canonical policy going forward:

    raw cross-N prior
      -> mandatory simplex projection
      -> canonical valid baseline
      -> optional frozen rank-1 s=0.20 correction
      -> retain independent report-313 process envelope.

## 9. Next task

317 — FROZEN-N10-PREDICTOR-AND-DIRECT-HOLDOUT.

Before constructing any N10 microscopic state, freeze:

- the N10 clock and shape law;
- the preparation map;
- mandatory simplex sanitation;
- both central predictions: canonical projected baseline and optional 316 correction;
- report-313 process envelope as the independent pass/fail certificate.

Then open direct N10 without changing any of these objects.
