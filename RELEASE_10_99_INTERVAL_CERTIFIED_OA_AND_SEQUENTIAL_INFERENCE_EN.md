# Release 10.99 — Interval-Certified OA Optimization, Composite Identifiability and Sequential Inference

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 10.99
- **Publication date:** 2026-08-31
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0

## Abstract

Release 10.99 reports six research rounds, ST1362–ST1451, upgrading the frozen
FIN Operational/State discrimination programme with outward interval
certificates, calibration-error budgets, full likelihoods and anytime-valid
sequential inference.

For an explicitly frozen seven-value decimal spectrum with multiplicities
\((1,2,2,2,2,2,1)\), interval arithmetic proves that the absolute return gap

\[
|Q(t)-C(t)|
\]

has one global maximum on \([0,10]\). Its location and value satisfy

\[
t_*\in[0.59453078,0.59453079],
\]

\[
|Q(t_*)-C(t_*)|
\in[0.4112409998076315,0.41124099980763185].
\]

A second complete interval cover proves four-time composite identifiability for
the declared energy-dephased quantum family with clock scale
\(\alpha\in[0.9,1.1]\) and dephasing \(\gamma\in[0,100]\):

\[
\mathrm{RSS}\ge0.046692011213846946.
\]

Therefore at least one of the four return probabilities differs from the
classical target by at least

\[
0.10804167160619894.
\]

The full twelve-vertex likelihood reduces without information loss to seven
cyclic-distance classes. Its Chernoff information is approximately 0.136480,
compared with 0.100286 for binary return. A conservative product Chernoff bound
gives 29 full-record shots for equal-prior error below one percent in the two
frozen simple hypotheses. The interval composite certificate yields a much
more conservative 1146 binary shots per time through a Hoeffding union bound.

For the simple seven-class models, the likelihood-ratio process is a
nonnegative martingale under the classical hypothesis, while its inverse is
one under the quantum hypothesis. Boundaries 100 and 0.01 therefore give
anytime-valid one-percent wrong-boundary control. This guarantee does not apply
to a nuisance likelihood refitted after each event.

Protocol 10.99 freezes the decimal spectrum, simple distributions, sequential
boundaries, composite nuisance box, certified lower bound, calibration budgets
and role separation. Synthetic fixtures test the validator only. No apparatus,
raw natural event, independent custody or physical state-category selection is
claimed.

## Principal results

\[
t_*\in[0.59453078,0.59453079],
\qquad
D_{\max}\in[0.4112409998076315,0.41124099980763185],
\]

\[
\inf_{\substack{0.9\le\alpha\le1.1\\0\le\gamma\le100}}
\sum_{t\in\{0.3,0.6,1.2,2.0\}}
(Q_{\alpha,\gamma}(t)-C(t))^2
\ge0.046692011213846946,
\]

\[
E_n=\prod_{i=1}^n\frac{q(X_i)}{p(X_i)},
\qquad
E_n\ge100\Rightarrow Q,
\qquad
E_n\le0.01\Rightarrow C.
\]

## Scientific boundary

The interval theorems are exact for the declared decimal spectrum and nuisance
box, not for symbolic kernel parameters or every open quantum channel. The
sequential e-process theorem is for two simple frozen models. No laboratory
evidence or strict FIN physics is exported.

## Included files

- English PDF and LaTeX source
- six research scripts and interval helper
- protocol 10.99, validator and sequential fixtures
- six result sets and 90 evidence packets
- certification figure, tests, source bundle and SHA-256 manifest
