# FIN ST2442–ST2531 — Strict-Shape Ternary Gibbs Law and Prism Renormalization

## Executive result

The strict loop tensor

\[
\tau_{ijk}=W_{ij}W_{jk}W_{ki}
\]

defines the strongest currently available non-Gaussian state candidate:

\[
p_\theta(x)=Z(\theta)^{-1}
\exp\!\left(\theta\sum_{i<j<k}\tau_{ijk}x_ix_jx_k\right),
\qquad x_i\in\{-1,+1\}.
\]

The interaction shape is strict-derived. The coupling \(\theta\), its sign,
the probability-law postulate, preparation and physical realization are not.

The most important new positive theorem is that the corresponding three-bit
family is exactly closed under symmetric triangular-prism marginalization.
The most important negative theorem is that closure introduces an additional
vertical coupling \(q\), so \(\theta\) is not the single missing scalar.

## Round I — existing nonlinear laws

The repository candidates were classified as follows:

| Candidate | Third-order source status |
|---|---|
| focusing/saturating DNLS | cubic equation but supplied fourth-order action coefficient; no strict three-site statistic |
| negative-information Landau model | active gain \(g\) is supplied; first strict angular anisotropy is order 12 |
| Hebb/Oja | covariance/second-order learning only |
| BCM/STDP | nonlinear/directional possibilities, but thresholds or delays are additional laws |
| operator potential \(V(K)\) | arbitrary until its higher jet is sourced |
| loop-tensor Gibbs law | strict interaction shape; free \(\theta\) |

For all 4096 binary configurations the conditional Gibbs law is finite and
positive. At \(\theta=0.1\), all 220 triangle moments are nonzero. At small
coupling,

\[
\frac{d}{d\theta}\langle x_ix_jx_k\rangle_{\theta=0}
=\tau_{ijk}.
\]

The numerical residual at \(\theta=10^{-3}\) is
\(1.20\times10^{-10}\).

## Round II — symmetry and polarity

The tensor is \(D_{12}\)-covariant and spatial-reflection even. A global bit
flip gives

\[
x\mapsto-x,
\qquad
S(x)\mapsto-S(x),
\qquad
\theta\mapsto-\theta.
\]

Therefore

\[
Z(\theta)=Z(-\theta),
\qquad
\langle S\rangle_{-\theta}=-\langle S\rangle_\theta.
\]

Entropy and Fisher information are even in \(\theta\). No scalar strict datum
currently selects one polarity. A finite 4096-state partition function is
analytic, so genuine spontaneous phase selection cannot occur at finite
volume without an additional limit or dynamical law.

## Round III — stationary-family theorem

The model is a regular one-parameter exponential family. Its log partition is
strictly convex:

\[
\frac{d^2}{d\theta^2}\log Z=\operatorname{Var}_\theta(S)>0.
\]

Thus a supplied value of \(\langle S\rangle\) identifies a unique \(\theta\).
This is identifiability, not provenance. Rescaling \(\tau\) can be absorbed
into \(\theta\), and a thermodynamic interpretation \(\theta=\beta J\)
requires an energy coupling and temperature.

## Round IV — exact prism renormalization closure

For bottom spins \(b_i\), top spins \(t_i\), triangle coupling \(\theta\) and
vertical pair coupling \(q\), use

\[
H=\theta(b_1b_2b_3+t_1t_2t_3)+q\sum_{i=1}^3b_it_i.
\]

After summing over the top layer, every one-body and pair-body Walsh
coefficient cancels. The measured maximum pair coefficient is
\(5.55\times10^{-17}\). The marginal remains a pure three-body law.

Let

\[
A(q)=e^{3q}+3e^{-q},
\qquad
B(q)=3e^q+e^{-3q}.
\]

Then the exact renormalized coupling is

\[
\boxed{
\theta_{\mathrm{eff}}
=\theta+rac12\log
\frac{e^\theta A(q)+e^{-\theta}B(q)}
     {e^{-\theta}A(q)+e^\theta B(q)}.}
\]

For \(\theta=0.3,q=0.2\),

\[
\theta_{\mathrm{eff}}=0.3022399512016273.
\]

This is a genuine exact renormalization theorem. It does not close FIN because
each level still needs \(q_n\), or a strict law determining the sequence.

## Round V — dual dynamics and observability

The strict unitary and heat channels share the generator \(A\), but neither is
automatically a Markov dynamics on the 4096 binary configurations. A Glauber
or other detailed-balance generator can realize \(p_\theta\), but its rates,
clock, bath and coupling are additional objects.

A joint parity instrument measures the three-body signal. Pair-only records
cannot reconstruct the fine coupling or its vertical refinement parameter.
An operational record must retain at least run ID, three bit outcomes, face,
layer, time and configuration.

## Round VI — gate

| Requirement | Result |
|---|---|
| strict loop tensor \(\tau\) | PASS |
| finite non-Gaussian family | PASS as mathematics |
| nonzero triangle response | PASS |
| strict \(\theta\) magnitude | FAIL |
| strict \(\theta\) polarity | FAIL |
| unique refinement \(q_n\) law | FAIL |
| strict dual-channel realization | FAIL |
| physical units | FAIL |
| OA platform and record | FAIL |
| laboratory evidence | FAIL |

Final score:

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

## Deepest interpretation

The strict pair kernel determines a natural shape for a non-Gaussian
three-body interaction. The triangular prism supplies an exact closed
renormalization map for that interaction. These are substantial mathematical
structures. They do not determine the statistical law, coupling polarity,
vertical refinement sequence, configuration-space dynamics, units or physical
realization.

The hypothesis that only one scalar \(\theta\) was missing is refuted: prism
refinement introduces at least the vertical coupling sequence \(q_n\).

## Highest-value next theorem

Derive or refute a strict configuration-space Markov generator whose
detailed-balance affinity simultaneously:

1. produces the loop-tensor Gibbs law,
2. fixes the sign and magnitude of \(\theta\),
3. supplies the vertical sequence \(q_n\),
4. intertwines under \(12\to24\to48\),
5. and preserves the exact prism renormalization formula.

Without such a generator, the candidate remains a conditional statistical
model rather than fundamental physics.
