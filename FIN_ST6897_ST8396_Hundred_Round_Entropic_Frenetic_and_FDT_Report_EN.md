# FIN ST6897–ST8396 — Hundred-Round Entropic–Frenetic and FDT Report

## Executive theorem

Every positive pair of transition rates has the unique decomposition

\[
\boxed{
k_{xy}=a_{xy}e^{F_{xy}/2},
\qquad
k_{yx}=a_{xy}e^{-F_{xy}/2},}

\]

where \(F_{xy}=-F_{yx}\) is the time-antisymmetric affinity and
\(a_{xy}=a_{yx}>0\) is the time-symmetric dynamical activity.

Local detailed balance determines the ratio

\[
\frac{k_{xy}}{k_{yx}}=e^{F_{xy}},

\]

but does not determine \(a_{xy}\). This is the exact entropic–frenetic split.

## Equilibrium no-section result

For a Gibbs state \(\pi\), equilibrium detailed balance fixes

\[
F_{xy}=\log\frac{\pi_y}{\pi_x}.

\]

Every positive symmetric conductance/activity still defines a different
reversible generator with the same \(\pi\). Its spectral gap, mixing time,
waiting-time law and response functions vary.

Therefore:

\[
\boxed{
\text{there is no natural map from an equilibrium state alone to one kinetics}.}

\]

This is the kinetic version of the realization-fibration no-section theorem.

## Path-space decomposition

The logarithmic ratio between a path and its time reversal determines entropy
flux. The time-symmetric part of path action contains escape rates, waiting
times and dynamical traffic. Fluctuation theorems constrain the former but not
the latter.

Two processes may therefore have:

- the same equilibrium state;
- the same local affinities;
- the same entropy-flux ratios;

while having different event frequencies, relaxation spectra and response.

## Nonequilibrium currents

Away from equilibrium, stationary flow decomposes into symmetric traffic and
antisymmetric current. Entropy production has the schematic form

\[
\sigma=\sum_{x<y}J_{xy}F_{xy}.

\]

It does not reconstruct symmetric traffic. Reversing current polarity can
preserve the stationary state and entropy-production magnitude. A physical
arrow needs a driven resource, orientation and clock.

## Fluctuation–dissipation boundary

In continuous form, FDT relates noise covariance and mobility, for example

\[
D=\mu k_BT.

\]

This relation reduces independent parameters but does not determine mobility
\(\mu\), temperature calibration, time unit or microscopic dynamics. Onsager
reciprocity constrains a transport matrix without fixing its entries. Green–
Kubo relations require the underlying time-correlation law that they are meant
to evaluate.

Thus FDT is a genuine cross-fibre theorem, but not a complete section.

## Information-geometric boundary

Relative entropy or free energy becomes a dynamics only after choosing an
Onsager operator, Wasserstein metric, Fisher/natural-gradient preconditioner or
reference process. Maximum caliber likewise needs a reference path law and
dynamical constraints.

The information functional supplies a landscape. The frenetic structure
supplies how the system moves on it.

## Application to FIN

FIN currently supplies several pieces of the entropic sector:

- a strict finite Dirichlet operator;
- equilibrium and ternary-Gibbs candidate shapes;
- exact prism equilibrium RG;
- Shannon and relative-entropy functionals;
- finite dual unitary and diffusive channels.

It does not source the complete frenetic sector:

- edge activities or mobilities;
- configuration-space clock scale;
- refinement-rate sequence;
- nonreversible current polarity;
- environment and bath;
- calibrated event record.

The common generator \(A\) organizes unitary and heat dynamics but does not
identify them with a unique stochastic configuration process.

## Hierarchical consequence

Coarse-graining a prism or Cantor refinement changes hidden traffic and may
create memory even when equilibrium marginals close exactly. A stationary RG
identity therefore does not imply dynamic RG closure. Each level needs an
activity/clock law or a non-Markovian memory kernel.

## Operational test

Affinity and activity require different records:

- forward/backward transition ratios estimate affinity;
- waiting times and total traffic estimate activity;
- signed cycle counts estimate nonequilibrium current;
- calibrated timestamps determine clock scale;
- temperature and energy calibration are separate.

A state histogram alone cannot recover these objects.

## Master gate

| Requirement | Result |
|---|---|
| rate-factorization theorem | PASS |
| entropic–frenetic decomposition | PASS |
| equilibrium kinetic no-section | PASS |
| FDT residual-mobility theorem | PASS |
| path-space factorization | PASS |
| operational record separation | PASS |
| strict activity/mobility source | FAIL |
| strict clock and temperature | FAIL |
| OA evidence | FAIL |
| SM/GR/\(L_{\rm total}\)/ToE closure | FAIL |

\[
\boxed{6/10\ \text{master mathematical rows},qquad
0/10\ \text{physical-closure rows}.}
\]

## Deepest conclusion

Equilibrium information determines thermodynamic affinities but not dynamical
traffic. FIN's missing dynamics is not merely a coupling constant: it is the
time-symmetric frenetic geometry plus clock and, out of equilibrium, a signed
current source.

The highest-value next result must derive a strict mobility/activity operator
from independent FIN structure and demonstrate that it lowers the realization
rank. Without it, FDT remains a conditional bridge rather than a physical
completion.

No physical temperature, time unit, laboratory evidence, Standard Model,
gravity, \(L_{\rm total}\), or Theory of Everything closure follows.
