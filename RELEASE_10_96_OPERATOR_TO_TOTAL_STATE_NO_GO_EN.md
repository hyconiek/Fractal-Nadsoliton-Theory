# Release 10.96 — Operator-to-Total-State Nonuniqueness and the A-Only No-Go

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 10.96
- **Publication date:** 2026-08-30
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0

## Abstract

This release reports six analytical and computational research rounds,
ST1092–ST1181, testing whether the frozen FIN strict operator uniquely
determines the state space and complete dynamics of the one total nadsoliton.

The first result strengthens the diffusion interpretation. With

\[
A=sI-W,\qquad s=1.660307278766099,
\]

the generator \(Q=-A\) has positive off-diagonal rates and zero row sums.
Therefore

\[
P_t=e^{-tA}
\]

is an irreducible reversible mass-preserving Markov semigroup for every
\(t\ge0\). Its edge currents satisfy a continuity equation and Shannon entropy
production is nonnegative.

The central result is a total no-go for the current A-only state-selection
lane. The same operator canonically supports at least three pairwise
inequivalent normalized models:

\[
(\Delta_{11},e^{-tA}),
\qquad
(\mathbb{CP}^{11},e^{-itA}),
\qquad
(\mathcal D_{12},\rho\mapsto e^{-itA}\rho e^{itA}).
\]

Their real dimensions are 11, 22 and 143. Their boundaries, convex structures,
stationary sets, invertibility and long-time dynamics differ. The inherited
cyclic symmetry acts on every candidate and does not select one. Exact
12-to-24 refinements preserve both heat and unitary functional calculi but
repeat the ambiguity and leave a free fiber rate.

A fixed-time CPTP dilation of heat exists, but no single finite closed
environment with fixed initial state and time-independent Hamiltonian can
realize the convergent heat semigroup exactly for all times; finite closed
unitary dynamics is almost periodic. Strict heat is a detailed-balance
gradient flow with zero stationary currents, not a self-sustaining internal
pump.

Consequently, annihilation protection is model-relative: classical mass,
quantum norm/trace, or conditional total-state closure. The operator alone
does not determine which protection is physical. At least two additional axiom
classes are unavoidable: a state/observable category and a physical channel
with clock/instrument semantics. Repeating A-only state searches is now closed.

## Principal formulas

\[
Q=-A,\qquad Q_{ij}>0\ (i\ne j),\qquad Q\mathbf1=0,
\]

\[
\frac{dH}{dt}
=\frac12\sum_{i,j}W_{ij}(p_i-p_j)(\log p_i-\log p_j)\ge0,
\]

\[
A_{24}(q)=A_{12}\otimes I_2+I_{12}\otimes
\begin{pmatrix}q&-q\\-q&q\end{pmatrix},
\]

\[
Ce^{-tA_{24}(q)}=e^{-tA_{12}}C,
\qquad
Ce^{-itA_{24}(q)}=e^{-itA_{12}}C.
\]

## Included files

- English PDF and LaTeX source
- six executable research scripts
- six aggregate result JSON files and CSV summaries
- 90 individual JSON evidence packets
- combined regression/integrity test
- source bundle and SHA-256 manifest

## Scientific boundary

The no-go is total only for unique state/dynamics selection from the current
operator datum and its tested symmetry, functional calculus and refinement.
It does not exclude a richer FIN theory with independently derived algebraic,
nonlinear, refinement-limit or operational structure. No laboratory evidence,
physical scale, Standard Model, gravity or Theory-of-Everything closure is
claimed.
