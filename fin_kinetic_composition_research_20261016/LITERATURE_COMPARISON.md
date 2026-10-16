# Primary-source comparison and duplicate audit

Date: 2026-10-16. User requested browser-harness. Public arXiv HTML was
opened in the real browser through that skill; saved DOM-text captures and
selected equation-preserving sections are in `evidence/`. Search summaries,
including Google's generated overviews, were discovery aids only.

This is a targeted six-source comparison, not an exhaustive literature review
or a claim that all papers published by 2026 were searched. Versions below
are the versions actually consulted. HTML conversion sometimes prints a
2026 date inside an older manuscript; the arXiv version/publication metadata,
not that rendering date, determines the bibliographic year.

## L01 — Hamiltonian source with a specified compression

Alexander Mielke, Mark A. Peletier, Johannes Zimmer,
[Deriving a GENERIC system from a Hamiltonian system](https://arxiv.org/html/2404.09284v2),
2024, v2. Read §1.2–1.6, Assumption 1.2, Lemma 3.2, Theorem 3.3,
Theorem 4.10 and Appendix C, especially Remark C.2.

The microscopic inputs include a Hamiltonian/Poisson system, a linear
infinite-dimensional bath, coupling and projection. The compression identity
and Gaussian preparation produce an OU process; an additional energy variable
has a microscopic interpretation. Positive-temperature bath states require
an enlarged space and have infinite bath energy in this idealization.

**FIN assessment:** this is the strongest *source-construction template*
found here, but its inputs do not come from FIN's Gibbs measure. A designed
dilation of a prescribed target D reverses the explanatory direction. The
paper is not a theorem about FIN's finite-state labelled heat-bath path law.
The regularity test in THEORETICAL_RESULTS §3 is an elementary quantitative
application of its domain issue, not a contradiction of the paper.

## L02 — Generalized fluctuation symmetry and a reversible current

Richard Kraaij, Alexandre Lazarescu, Christian Maes, Mark Peletier,
[Deriving GENERIC from a generalized fluctuation symmetry](https://arxiv.org/html/1706.10115v1),
2017, v1. Read §II.2 and §III.2, equations (22)–(27), (69)–(71), and conclusions.

A shifted fluctuation symmetry identifies a nondissipative current in addition
to gradient relaxation. The paper separates pre-GENERIC from full GENERIC and
discusses auxiliary extensions. Its inertial examples contain momentum or
velocity variables and an appropriate reversal.

**FIN assessment:** applying the shifted-symmetry criterion to the admitted
TH occupation action yields W=0 when its action and entropy are retained.
This is the concrete new comparison performed here. It does not prohibit
hidden currents or extended carriers, and the general published method is
not claimed as a FIN discovery.

## L03 — Gradient structure from a supplied path action

Alexander Mielke, Mark A. Peletier, D. R. Michiel Renger,
[On the relation between gradient flows and the large-deviation principle,
with applications to Markov chains and diffusion](https://arxiv.org/html/1312.7591v1),
2013 preprint, published 2014. Read Theorems 2.3 and 3.3, §4.1.

The uniqueness statement is conditional on the given L-function and the
specified gradient representation. The finite-chain example starts from
a Markov generator; it does not recover a microscopic generator from static
entropy alone. The entropy factors and current conventions must be kept.

**FIN assessment:** this is already in the repository's theory synthesis,
TH and `fin_jump_fluctuations`. The current comparison reuses TH's interacting
application; it does not replace it by the paper's independent-particle
theorem. No new source dependency is paid by writing another gradient formula.

## L04 — A genuine composition law, with different object types

Blake Stephen Swistock Pollard,
[Open Markov Processes and Reaction Networks](https://arxiv.org/html/1709.09743v2),
2017, v2. Read Chapter 2 §5, Chapter 4 §§1,3,4 and Chapter 5 §4.
Chapters 4–5 develop joint work with Baez and Fong.

Open processes have explicit boundary maps. Composition identifies boundary
states and matches flows; the monoidal operation uses disjoint union.
Black-boxing describes steady boundary populations and currents and respects
that composition. Rates and boundary interfaces are input data.

**FIN assessment:** neither this category nor its energy labels supply
simultaneous FIN subsystems. Disjoint union and a joint product configuration
are different constructions. A physically typed interface and its source
are still required. This clarifies the proposed object's type instead of
reopening the exhausted tagging/incidence search.

## L05 — Conditional bridges as dynamic constraints

Giovanni Conforti, Christian Léonard,
[Reciprocal classes of random walks on graphs](https://arxiv.org/html/1505.01323v1),
2015 preprint. Read introduction and the discussion of the reference path
measure, jump intensities and reciprocal characteristics. This is a concept
screen; no full theorem application to FIN is claimed.

Sharing endpoint-conditioned bridges is a path-law constraint defined
relative to a reference walk. Reciprocal characteristics encode dynamic
information unavailable in a one-time equilibrium marginal.

**FIN assessment:** an independently supplied bridge law might constrain a
kinetic source. Computing the characteristics from the already chosen FIN
rates merely encodes the same rates. This candidate is not admitted as K07.

## L06 — Thermal time does not remove the observable-algebra premise

Alain Connes, Carlo Rovelli,
[Von Neumann algebra automorphisms and time-thermodynamics relation in
general covariant quantum theories](https://arxiv.org/html/gr-qc/9406019v1),
1994, v1. Read the stated thermal-time postulate, algebraic setup and §4.2.

The proposal uses an observable algebra and a state; its classical limit
already has a Poisson structure. It is not a universal state-to-kinetics
construction for arbitrary classical probability models.

**FIN assessment (elementary inference):** on the finite passive algebra
C^X with faithful diagonal density rho, rho^(it) f rho^(-it)=f. Its modular
action is trivial. Enlarging to noncommuting operations or supplying a
classical Poisson carrier is additional input, not a consequence of NT's
Hilbert-space notation. The older FIN monograph already mentions thermal
time as a conditional tool; this route is not selected.

## What was already present in FIN

| Prior result actually consulted | Consequence for this search |
|---|---|
| AGENTS entropic/frenetic/FDT guardrail; ST8547/8548 | Rates and current costs are not obtained from equilibrium alone; gradient reformulation is old. |
| KINETIC_NATURALITY_43 and FISHER_KINETIC_METRIC_44 | Refinement fixes some scaling shapes; a statistical metric is not a selected physical kinetic metric. |
| TH-002 and its intake | Use the admitted interacting path action and ordinary reversal directly. |
| NT intake | A transfer generator is not an operational quantum Hamiltonian. |
| HM/MM intakes and completions | Preserve nonunique labelled reset lifts, native observation memory and the actual diagnostic initialization. |
| FIN_Ten_Research_Programs_Monograph, Program 15 | Cauchy dilation and the finite-bath recurrence obstruction are already known. |
| AGENTS/post-SRN summaries of reports 148/154/200/202/203 | Copy labels, state transitions and observation partitions do not supply spatial incidence. No new inventory was launched. |
| MULTIUNIT_COMPOSITION_REENTRY_VERDICT_244 and post-255 corrections | SWAP transport requires its stated gate, conservation and incidence premises. |

The targeted duplicate search found the dilation, gradient/FDT, Fisher and
composition boundaries above. It did not locate the specific same-action
pre-GENERIC W=0 application or bound (4) in the searched current FIN proofs
and selected monographs. That bounded search is not a global novelty theorem.

## Access and provenance limitations

Firecrawl was authenticated but returned HTTP 402 on the search request;
browser-harness supplied the actual source access. The browser initially
needed its normal runtime write access outside the sandbox; escalation
succeeded. One navigation timed out while a long thesis was loading;
subsequent extraction succeeded. An attempted nonexistent MPR v3 returned
no document; arXiv metadata identified v1, whose full text was then read.
No login, purchase, remote browser rental or external message was performed.

Captures are rendered text/DOM-section artifacts, not original PDF bytes or
proofs of the cited results. Their hashes bind what was actually captured.
The MPR text companion has trailing whitespace normalized; its unmodified
DOM text is retained in `evidence/gradient_ldp_browser_raw.json`.
No missing producer ZIP or historical original is retrospectively verified.
