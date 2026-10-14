# FIN theoretical literature database

Checked online: 2026-10-14. Fifteen selected primary records; not an exhaustive
survey or a claim that no newer relevant work exists. Recent verified examples
include2024 and2025 results. Search also covered2026; unrelated results were
not added just to make the dates look newer. Bibliographic metadata, reading
scope, assumptions and FIN role are in LITERATURE_REGISTER.json.

Raw captures stay in ignored `.firecrawl/fin_theory_20261014/`. This database
contains original notes and direct references, not redistributed full papers.
Publisher abstracts behind a paywall were not treated as full proof access.

## 1. Primary methods for the proposed campaign

| ID | Source | What it offers FIN | Necessary limitation |
|---|---|---|---|
| L01 | [Kraaij,2015: interacting Markov-jump LDP](https://arxiv.org/html/1502.05185v3) | A precise theorem framework for a path cost derived from the actual count generator | Verify positivity, boundary regularity, rate convergence and initial-law LDP; visible contraction needs its own proof |
| L02 | [Mielke–Peletier–Renger,2014: gradient structures and LDP](https://link.springer.com/article/10.1007/s11118-014-9418-5), [preprint](https://arxiv.org/abs/1312.7591) | Links microscopic path fluctuations and generalized dissipation, beyond fitting the mean flow | Independent-particle examples do not constitute the interacting FIN theorem; conditions must be checked |
| L03 | [Cuff et al.,2012: mean-field Potts Glauber dynamics](https://arxiv.org/abs/1204.4503) | Spinodal/critical scaling benchmark for a first-order mean-field system | Potts and rank7 FIN have different matrices; single-update versus per-copy clocks differ |
| L09 | [Dawson,1983: critical mean-field fluctuations](https://link.springer.com/article/10.1007/BF01010922) | Demonstrates why Gaussian fixed-time approximations can fail at a critical point and a slower time scale matters | Oscillator/diffusion model, not the finite-label generator; coefficients and powers are not transferred |
| L10 | [Katzenberger,1991: strong-drift manifold limits](https://projecteuclid.org/journals/annals-of-probability/volume-19/issue-4/Solutions-of-a-Stochastic-Differential-Equation-Forced-Onto-a-Manifold/10.1214/aop/1176990225.full) | Possible tool for eliminating stable directions after rescaling | Stable-manifold/tightness assumptions need a local FIN proof, especially near exit/explosion |

Kraaij's full Theorem2.3 assumptions were inspected. The FIN rate limit
p_i softmax_j(gAp) has the natural positive boundary factor p_i, suggesting
a direct application after the remaining hypotheses are paid. This is a
**theorem application to be completed**, not a declaration that this turn
proved FIN's full path LDP. [Kraaij](https://arxiv.org/html/1502.05185v3)

The Potts paper distinguishes a dynamical spinodal from the thermodynamic
transition and has an n^(-2/3) crossover window. That motivates the FIN
double-scaling problem, but not inheritance of Potts mixing constants or
global switching results. [Cuff et al.](https://arxiv.org/abs/1204.4503)

## 2. Current research that prevents overclaiming

| ID | Source | Use | Assumption that cannot be skipped |
|---|---|---|---|
| L04 | [Cohen–Huffman,2025, revised November: uniform-in-time mean-field jump convergence](https://arxiv.org/abs/2502.20262) | A recent target for stable-regime weak error estimates | Unique exponentially stable nonlinear stationary law. Multiple FIN stationary points or a fold violate that premise |
| L05 | [Aristoff–Johnson–Perez,2023: Markov renewal/Mori–Zwanzig coarse graining](https://arxiv.org/abs/2305.20083) | Memory-aware alternatives to forcing every coarse output into a CTMC | A supplied microscopic process, observation and memory/renewal representation; no source of fundamental memory |
| L14 | [Halberstam–Hutchcroft,2024: partial-exchangeability conjecture](https://arxiv.org/abs/2405.20276) | Boundary benchmark for history/transition-count sufficiency | Transience and partial exchangeability differ from copy-permutation symmetry; finite irreducible FIN is recurrent |
| L15 | [Yasodharan–Sundaresan,2020/2021: mean field with a fast environment](https://arxiv.org/abs/2008.06855) | Comparison framework for multiscale path contraction | Do not add a new environment or speed its clock by hand; derive fast variables/scales from FIN first |

A unique **global minimum of F** is not the same as a unique stationary
solution of the nonlinear mean-field equation. Consequently the2025 theorem
cannot be applied at g3.7 merely from the MP7 global-minimum result while
other stationary branches exist. Local stopped estimates are a separate
question. [Cohen–Huffman](https://arxiv.org/abs/2502.20262)

The Mori–Zwanzig/renewal paper treats loss of dynamical information under
projection; it does not identify that memory with inertia, particles or a
physical clock. Its canonical arXiv dates, not a rendered TeX “today” date,
are used in this register. [Aristoff et al.](https://arxiv.org/abs/2305.20083)

## 3. Theory already present in FIN: credit and reuse

| ID | Source | Existing FIN connection | Do not repeat/promote |
|---|---|---|---|
| L06 | [Maas,2011: entropy flow for reversible finite chains](https://arxiv.org/abs/1102.5238) | Related to ST8547 entropy/refinement/metric discussion | A metric on probability laws is not spacetime, and a given generator is not its source |
| L07 | [Diaconis–Freedman,1980: finite exchangeable sequences](https://projecteuclid.org/journals/annals-of-probability/volume-8/issue-4/Finite-Exchangeable-Sequences/10.1214/aop/1176994663.full) | Conditional-chaos/phase-mixture and finite-marginal limits | Exchangeability alone does not make the different-N Gibbs measures exactly consistent |
| L08 | [Kurtz,1970: ODE limits of pure jump processes](https://www.cambridge.org/core/journals/journal-of-applied-probability/article/solutions-of-ordinary-differential-equations-as-limits-of-pure-jump-markov-processes/E2C5068B4B6EFF00E2FC4948FF95B152) | Background for already admitted mean-field dynamics131/138 | Deriving the same leading ODE again is not the new campaign |
| L13 | [Moreau,1965: proximal/duality framework](https://www.numdam.org/articles/10.24033/bsmf.1625/) | Context for an existing argmin-plus-square-penalty proposal | Convex-analysis hypotheses or a new cost cannot be inserted to produce a desired world |

Finite exchangeable marginals admit quantitative approximations by mixtures
of iid laws. In FIN, applying this alongside Gibbs concentration can clarify
phase mixtures versus phase-conditioned independence; neither is quantum
entanglement or a derived population of spatial cells.
[Diaconis–Freedman](https://projecteuclid.org/journals/annals-of-probability/volume-8/issue-4/Finite-Exchangeable-Sequences/10.1214/aop/1176994663.full)

## 4. Distant bridges: useful benchmarks, not active physics promotions

| ID | Source | What is permitted | What is not derived |
|---|---|---|---|
| L11 | [Chetrite–Touchette,2015: conditioned Markov processes](https://arxiv.org/abs/1405.5157) | Analyse rare-event ensembles or a Doob-transform representation of a specified record condition | A conditioned/driven process is not automatically FIN's native law or a physical drive |
| L12 | [Henley,2003: classical/quantum RK correspondence](https://arxiv.org/abs/cond-mat/0311345) | Construct an operator representation from a detailed-balance generator | Physical unitary time, quantum measurement, nonclassical operational structure or Born physics |

For FIN's reversible generator S, the matrix -S and the vector sqrt(pi)
give a mathematical positive-parent construction. Changing exp(tS) to
exp(i tS), however, changes the evolution and its interpretation. It is not
a consequence that turns the classical record law into physical QM. The
RK construction is therefore a comparison/representation lane, not a
priority programme under the no-new-assumptions requirement.
[Henley](https://arxiv.org/abs/cond-mat/0311345)

## 5. How the execution agent must use this database

Select the exact theorem needed, read its full statement/proof assumptions,
and export a condition-by-condition FIN application table. Do not claim a
theorem from an abstract, transfer a diffusion result to a jump process,
use uniform stability at a fold, equate history partial exchangeability with
copy symmetry, or count similar vocabulary as support for reality generation.

The first useful references are L01/L02/L03/L10. L04 defines a sharp scope
boundary. L06/L08 are reused background, not newly discovered FIN content.
