# FIN post-PHY task catalogue

Date: 2026-09-24. Reviewed local HEAD: `01c89aa3c24ce5840ef86e950530b6ae8caafd24`.

**Design only.** Every task is PLANNED_NOT_EXECUTED. See
[programme diagnosis](FIN_POST_PHY_RESEARCH_PROGRAM.md),
[de-dup matrix](FIN_POST_PHY_DEDUP_MATRIX.md),
[kill/source gates](FIN_POST_PHY_KILL_SOURCE_GATES.md) and
[machine-readable tasks](FIN_POST_PHY_TASKS.json).

The catalogue has 37 substantive research/protocol atoms, three audit/synthesis
nodes and two dormant promotion gates. It does not manufacture 42 new theorems.
New means a proposed extension relative to the searched FIN snapshot, not a
world-priority claim. Literature methods are benchmarks with hypotheses to pay.

## Common execution contract

Write outputs under `fin_post_phy_campaign/<ID>/`; never edit predecessor
payloads. Before research record exact source hashes, mathematical domains,
PM assumptions, observation/reset laws and numerical acceptance inequalities.
Keep PROVEN, STRONG_NUMERICAL_EVIDENCE, CONDITIONAL, HYPOTHESIS and SPECULATION
separate from execution status. A negative theorem is a successful outcome.

S: normally <=5 min/run and 2 GiB. M: <=30 min/run and 4 GiB.
L: <=2 h/task in <=30-minute checkpointed jobs, <=6 GiB.
X: separate explicit campaign/budget approval; do not launch automatically.
Resource ceilings are not proof-time promises. Use one heavy worker by default.
Stop after two failed attempts with the same cause without a new mathematical
input. A timeout or failed sufficient bound is not a refutation.

Source tasks must state the pretarget information and freeze the law/grammar
before inspecting held-out targets. Do not pick a law because it returns rank
seven, dimension three or a desired constant. A missing independently motivated
source object is a valid STOP outcome, not permission for unlimited inverse design.

Each task returns REPORT.md, a theorem/counterexample note, results.json, replay.py
and a manifest. An analytic-only replay may check exact symbolic identities;
do not invent a numerical certificate for a conceptual proposition. Protocol
tasks additionally return preregistration and raw-record schemas; they collect
no empirical data without authorization. Numeric targets require outward bounds,
not tolerances alone. Never use a local likelihood maximum as an upper bound.

Priority scores are ordinal judgments on eight explicit dimensions. They rank
scientific value, not immediate executability. Dependencies and dormant/external
gates override scores. SYN-001 consumes outcomes, including honest skips, not
an assumption that every upstream hypothesis was proved.

## AUD-001 — Admissibility of the numerical design constants

Priority: **P1**, score 38; wave 0.

- **Dependencies:** None.
- **Central question:** Which latest OCB numbers support rigorous error or cost claims?
- **Why it matters:** Prevents population design points and floating optimizers becoming power certificates.
- **Input status:** Accepted master intake; C13 supersedes C12; numerical replay variation remains.
- **Claim type:** NUMERICAL_CERTIFICATE_TARGET
- **Exact mathematical object:** Fixed matrices, singular values, calibration radii and probability intervals in C13.
- **Exact method:** Freeze hashes; recompute only decisive constants with rational/interval bounds. Distinguish r<sigma from a power proof requiring an additional deviation allowance. Check wave subprobability conventions and fitted-null membership.
- **Acceptance criterion:** Each retained number has outward endpoints, dependencies and a correct theorem/design label.
- **Failure/refutation criterion:** A supposed probability bound is negative, relies on an invalid readout, or floating error crosses a claimed sign.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Perturb final bits, replace EM likelihood by a provably suboptimal value, and test boundary singularities.
- **Stop rule:** Stop after locating a decisive defective gate; quarantine its descendants, not unrelated carrier theorems.
- **Computational class:** M
- **What may be concluded:** Admit exact constants in the declared fixture or retain only numerical-design status.
- **What must NOT be concluded:** Exact independent replay of all campaigns, global EM optimality or empirical evidence.
- **Repo de-dup / prior:** [C13:HANKEL03B_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b76d84bfbbe9f7368516c053c6ab6105f2d5ae5e5ffe60fe19b4090f52f5ea96); [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010); [C13:CHANNEL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/cf/cfab131c8cca3a29ca26cb45b1451974d08a8e8e1fa7dd799eb20719aab7af83); [fin_physics_review/MASTER_INTAKE_20260923.md](fin_physics_review/MASTER_INTAKE_20260923.md).
- **Genuinely new atom:** New outward bounds and semantic gates, not another content-hash or producer replay.
- **Artifacts:** `fin_post_phy_campaign/AUD-001/REPORT.md`, `fin_post_phy_campaign/AUD-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/AUD-001/results.json`, `fin_post_phy_campaign/AUD-001/replay.py`, `fin_post_phy_campaign/AUD-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 5, 2, 1, 5, 3].

## SRC-001 — Seven-feature observable-algebra obstruction

Priority: **P0**, score 42; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Can the retained rank-seven feature space be a closed classical observable algebra?
- **Why it matters:** Identifies the extra state/readout structure required before interpreting rank seven ontologically.
- **Input status:** Rank/source nonidentifiability and MP7 Fourier support closure are already established.
- **Claim type:** NO_GO_TARGET
- **Exact mathematical object:** B=span{1,C3,S3,C4,S4,C5,S5,C6} in functions on Z12.
- **Exact method:** Multiply character generators exactly. Classify translation-equivariant unital subalgebras as partition algebras; test whether dimension eight is possible under the transitive Z12 action. Contrast observable-product closure with log-density addition in the valid exponential family.
- **Acceptance criterion:** A closure theorem or explicit missing modes, with the minimum algebraic enlargement and scope proved.
- **Failure/refutation criterion:** Nonclosure kills the proposed exact eight-observable algebra, not the rank-seven exponential-family model.
- **Pretarget information:** Typed cyclic action and multiplication of observable functions; no desired rank enters the general classification.
- **Held-out test:** Repeat the classification on Z18 and Z24 without choosing a desired retained rank.
- **Stop rule:** Stop once an exact obstruction is found; do not adjust the required algebra to retain seven by definition.
- **Computational class:** S
- **What may be concluded:** Rank-seven statics need not be a complete classical observable algebra.
- **What must NOT be concluded:** Refutation of MP7, derivation of a different physical rank, or spatial dimension.
- **Repo de-dup / prior:** [C00:PHY_CAMPAIGN_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/e6/e619e95ca6274dc84872d358dbf5766f3c6b27e099fb1a819043c515351e8149); [C05:CARRIER_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/73/7388582782c815184dadbd20bdecf349b1c0af9c0d73cfd61dad783cc75d2572); [fin_rank7_mathphysics_next/proofs/MP7-013_boundary_support_classification.md](fin_rank7_mathphysics_next/proofs/MP7-013_boundary_support_classification.md).
- **Genuinely new atom:** Tests pointwise algebra/coarse-partition closure, not PSD symmetry or largest-gap selection again.
- **Artifacts:** `fin_post_phy_campaign/SRC-001/REPORT.md`, `fin_post_phy_campaign/SRC-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/SRC-001/results.json`, `fin_post_phy_campaign/SRC-001/replay.py`, `fin_post_phy_campaign/SRC-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 5, 1, 1, 5, 4].

## SRC-002 — Compositional spectral-cut rigidity

Priority: **P1**, score 38; wave 1.

- **Dependencies:** SRC-001.
- **Central question:** Which spectral selections commute exactly with independent composition and refinement?
- **Why it matters:** Tests whether naturality can source a nontrivial rank cut or instead exposes an unavoidable resolution premise.
- **Input status:** C04 already derives quadratic-law rigidity; spectral eigengap and penalized-rank rules are not sources.
- **Claim type:** NO_GO_TARGET
- **Exact mathematical object:** Spectral filters P(A), including the constant sector, under direct sums and Kronecker-sum embeddings.
- **Exact method:** State each naturality axiom separately. Derive its scalar functional equation; distinguish exact idempotent filters from approximate filters. Search for countermodels when an axiom is removed, not for a cut producing seven.
- **Acceptance criterion:** Exact classification or a finite counterexample showing which composition axiom forbids a nontrivial cut.
- **Failure/refutation criterion:** If only trivial filters survive, reject that rank-source class; do not declare all possible source laws impossible.
- **Pretarget information:** Composition maps and spectral calculus axioms fixed before inspecting the FIN cut.
- **Held-out test:** Unequal factor dimensions and spectra with different multiplicities.
- **Stop rule:** No arbitrary new filter penalty after failure; export the missing resolution or observation datum.
- **Computational class:** M
- **What may be concluded:** A premise-level rank-source obstruction or conditional filter theorem.
- **What must NOT be concluded:** That the chosen naturality axioms are themselves laws of nature.
- **Repo de-dup / prior:** [C02:RANK_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/41/4167e4d4ae66b0caa49c2073cb9dacc3de8d9c7063fc7ca9b3ca8498317f487d); [C04:PROTECT_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/5c/5c07d04f7b05ee444a810a04bc80569fa46737583b48d95dcfa355a8d240aa65); [C06:results/OPERATOR-REFINE-07_intertwining_moduli.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/69/69560c85308f541c472d9433259178afe03d12ab6fb005af968a09576b6296c6).
- **Genuinely new atom:** Classifies projector/filter functors rather than repeating the nonlinear-energy Cauchy equation.
- **Artifacts:** `fin_post_phy_campaign/SRC-002/REPORT.md`, `fin_post_phy_campaign/SRC-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/SRC-002/results.json`, `fin_post_phy_campaign/SRC-002/replay.py`, `fin_post_phy_campaign/SRC-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 5, 5, 5, 2, 2, 5, 4].

## SRC-003 — Minimal operational data selecting a state cone

Priority: **P1**, score 37; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** What extra preparation/effect data remove the abstract-operator simplex ambiguity?
- **Why it matters:** Separates supplied vertex ontology from a genuinely identifiable operational state space.
- **Input status:** C05 supplies isospectral Markov-cone ambiguity; full-state tomography is not yet a physical source.
- **Claim type:** IDENTIFIABILITY_TARGET
- **Exact mathematical object:** An operator A together with a finite preparation/effect algebra and its residual similarity gauge.
- **Exact method:** Compute the gauge preserving probabilities and A. Add candidate rank-one preparations/effects one at a time; identify remaining equivalences and require held-out predictions. Distinguish a measured algebra from an axiomatically imposed basis.
- **Acceptance criterion:** A necessary/sufficient finite observation package within a declared class, or indistinguishable alternative cones.
- **Failure/refutation criterion:** Distinct cones still predict identical admitted records; ontology remains unidentifiable in that experiment.
- **Pretarget information:** Available preparations/effects and their independent calibration, not a preferred simplex chosen afterward.
- **Held-out test:** Relabelings, degenerate effects and a noncirculant operator perturbation.
- **Stop rule:** Stop at the first proved nonidentifiability class; do not fill missing observations with assumed vertex labels.
- **Computational class:** M
- **What may be concluded:** Conditional operational identifiability up to a stated gauge.
- **What must NOT be concluded:** Unique fundamental ontology or QW-2191 discharge.
- **Repo de-dup / prior:** [C05:STATE_DYN_PROX_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c5/c588238d4f51da564e091f141c8119d52240c3e58b08041c5cceccb88f7eb398); [C04:DYN_OBS_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/5c/5ca93eb709498437e13035f47bf350cf2c8e335718db618c0d6e01ff5ca99d65); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Quantifies the minimal stabilizer-breaking observation package, not another abstract-A counterexample.
- **Artifacts:** `fin_post_phy_campaign/SRC-003/REPORT.md`, `fin_post_phy_campaign/SRC-003/theorem_or_counterexample.md`, `fin_post_phy_campaign/SRC-003/results.json`, `fin_post_phy_campaign/SRC-003/replay.py`, `fin_post_phy_campaign/SRC-003/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 4, 5, 4, 2, 2, 5, 4].

## SRC-004 — Resolution-relative rank without ontology claims

Priority: **P2**, score 32; wave 2.

- **Dependencies:** SRC-001, SRC-003.
- **Central question:** Can an observable error budget determine a useful effective rank across carriers?
- **Why it matters:** Replaces unsourced rank numerology by a conditional, falsifiable approximation theorem.
- **Input status:** Eckart–Young conditional truncation and disagreeing selection criteria are old results.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** A family of observable maps and truncations with a specified error tolerance.
- **Exact method:** Bound discarded-mode effects on responses and phase observables using gaps/derivative bounds. Establish when spectral-norm error is insufficient near folds. Treat tolerance, preparation family and gain window as inputs.
- **Acceptance criterion:** A rank/error tradeoff that predicts required rank on a held-out family without changing its tolerance.
- **Failure/refutation criterion:** Selected rank changes arbitrarily with unrecorded observables or cannot control the promised physical statistic.
- **Pretarget information:** Observable accuracy and preparations fixed independently of the rank eventually selected.
- **Held-out test:** q=18,24 and a gain interval not used to set the tolerance.
- **Stop rule:** If no stable rank window exists, retain an adaptive effective model; never tune tolerance to seven.
- **Computational class:** M
- **What may be concluded:** Effective resolution-dependent rank with certified errors.
- **What must NOT be concluded:** An intrinsic rank-seven source or a universal number of physical degrees of freedom.
- **Repo de-dup / prior:** [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db); [C02:RANK_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/41/4167e4d4ae66b0caa49c2073cb9dacc3de8d9c7063fc7ca9b3ca8498317f487d); [C00:PHY_CAMPAIGN_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/e6/e619e95ca6274dc84872d358dbf5766f3c6b27e099fb1a819043c515351e8149).
- **Genuinely new atom:** Links retained rank to held-out nonlinear response errors, not matrix approximation error alone.
- **Artifacts:** `fin_post_phy_campaign/SRC-004/REPORT.md`, `fin_post_phy_campaign/SRC-004/theorem_or_counterexample.md`, `fin_post_phy_campaign/SRC-004/results.json`, `fin_post_phy_campaign/SRC-004/replay.py`, `fin_post_phy_campaign/SRC-004/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 4, 4, 2, 2, 5, 3].

## SRC-005 — Admit one genuinely new parent source law

Priority: **P2**, score 32; wave 2.

- **Dependencies:** SRC-001, SRC-002, SRC-003.
- **Central question:** Does any independently motivated parent law explain a held-out FIN operator feature?
- **Why it matters:** Keeps SOURCE-R7-01 open without disguising inverse design as derivation.
- **Input status:** Passive connected Laplacians, generic symmetry and prior rank selectors fail; no qualifying parent law is currently exported.
- **Claim type:** CONDITIONAL_MODEL
- **Exact mathematical object:** One declared constrained/signed/hidden-variable parent action and its elimination map.
- **Exact method:** Before unblinding targets, freeze primitive variables, permitted terms and dimensionless parameters. Prove positivity/constraints and derive the effective operator. Expose any embedded projector or target coefficient. If no independent candidate is provided, issue NO_ADMISSIBLE_SOURCE rather than inventing one.
- **Acceptance criterion:** A new typed law passes target-blind admission and predicts an unused carrier observable; source status of its premises stays explicit.
- **Failure/refutation criterion:** The desired operator or rank is encoded in coefficients, constraints, topology or calibration.
- **Pretarget information:** Only the recovered relation type, independently justified parent primitives and predeclared law grammar; hide rank target and held-out spectrum.
- **Held-out test:** A predeclared q and fine-mode response outside the construction data, with no refit.
- **Stop rule:** At most one new law; class expansion requires new justification and human review.
- **Computational class:** X
- **What may be concluded:** A conditional source candidate or a precise absence-of-source result.
- **What must NOT be concluded:** Strict derivation from a freely chosen action; legacy-role transfer.
- **Repo de-dup / prior:** [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db); [C05:CARRIER_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/73/7388582782c815184dadbd20bdecf349b1c0af9c0d73cfd61dad783cc75d2572); [C13:NEXT_ATOMS.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/ad/ad3e78c61ca4ce5d4540bf9f3c922968cdeaa7aef8363b0b07d1a870e8590322); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** A finite new-law admission challenge, not an exhaustive replay of rejected source classes.
- **Artifacts:** `fin_post_phy_campaign/SRC-005/REPORT.md`, `fin_post_phy_campaign/SRC-005/theorem_or_counterexample.md`, `fin_post_phy_campaign/SRC-005/results.json`, `fin_post_phy_campaign/SRC-005/replay.py`, `fin_post_phy_campaign/SRC-005/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately; separate bounded campaign/budget approval required.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 5, 5, 2, 5, 3, 5, 5].

## REF-001 — Finite-noisy-data continuum-order impossibility

Priority: **P0**, score 41; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Can finite precision, size and observation time ever identify the asymptotic operator order without a tail class?
- **Why it matters:** Turns the old exact finite-prefix no-go into a quantitative experimental limit.
- **Input status:** C06 already gives identical finite shells with different low-k exponents.
- **Claim type:** NO_GO_TARGET
- **Exact mathematical object:** Positive tail families with bounded observable differences but distinct asymptotic symbols.
- **Exact method:** Construct delayed/small tails; propagate operator differences to finite-time records with semigroup bounds, then derive Le Cam-type limits. Separate a tail amplitude/onset lower bound from positivity and monotonicity.
- **Acceptance criterion:** An explicit experiment-dependent indistinguishability bound and the minimum additional asymptotic restriction needed.
- **Failure/refutation criterion:** Uniform exponent identification fails when crossover can be moved beyond every admitted scale.
- **Pretarget information:** Finite record budget, errors and observation windows; no preferred continuum exponent.
- **Held-out test:** Two observation times and two carrier sizes not used in the construction.
- **Stop rule:** Stop at a constructive indistinguishability theorem; no brute-force far-shell search.
- **Computational class:** M
- **What may be concluded:** A scoped impossibility theorem or an identifiable restricted-tail experiment.
- **What must NOT be concluded:** That finite-range or fractional behavior is physically selected.
- **Repo de-dup / prior:** [C06:results/OPERATOR-REFINE-03_universal_no_go.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/3f/3f6b97b2aad8ed8c11888c55586a53d5133dbb43213f769edbeac91be33c0644); [C06:results/OPERATOR-REFINE-05_fixed_stencil.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/73/730bccfabe4bb169cb59d97fdf974876184161c857984763a110fee0f9d7f8a0); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Adds finite-window statistical indistinguishability and crossover-scale obligations; does not append another shell.
- **Artifacts:** `fin_post_phy_campaign/REF-001/REPORT.md`, `fin_post_phy_campaign/REF-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/REF-001/results.json`, `fin_post_phy_campaign/REF-001/replay.py`, `fin_post_phy_campaign/REF-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 5, 2, 1, 5, 4].

## REF-002 — Dynamic subdivision forces a memory object

Priority: **P0**, score 41; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Can static subdivision-compatible stiffness and additive mass close dynamically with finitely many edge parameters?
- **Why it matters:** Links refinement to hidden-state realization rather than treating memory as optional.
- **Input status:** c(ell)=kappa0/ell and static elimination are known; dynamic Schur memory is known in finite controls.
- **Claim type:** NO_GO_TARGET
- **Exact mathematical object:** Two-terminal dynamic Dirichlet-to-Neumann maps and their Schur composition.
- **Exact method:** Eliminate one interior massive node symbolically; track poles under repeated refinement. Test closure of stiffness-plus-lumped-mass ansatz for all frequencies. Identify cancellations and minimal passive memory order; do not infer hidden physical states from linear realization order.
- **Acceptance criterion:** Exact closure classification or a pole-growth obstruction plus a minimal enlarged response type.
- **Failure/refutation criterion:** The frequency-independent ansatz cannot reproduce the boundary transfer on an open frequency interval.
- **Pretarget information:** Static edge law, declared kinetic premise and composition rule, not a desired wave equation.
- **Held-out test:** Unequal and mixed binary/ternary subdivisions.
- **Stop rule:** Stop the constant-parameter ansatz at the first exact obstruction; move only to the specified memory extension.
- **Computational class:** M
- **What may be concluded:** Refinement requires a response/memory object in the declared inertial model.
- **What must NOT be concluded:** A sourced physical inertia, clock or literal microscopic memory ontology.
- **Repo de-dup / prior:** [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db); [C06:results/OPERATOR-REFINE-06_local_schur_no_go.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/1f/1f7b2707b0f3cd2f4d5448fb4843ae74d3d5ba952916cef451b38dbca47a2fb2); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Classifies exact frequency-dependent closure across every subdivision, beyond one C8 example.
- **Artifacts:** `fin_post_phy_campaign/REF-002/REPORT.md`, `fin_post_phy_campaign/REF-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/REF-002/results.json`, `fin_post_phy_campaign/REF-002/replay.py`, `fin_post_phy_campaign/REF-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 4, 2, 1, 5, 5].

## REF-003 — Classify passive two-port refinement semigroups

Priority: **P1**, score 38; wave 2.

- **Dependencies:** REF-002.
- **Central question:** What is determined by homogeneous segment composition once memory is retained?
- **Why it matters:** Tests whether added naturality reduces operator freedom or only repackages it.
- **Input status:** REF-002 identifies the required boundary-response type; static Cauchy rigidity alone is insufficient.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Reciprocal passive 2x2 transfer matrices T_ell(z) with T_(a+b)=T_a T_b.
- **Exact method:** Derive the regular semigroup generator, impose passivity/reciprocity and recover the static limit. Classify which frequency functions remain free. Use a uniform string as a conditional example, not the selected answer.
- **Acceptance criterion:** Necessary/sufficient generator conditions and a precise list of unsourced temporal/material functions.
- **Failure/refutation criterion:** Multiple passive generators share the static law but predict different dynamic continua.
- **Pretarget information:** Boundary reciprocity, passivity, concatenation and stated regularity before selecting a dispersion target.
- **Held-out test:** Irregular meshes and frequencies withheld from fitting.
- **Stop rule:** Do not add regularity or constitutive axioms merely to force the standard wave generator.
- **Computational class:** L
- **What may be concluded:** A conditional uniqueness theorem or sharper dynamic-source nonidentifiability.
- **What must NOT be concluded:** Selection of heat/wave/unitary physics or SI speed.
- **Repo de-dup / prior:** [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db); [C04:PROTECT_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/5c/5c07d04f7b05ee444a810a04bc80569fa46737583b48d95dcfa355a8d240aa65); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Joint dynamic transfer/positivity classification, not the already proved static resistance equation.
- **Artifacts:** `fin_post_phy_campaign/REF-003/REPORT.md`, `fin_post_phy_campaign/REF-003/theorem_or_counterexample.md`, `fin_post_phy_campaign/REF-003/results.json`, `fin_post_phy_campaign/REF-003/replay.py`, `fin_post_phy_campaign/REF-003/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 3, 3, 2, 5, 5].

## REF-004 — Circular-planar realization challenge for the full strict operator

Priority: **P1**, score 39; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Can all six strict shells be a boundary response of a positive circular-planar network?
- **Why it matters:** Tests a concrete local hidden-geometry alternative left open by the one-dimensional chain no-go.
- **Input status:** The unlabeled cycle is recovered; pure edge subdivision cannot produce nonnearest couplings.
- **Claim type:** NUMERICAL_CERTIFICATE_TARGET
- **Exact mathematical object:** The full strict Laplacian as a candidate boundary response matrix in the recovered cyclic order.
- **Exact method:** Fix sign/zero-mode conventions; derive or apply the exact circular-minor characterization. Enclose all decisive minors from strict spectral/kernel inputs. If feasible, construct a certified network in the admitted class; if infeasible export one signed minor witness.
- **Acceptance criterion:** A realization certificate or an exact circular-planarity obstruction with boundary-order dependence explicit.
- **Failure/refutation criterion:** A required circular minor has the forbidden sign; circular-planar passive interpretation fails.
- **Pretarget information:** Only the full typed W and its intrinsically recovered cyclic order; no target spatial dimension or chosen hidden circuit.
- **Held-out test:** A perturbed strict tuple and a withheld boundary-current experiment; reconstruction is not uniqueness.
- **Stop rule:** Stop on a verified forbidden minor; do not insert nonplanar edges and keep the same claim.
- **Computational class:** L
- **What may be concluded:** Compatibility or incompatibility with a specified hidden passive geometry.
- **What must NOT be concluded:** Physical dimension two, unique network, spacetime or an A7 rank source.
- **Repo de-dup / prior:** [C05:CARRIER_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/73/7388582782c815184dadbd20bdecf349b1c0af9c0d73cfd61dad783cc75d2572); [C06:results/OPERATOR-REFINE-06_local_schur_no_go.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/1f/1f7b2707b0f3cd2f4d5448fb4843ae74d3d5ba952916cef451b38dbca47a2fb2).
- **Genuinely new atom:** Tests circular-planar response constraints, not arbitrary graph fitting or another chain subdivision.
- **Artifacts:** `fin_post_phy_campaign/REF-004/REPORT.md`, `fin_post_phy_campaign/REF-004/theorem_or_counterexample.md`, `fin_post_phy_campaign/REF-004/results.json`, `fin_post_phy_campaign/REF-004/replay.py`, `fin_post_phy_campaign/REF-004/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 4, 3, 1, 5, 4].

## REF-005 — Fine-response experiment breaking refinement fibers

Priority: **P1**, score 37; wave 2.

- **Dependencies:** REF-002, REF-004.
- **Central question:** Which smallest new observable distinguishes refinements invisible to all coarse records?
- **Why it matters:** Converts nonidentifiability into an executable observation obligation.
- **Input status:** Exact coarse intertwining leaves fine rates free; boundary realization may add hidden modes.
- **Claim type:** IDENTIFIABILITY_TARGET
- **Exact mathematical object:** Fine-sector preparation/effect maps and resolvent or transient transfer functions.
- **Exact method:** Compute observable/controllable quotients; minimize added probe rank or frequency information needed to separate the surviving parameters. Provide equivalent-model counterexamples when probes are insufficient.
- **Acceptance criterion:** A finite minimal-probe theorem within a specified realization class and a held-out transfer prediction.
- **Failure/refutation criterion:** All admitted probes remain constant on the refinement fiber.
- **Pretarget information:** Available probe algebra and cost, fixed before choosing the refinement law.
- **Held-out test:** An unused refinement level and nonuniform split.
- **Stop rule:** No further coarse-only statistics after exact blindness is proved.
- **Computational class:** M
- **What may be concluded:** Conditional observability of refinement parameters.
- **What must NOT be concluded:** Physical existence or unique microscopic interpretation of hidden modes.
- **Repo de-dup / prior:** [C06:results/OPERATOR-REFINE-07_intertwining_moduli.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/69/69560c85308f541c472d9433259178afe03d12ab6fb005af968a09576b6296c6); [C03:HDIM_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/03/038ff051203ba1da89646f17562409a79c2d98004606020019038cbbea296ee4); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Designs probes for the new response/refinement family, not rediscovery of the old free fiber.
- **Artifacts:** `fin_post_phy_campaign/REF-005/REPORT.md`, `fin_post_phy_campaign/REF-005/theorem_or_counterexample.md`, `fin_post_phy_campaign/REF-005/results.json`, `fin_post_phy_campaign/REF-005/replay.py`, `fin_post_phy_campaign/REF-005/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 4, 5, 4, 2, 2, 5, 4].

## REF-006 — Joint operator/refinement convergence contract

Priority: **P2**, score 31; wave 3.

- **Dependencies:** REF-003, REF-005.
- **Central question:** Which response-consistent families have a controlled continuum rather than only a convergent metric?
- **Why it matters:** Prevents transfer of the backbone circle limit to the dense full operator.
- **Input status:** Metric subdivision and fixed-stencil convergence are conditional known results.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Dirichlet forms, measures and frequency-dependent boundary maps across a declared refinement tower.
- **Exact method:** Specify topology and time/mesh scaling; prove Mosco/resolvent convergence or exhibit a defect. Keep spatially fractional and temporally nonlocal alternatives distinct. Treat fixed finite rank as an internal truncation, not a spatial differential operator.
- **Acceptance criterion:** Uniform convergence/error estimates on a declared observable/frequency class, or a scoped no-go.
- **Failure/refutation criterion:** Metric convergence coexists with incompatible operators or uncontrolled memory.
- **Pretarget information:** Admitted refinement law and probes, not the desired limiting symbol.
- **Held-out test:** Mixed subdivisions and fine modes outside the development window.
- **Stop rule:** No continuum claim from fitted slopes alone; X-scale proof requires a separate budget.
- **Computational class:** X
- **What may be concluded:** A conditional continuum family with explicit residual source assumptions.
- **What must NOT be concluded:** Uniqueness selected by FIN or physical space/time.
- **Repo de-dup / prior:** [C05:CARRIER_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/73/7388582782c815184dadbd20bdecf349b1c0af9c0d73cfd61dad783cc75d2572); [C06:results/OPERATOR-REFINE-05_fixed_stencil.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/73/730bccfabe4bb169cb59d97fdf974876184161c857984763a110fee0f9d7f8a0); [C06:results/OPERATOR-REFINE-07_intertwining_moduli.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/69/69560c85308f541c472d9433259178afe03d12ab6fb005af968a09576b6296c6).
- **Genuinely new atom:** Checks joint response, measure, memory and operator convergence under the newly admitted law.
- **Artifacts:** `fin_post_phy_campaign/REF-006/REPORT.md`, `fin_post_phy_campaign/REF-006/theorem_or_counterexample.md`, `fin_post_phy_campaign/REF-006/results.json`, `fin_post_phy_campaign/REF-006/replay.py`, `fin_post_phy_campaign/REF-006/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately; separate bounded campaign/budget approval required.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 5, 3, 5, 2, 5, 4].

## PHA-001 — Uniform transverse phase-locking gap

Priority: **P1**, score 34; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Is qualitative all-orders locking quantitatively stable across amplitudes and carrier size?
- **Why it matters:** A uniform gap is needed before eliminating relative phases or treating a translation coordinate as slow.
- **Input status:** C04/C05 already prove relation lattices and equality sets; they do not by themselves give every required uniform gap.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Phase Hessians on a lattice quotient with explicit amplitude windows.
- **Exact method:** Lower-bound selected positive Bessel coefficients and the relation-frame quadratic form. Separate transverse and common-shift eigenvalues; carry the metric. Determine how bounds deteriorate when a mode vanishes.
- **Acceptance criterion:** Uniform bounds on a declared compact nonzero-amplitude family, with exact failure strata.
- **Failure/refutation criterion:** No uniform gap survives the proposed window or the selected quotient changes at a boundary.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** q values not used in discovery and unequal modal amplitudes.
- **Stop rule:** Do not remove weak modes silently or interpret loss of a sufficient bound as loss of all locking.
- **Computational class:** M
- **What may be concluded:** Controlled phase reduction in the specified model family.
- **What must NOT be concluded:** Additional spatial dimensions, Goldstone particles or spontaneous protection of all radii.
- **Repo de-dup / prior:** [C04:RESONANCE_DIMENSION_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/76/762a8878309629e01aa6d3bb01e527971f27c4d143036e1f49150e56a079ce65); [C05:RESONANCE_CERT_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/6e/6ee7134de2107784707a61503619b7b5b577113afd39d4deaafb5417bc328655); [C03:HDIM_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/03/038ff051203ba1da89646f17562409a79c2d98004606020019038cbbea296ee4).
- **Genuinely new atom:** Pays quantitative conditioning and amplitude-zero strata, not another alignment/equality proof.
- **Artifacts:** `fin_post_phy_campaign/PHA-001/REPORT.md`, `fin_post_phy_campaign/PHA-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/PHA-001/results.json`, `fin_post_phy_campaign/PHA-001/replay.py`, `fin_post_phy_campaign/PHA-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 4, 5, 2, 2, 5, 4].

## PHA-002 — Finite-q pinning versus finite-N resolution

Priority: **P1**, score 35; wave 2.

- **Dependencies:** PHA-001.
- **Central question:** Which joint q,N regime resolves local phases while washing out Umklapp pinning?
- **Why it matters:** Connects the exact finite carrier with the continuum translation manifold quantitatively.
- **Input status:** The first Umklapp order and q/M² collision gate are already known.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** The common-shift free-energy barrier Delta_q and the scaled quantity N Delta_q.
- **Exact method:** Bound alias tails uniformly using complex-strip/Bessel estimates. Compare N Delta_q with transverse fluctuation and phase-collision scales; prove a joint window or incompatible requirements. Fix amplitude/gain hypotheses before optimizing q.
- **Acceptance criterion:** An explicit joint-limit theorem with remainder and a separate pinned regime.
- **Failure/refutation criterion:** The chosen scaling cannot simultaneously satisfy phase coherence, diffuse geometry and pinning control.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Mixed refinement schedules and q values outside powers of two.
- **Stop rule:** No continuum-symmetry claim from truncating every Umklapp term at a fixed order.
- **Computational class:** L
- **What may be concluded:** Conditional quantitative finite-to-continuum phase correspondence.
- **What must NOT be concluded:** A physical continuum selected by FIN or disappearance of finite-q symmetry at fixed q.
- **Repo de-dup / prior:** [C04:RESONANCE_DIMENSION_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/76/762a8878309629e01aa6d3bb01e527971f27c4d143036e1f49150e56a079ce65); [C05:RESONANCE_CERT_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/6e/6ee7134de2107784707a61503619b7b5b577113afd39d4deaafb5417bc328655); [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db).
- **Genuinely new atom:** Controls the full pinning barrier and noncommuting limits, not only its first Taylor order.
- **Artifacts:** `fin_post_phy_campaign/PHA-002/REPORT.md`, `fin_post_phy_campaign/PHA-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/PHA-002/results.json`, `fin_post_phy_campaign/PHA-002/replay.py`, `fin_post_phy_campaign/PHA-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 4, 4, 4, 3, 1, 5, 4].

## PHA-003 — Observable mesoscopic phase-symmetry window

Priority: **P1**, score 38; wave 2.

- **Dependencies:** PHA-001, SRC-003.
- **Central question:** Can approximate T5 symmetry coexist with statistically resolvable phases?
- **Why it matters:** Tests a new bridge between entropy-breaking estimates and observation information.
- **Input status:** C04 gives O(r³) entropy breaking; C03/C04 establish observable relative phases and nonprotected radii.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** A declared fixed-radius or annular preparation family with copy scale N, readout count n and amplitude r.
- **Exact method:** Test N r³→0 against n r²→infinity with exact normalization. The illustrative n=N window is N^-1/2 << r << N^-1/3. Prove or reject uniform likelihood and angular-measure bounds; record preparation, coherence lifetime and radial conditioning separately.
- **Acceptance criterion:** A nonempty, observable conditional window with explicit error and no hidden identification of N and n.
- **Failure/refutation criterion:** Noise, angular Jacobians, radial fluctuations or short lifetime destroy the simultaneous window.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Nonuniform amplitudes, independent n,N scaling and a held-out phase contrast.
- **Stop rule:** If annular conditioning is essential, keep it PM; do not call the window a spontaneous stable HD phase.
- **Computational class:** M
- **What may be concluded:** A mesoscopic operational symmetry regime or a precise incompatibility theorem.
- **What must NOT be concluded:** Four/five physical dimensions, exact nonlinear symmetry or a physical information source.
- **Repo de-dup / prior:** [C04:RESONANCE_DIMENSION_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/76/762a8878309629e01aa6d3bb01e527971f27c4d143036e1f49150e56a079ce65); [C04:DYN_OBS_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/5c/5ca93eb709498437e13035f47bf350cf2c8e335718db618c0d6e01ff5ca99d65); [C03:HDIM_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/03/038ff051203ba1da89646f17562409a79c2d98004606020019038cbbea296ee4).
- **Genuinely new atom:** Combines symmetry-breaking scale and Fisher/readout scale, not a repeat of shrinking-amplitude scans.
- **Artifacts:** `fin_post_phy_campaign/PHA-003/REPORT.md`, `fin_post_phy_campaign/PHA-003/theorem_or_counterexample.md`, `fin_post_phy_campaign/PHA-003/results.json`, `fin_post_phy_campaign/PHA-003/replay.py`, `fin_post_phy_campaign/PHA-003/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 4, 2, 1, 5, 4].

## GEO-001 — Positive-weight consistency cone on RNG support

Priority: **P1**, score 39; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Can scalar conductances on the existing sparse support repair its local anisotropy?
- **Why it matters:** Kills impossible weighting programmes before expensive continuum simulation.
- **Input status:** Raw RNG anisotropy and structured-grid metric obstruction already exist.
- **Claim type:** NO_GO_TARGET
- **Exact mathematical object:** Weights satisfying drift cancellation and second-moment tensor consistency on fixed incident directions.
- **Exact method:** Form exact/interval linear feasibility for sum w_e d_e=0 and sum w_e d_e d_e^T=target tensor. Use Farkas dual witnesses. Compare pointwise consistency with weaker homogenization; explicitly supply volume measure if needed.
- **Acceptance criterion:** A positive-weight construction with bounds or a certified separating dual on a declared generic class.
- **Failure/refutation criterion:** The target tensor lies outside the support cone; pointwise repair is impossible on that support.
- **Pretarget information:** Metric sample and measure supplied before selecting weights; no desired dimension is inferred from the consistency target.
- **Held-out test:** Rectangular tori with unseen aspect ratios and iid point populations.
- **Stop rule:** Do not change the graph and still claim the fixed-support weighting problem passed.
- **Computational class:** M
- **What may be concluded:** A feasibility theorem or support-specific obstruction.
- **What must NOT be concluded:** Physical dimension or unique isotropic geometry from Haar density alone.
- **Repo de-dup / prior:** [C04:PROX04_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/ed/edeb612adef3d894b9392f15ed3af7145e7b47cfcc8a0b81402d85dcb744a0e6); [C05:STATE_DYN_PROX_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c5/c588238d4f51da564e091f141c8119d52240c3e58b08041c5cceccb88f7eb398).
- **Genuinely new atom:** Tests all positive weights by local moment-cone feasibility, not another unweighted random scout.
- **Artifacts:** `fin_post_phy_campaign/GEO-001/REPORT.md`, `fin_post_phy_campaign/GEO-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/GEO-001/results.json`, `fin_post_phy_campaign/GEO-001/replay.py`, `fin_post_phy_campaign/GEO-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 5, 2, 1, 5, 4].

## GEO-002 — Homogenization of a declared relational population

Priority: **P2**, score 29; wave 3.

- **Dependencies:** GEO-001, PHA-003.
- **Central question:** Can weak continuum behavior survive when pointwise RNG consistency fails?
- **Why it matters:** Keeps the correct probabilistic continuum question separate from an impossible local tensor limit.
- **Input status:** C05 proves Haar empirical density alone does not fix the operator; no iid/Haar continuum theorem is accepted.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Random conductance forms on a specified RNG or minimally enlarged proximity graph.
- **Exact method:** State sampling and ellipticity/moment assumptions. Seek a corrector/invariance-principle or form-convergence proof. Test whether FIN-induced correlated phase populations satisfy them; do not substitute iid sampling unnoticed.
- **Acceptance criterion:** A conditional limit theorem or a trap/anisotropy counterexample within the declared population class.
- **Failure/refutation criterion:** Missing mixing, unbounded traps or degeneracy defeats the proposed scaling.
- **Pretarget information:** Population law and conductance rules fixed without using a desired continuum exponent or dimension.
- **Held-out test:** Correlated populations outside the proof class as negative controls; new aspect ratios inside it.
- **Stop rule:** No larger iid simulations after a structural assumption fails; identify the missing population datum.
- **Computational class:** X
- **What may be concluded:** Conditional homogenized response, not pointwise isotropy.
- **What must NOT be concluded:** A sourced sampling law or unique physical metric.
- **Repo de-dup / prior:** [C04:PROX04_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/ed/edeb612adef3d894b9392f15ed3af7145e7b47cfcc8a0b81402d85dcb744a0e6); [C05:STATE_DYN_PROX_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c5/c588238d4f51da564e091f141c8119d52240c3e58b08041c5cceccb88f7eb398); [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db).
- **Genuinely new atom:** Adds a precise stationary/mixing population law and quenched/annealed target instead of extrapolating density plots.
- **Artifacts:** `fin_post_phy_campaign/GEO-002/REPORT.md`, `fin_post_phy_campaign/GEO-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/GEO-002/results.json`, `fin_post_phy_campaign/GEO-002/replay.py`, `fin_post_phy_campaign/GEO-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately; separate bounded campaign/budget approval required.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 4, 3, 5, 2, 5, 4].

## GEO-003 — Reconcile operator geometry with path geometry

Priority: **P2**, score 32; wave 3.

- **Dependencies:** GEO-001, REF-006.
- **Central question:** Does a response-defined intrinsic metric converge where shortest-path distance fails?
- **Why it matters:** Offers a testable repair of the known grid L1-versus-Euclidean mismatch.
- **Input status:** C04 already separates Laplacian and shortest-path limits.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Energy-measure intrinsic distances and small-time response distances on a convergent graph family.
- **Exact method:** Compare the intrinsic metric of the limiting quadratic form with finite graph metrics. Analyze mesh→0 before time→0 versus the reverse; fixed finite-graph heat kernels need not have the continuum short-time distance law.
- **Acceptance criterion:** A metric/operator compatibility theorem with explicit limits and measurable response map, or a no-go.
- **Failure/refutation criterion:** Form convergence alone fails to control distance or the proposed finite-data estimator is nonidentifiable.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Rectangular and mildly irregular grids plus the admitted random family if available.
- **Stop rule:** Stop physical-space language until the observation map and scale regime are paid.
- **Computational class:** L
- **What may be concluded:** A conditional response geometry with a controlled limit.
- **What must NOT be concluded:** That changing the metric definition derives spacetime.
- **Repo de-dup / prior:** [C04:PROX04_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/ed/edeb612adef3d894b9392f15ed3af7145e7b47cfcc8a0b81402d85dcb744a0e6); [C05:STATE_DYN_PROX_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c5/c588238d4f51da564e091f141c8119d52240c3e58b08041c5cceccb88f7eb398); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Derives a different operational metric and its order-of-limits conditions; does not rename the old graph distance.
- **Artifacts:** `fin_post_phy_campaign/GEO-003/REPORT.md`, `fin_post_phy_campaign/GEO-003/theorem_or_counterexample.md`, `fin_post_phy_campaign/GEO-003/results.json`, `fin_post_phy_campaign/GEO-003/replay.py`, `fin_post_phy_campaign/GEO-003/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 4, 4, 3, 2, 5, 4].

## DYN-001 — Causality under state-dependent local rewiring

Priority: **P1**, score 36; wave 3.

- **Dependencies:** PHA-002, GEO-001.
- **Central question:** Can an evolving phase-order graph retain an influence bound without a uniform dwell time?
- **Why it matters:** Addresses the actual state-sourced graph, not a frozen support that assumes away its main difficulty.
- **Input status:** Adjacent crossings alter finitely many edges, but Haar minimum gaps scale as M^-2 and uniform dwell is absent.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Time-ordered local generators on a phase-sorted graph with explicit collision/update rules.
- **Exact method:** Define event supports and pathwise coupling. Bound chronological influence using angular transport and update intensity; examine whether singular conductances at small gaps destroy uniformity. Charge every collision rule or hard-core premise.
- **Acceptance criterion:** A size-uniform conditional bound in a stated metric or an explicit accelerated-influence construction.
- **Failure/refutation criterion:** Rewiring/global sorting creates uncontrolled remote influence or requires an unrecorded collision regularizer.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Near-collision initial data and mixed update schedules.
- **Stop rule:** No uniform-dwell assumption may be silently reintroduced after its known failure.
- **Computational class:** L
- **What may be concluded:** Conditional causal control for a specified evolving network.
- **What must NOT be concluded:** A physical clock, exact Lorentz cone or sourced motion law.
- **Repo de-dup / prior:** [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db); [C00:RELATIONAL_CONTINUATION_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/ca/ca126bfaf26be40082ca597e3f70868aeea90098f82bb04aadb2f746dd13c9cb).
- **Genuinely new atom:** Proves a spacetime influence estimate for the hybrid rewiring law rather than repeating the frozen-graph cone.
- **Artifacts:** `fin_post_phy_campaign/DYN-001/REPORT.md`, `fin_post_phy_campaign/DYN-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/DYN-001/results.json`, `fin_post_phy_campaign/DYN-001/replay.py`, `fin_post_phy_campaign/DYN-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 3, 3, 1, 5, 4].

## DYN-002 — Derive coarse mobility and noise from a declared microscopic update

Priority: **P1**, score 34; wave 4.

- **Dependencies:** REF-002, PHA-001, DYN-001.
- **Central question:** Which kinetic coefficients follow once one local update law is actually fixed?
- **Why it matters:** Can replace arbitrary coarse mobility by a theorem while retaining the microscopic PM boundary.
- **Input status:** Same statics/different kinetics is already proved; exact elimination can generate memory.
- **Claim type:** CONDITIONAL_MODEL
- **Exact mathematical object:** A fixed reversible local jump law and its projected phase/amplitude generator.
- **Exact method:** Use generator expansion and averaging on a certified slow manifold. Derive mobility/noise and residual memory; state spectral separation and coherence times. Compare a second attempt-rate law to isolate invariant versus kinetic factors.
- **Acceptance criterion:** A controlled projected generator over a declared time/size window.
- **Failure/refutation criterion:** No time-scale separation, unresolved radial instability or persistent memory invalidates Markov closure.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** An unused attempt-rate law and carrier size.
- **Stop rule:** If memory is non-negligible, retain it; do not improve agreement by fitting a mobility afterward.
- **Computational class:** L
- **What may be concluded:** Coarse coefficients derived conditional on the microscopic update.
- **What must NOT be concluded:** Selection of a physical kinetic category or seconds from equilibrium.
- **Repo de-dup / prior:** [C00:PHY_CAMPAIGN_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/e6/e619e95ca6274dc84872d358dbf5766f3c6b27e099fb1a819043c515351e8149); [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db); [C02:DYNAMICAL_SOURCE_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/17/1745fdf6e0b84ab2927880143ada2ba7cf1520957c6e08e3b57f5e1aae47edb2).
- **Genuinely new atom:** Derives the coarse generator and error from microscopic updates, not another catalogue of supplied flows.
- **Artifacts:** `fin_post_phy_campaign/DYN-002/REPORT.md`, `fin_post_phy_campaign/DYN-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/DYN-002/results.json`, `fin_post_phy_campaign/DYN-002/replay.py`, `fin_post_phy_campaign/DYN-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 5, 3, 3, 2, 5, 5].

## DYN-003 — Derived nonlinear phase potential and drift obstruction

Priority: **P1**, score 33; wave 4.

- **Dependencies:** PHA-002, REF-003.
- **Central question:** Does the controlled FIN phase reduction support a transport mechanism rather than only static walls?
- **Why it matters:** Separates known linear failure from a genuinely nonlinear, microscopically connected question.
- **Input status:** All-orders locking and Umklapp exist; current linear collective class lacks stable localized transport.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** An effective phase action including leading pinning harmonic, higher harmonics and hard-mode corrections.
- **Exact method:** Project the actual potential using PHA bounds. Prove the unbiased-gradient travelling-front energy identity before searching. Identify whether inertia, drive or unequal well values are necessary and label them explicitly.
- **Acceptance criterion:** A controlled nonlinear reduction and a transport/no-drift theorem in its declared dynamics.
- **Failure/refutation criterion:** The proposed autonomous moving structure contradicts energy balance or relies on neglected leading terms.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Unbiased/bias-reversed controls and q values outside the fitted reduction.
- **Stop rule:** Stop gradient travelling-wave searches when the exact balance forbids nonzero speed.
- **Computational class:** L
- **What may be concluded:** A conditional nonlinear transport mechanism or obstruction.
- **What must NOT be concluded:** Fundamental matter, a sourced pump or an already stable soliton.
- **Repo de-dup / prior:** [C00:PHY_CAMPAIGN_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/e6/e619e95ca6274dc84872d358dbf5766f3c6b27e099fb1a819043c515351e8149); [C04:RESONANCE_DIMENSION_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/76/762a8878309629e01aa6d3bb01e527971f27c4d143036e1f49150e56a079ce65); [C05:RESONANCE_CERT_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/6e/6ee7134de2107784707a61503619b7b5b577113afd39d4deaafb5417bc328655); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Derives a periodic phase potential with remainders and an energy-balance transport test; does not assume sine-Gordon or DNLS.
- **Artifacts:** `fin_post_phy_campaign/DYN-003/REPORT.md`, `fin_post_phy_campaign/DYN-003/theorem_or_counterexample.md`, `fin_post_phy_campaign/DYN-003/results.json`, `fin_post_phy_campaign/DYN-003/replay.py`, `fin_post_phy_campaign/DYN-003/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 5, 4, 3, 3, 2, 5, 4].

## DYN-004 — Persistence and orbital stability in the full coupled model

Priority: **P3**, score 25; wave 4.

- **Dependencies:** DYN-002, DYN-003, REF-006.
- **Central question:** Does any reduced moving pattern survive radial modes, lattice pinning and perturbations?
- **Why it matters:** A reduced analytic kink is not yet a FIN excitation.
- **Input status:** Existing nonlinear DNLS/VK and linear recurrence results are not proofs for this newly derived phase model.
- **Claim type:** NUMERICAL_CERTIFICATE_TARGET
- **Exact mathematical object:** A specified travelling/relative-periodic branch, conserved quantities and its linearized operator.
- **Exact method:** Attempt validated continuation, localization and spectral/orbital stability with radiation and pinning controls. Include finite-energy tails and collision-free validity of the dynamic geometry; use the actual conditional kinetic law.
- **Acceptance criterion:** A complete existence/localization/stability certificate or a bounded counterexample to persistence.
- **Failure/refutation criterion:** Only tuned transient copying, uncontrolled radiation or nonrobust pinned states survive.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Two larger meshes, perturbed radial coordinates and a non-fitting initial condition.
- **Stop rule:** One branch only; an X campaign requires separate approval and durable checkpoints.
- **Computational class:** X
- **What may be concluded:** A mathematically stable transported pattern in a specified conditional model.
- **What must NOT be concluded:** Physical particles, universal solitons, immortality or a collision theory.
- **Repo de-dup / prior:** [C00:PHY_CAMPAIGN_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/e6/e619e95ca6274dc84872d358dbf5766f3c6b27e099fb1a819043c515351e8149); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Tests persistence from DYN-003 in the full multi-mode/refined model, rather than repeating older nonlinear kicks.
- **Artifacts:** `fin_post_phy_campaign/DYN-004/REPORT.md`, `fin_post_phy_campaign/DYN-004/theorem_or_counterexample.md`, `fin_post_phy_campaign/DYN-004/results.json`, `fin_post_phy_campaign/DYN-004/replay.py`, `fin_post_phy_campaign/DYN-004/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately; separate bounded campaign/budget approval required.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 4, 2, 5, 3, 4, 3].

## MEM-001 — Type the null before testing memory dimension

Priority: **P0**, score 42; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Which inclusions hold between the actual finite-block, Hankel and positive-process null classes?
- **Why it matters:** A rank-three block factorization is not automatically a three-state all-horizon HMM.
- **Input status:** C09 gives exact rank for explicit families; C13 distinguishes block-latent and single-jump challengers.
- **Claim type:** IDENTIFIABILITY_TARGET
- **Exact mathematical object:** Rank of the featured 4x4 H, rank/nonnegative rank of the raw 9x9 block table, complete path-Hankel rank, state-/edge-emitting HMMs and reset policies.
- **Exact method:** Define stationarity, homogeneity, emission convention and reset independence separately. Distinguish rank(T_F Q T_G^T)<=3 from rank(Q)<=3: the featured null can be strictly larger. Prove inclusions/counterexamples and state which null admits each reused bound.
- **Acceptance criterion:** A typed inclusion diagram and explicit legal scope of every reused rank/sample bound.
- **Failure/refutation criterion:** A claimed hidden-state lower bound uses a challenger outside its declared process class.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Five- and six-jump consistency plus an edge-emitting versus state-emitting control.
- **Stop rule:** No optimizer run until the null class and allowed resets are fixed.
- **Computational class:** S
- **What may be concluded:** Correct finite-record dimension claims and their process limitations.
- **What must NOT be concluded:** Unique hidden ontology or complete rank certification from one block.
- **Repo de-dup / prior:** [C09:HANKEL_MEMORY_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/1f/1f0fa673d18c02f309111b529b0ed0722b8da1b519bb590b1228e11e714e412b); [C10:MEM03_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/15/159582e5baf97850e1afde1f484f2cf5fc9845e44b69a3f32986d3e7670b8222); [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010); [C13:HANKEL03B_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b76d84bfbbe9f7368516c053c6ab6105f2d5ae5e5ffe60fe19b4090f52f5ea96).
- **Genuinely new atom:** Builds certified class-inclusion/separation witnesses for the present design, not another rank-four minor.
- **Artifacts:** `fin_post_phy_campaign/MEM-001/REPORT.md`, `fin_post_phy_campaign/MEM-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-001/results.json`, `fin_post_phy_campaign/MEM-001/replay.py`, `fin_post_phy_campaign/MEM-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 5, 1, 1, 5, 4].

## MEM-002 — Certify the explicit closest-known null witness

Priority: **P1**, score 36; wave 1.

- **Dependencies:** MEM-001.
- **Central question:** What rigorous information lower bound follows from the stored MEM-05 challenger alone?
- **Why it matters:** A single verified close null suffices for a lower bound; global optimality is unnecessary.
- **Input status:** The stored EM witness has numerical KL near 6.10e-6; exact optimum and byte-identical replay are absent.
- **Claim type:** NUMERICAL_CERTIFICATE_TARGET
- **Exact mathematical object:** Rational mixture weights/block laws and interval KL, TV and Hellinger coefficients.
- **Exact method:** Produce exactly normalized nonnegative factors; reconstruct their rank bound algebraically. Enclose the true alternative and logs, then pay rounding error in the 11-symbol lift. Re-derive necessary counts for the precise error criterion.
- **Acceptance criterion:** An explicit feasible null with certified distance and legal minimax lower bound.
- **Failure/refutation criterion:** Rounding destroys feasibility or the distance interval cannot support the quoted lower count.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Independent category-to-symbol reconstruction and a perturbed admissible delta.
- **Stop rule:** Stop after one sufficient witness; do not call it the globally nearest process.
- **Computational class:** M
- **What may be concluded:** A rigorous lower bound for the null class containing that witness.
- **What must NOT be concluded:** Sufficient sample size, global projection optimality or a physical hidden state.
- **Repo de-dup / prior:** [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010); [C13:results/OCB-MEM-05_broader_rank3.json](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/1a/1a0e45b0a80f06a70fe46b75b043b170cb1d432df761b861a01fd55bfe5a8235).
- **Genuinely new atom:** Rationalizes and verifies one admissible witness with outward information bounds, not more EM starts.
- **Artifacts:** `fin_post_phy_campaign/MEM-002/REPORT.md`, `fin_post_phy_campaign/MEM-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-002/results.json`, `fin_post_phy_campaign/MEM-002/replay.py`, `fin_post_phy_campaign/MEM-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 4, 4, 5, 2, 1, 5, 3].

## MEM-003 — Convex-hull obstruction to one-reset universal betting

Priority: **P0**, score 41; wave 1.

- **Dependencies:** MEM-001.
- **Central question:** Can a nontrivial one-reset e-factor be fair under every broad rank-three null?
- **Why it matters:** May rule out a whole statistical architecture before implementation.
- **Input status:** Determinant batching and checkpoint matrix tests exist; unrestricted one-reset betting was not established.
- **Claim type:** NO_GO_TARGET
- **Exact mathematical object:** Nonnegative e(Y) with sup_Q E_Q e<=1 over the specified one-reset null.
- **Exact method:** Use rank-one point masses, or their full-support limits, to compute the convex hull. Test whether fairness forces e<=1 pointwise. Treat pair-matched constraints separately. Distinguish a test supermartingale, an e-process dominated under each fixed null, and a batch e-value.
- **Acceptance criterion:** An exact no-go or nontrivial construction with the correct conditional/reset quantifiers.
- **Failure/refutation criterion:** A claimed universal single-reset bet grows only by violating fairness for an admissible null.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Pair-matched and unrestricted nulls; histories with positive null probability.
- **Stop rule:** If the hull obstruction applies, pivot to batches, confidence sets or extra structure, not a tuned betting fraction.
- **Computational class:** S
- **What may be concluded:** A validity architecture constraint and minimum new statistical structure.
- **What must NOT be concluded:** Impossibility of all anytime tests for a fixed iid composite null.
- **Repo de-dup / prior:** [C11:HANKEL03_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b708714877f0413b12012aea189db1c62baa279941d3ed93eeadec9b0834e918); [C13:HANKEL03B_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b76d84bfbbe9f7368516c053c6ab6105f2d5ae5e5ffe60fe19b4090f52f5ea96).
- **Genuinely new atom:** Tests the geometry of the null and conditional validity, not another determinant e-process.
- **Artifacts:** `fin_post_phy_campaign/MEM-003/REPORT.md`, `fin_post_phy_campaign/MEM-003/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-003/results.json`, `fin_post_phy_campaign/MEM-003/replay.py`, `fin_post_phy_campaign/MEM-003/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 5, 1, 1, 5, 5].

## MEM-004 — All-reset-time variance-adaptive matrix confidence sets

Priority: **P1**, score 34; wave 2.

- **Dependencies:** MEM-001, MEM-003.
- **Central question:** Can direct rank evidence be made valid at every reset time without the worst-case fixed-grid cost?
- **Why it matters:** Extends completed HANKEL-03B to a genuinely stronger inferential object.
- **Input status:** C13 controls a frozen nine-checkpoint grid; exact pair matching is an extra premise.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Time-uniform matrix-mean confidence sets intersected with the rank<=3 variety.
- **Exact method:** Use maximal/stitching or mixture matrix concentration with valid variance bounds. Handle estimated whitening by independent data or joint confidence sets. Prove crossing control at arbitrary reset stopping times; distinguish an e-process from a supermartingale.
- **Acceptance criterion:** Uniform Type-I theorem plus nonasymptotic power bounds under an explicit singular-gap alternative.
- **Failure/refutation criterion:** Variance/whitening reuse invalidates the bound or improvement assumes exact pair matching not supplied.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Nulls near rank boundaries and alternatives outside the frozen period-two family.
- **Stop rule:** Stop if the new bound cannot improve the certified benchmark in a declared regime; preserve validity result without a practicality claim.
- **Computational class:** L
- **What may be concluded:** An all-reset-time composite test with stated power/cost.
- **What must NOT be concluded:** Physical-time validity or universal optimal rank discrimination.
- **Repo de-dup / prior:** [C13:HANKEL03B_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b76d84bfbbe9f7368516c053c6ab6105f2d5ae5e5ffe60fe19b4090f52f5ea96); [C10:HANKEL02_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/f2/f22d78402f46ddb55cd75e98626be7bbce8cb85d57a2948ef116ccb965d72aec); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** New all-time/unknown-variance control and power analysis, not relabeling finite checkpoints as a continuous e-process.
- **Artifacts:** `fin_post_phy_campaign/MEM-004/REPORT.md`, `fin_post_phy_campaign/MEM-004/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-004/results.json`, `fin_post_phy_campaign/MEM-004/replay.py`, `fin_post_phy_campaign/MEM-004/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 4, 4, 4, 3, 2, 5, 4].

## MEM-005 — Safe likelihood denominators without an exact optimizer

Priority: **P1**, score 38; wave 2.

- **Dependencies:** MEM-001, MEM-002, MEM-003.
- **Central question:** Can certified upper relaxations of null likelihood make a practical split/batch e-test?
- **Why it matters:** The denominator needs an upper bound, not necessarily an attained global MLE.
- **Input status:** C13 correctly rejects local EM maxima as denominator certificates.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** L_alt(D_test)/U_null(D_test), where U_null >= sup_Q L_Q on the chosen null, including the featured-rank null if that is the tested hypothesis.
- **Exact method:** Learn numerator on independent/past data. Build certified outer relaxations using simplex constraints, rank-distance information and bounded branch/SDP refinement. Optimizing only raw/nonnegative-rank-three tables does not upper-bound the larger featured-rank null. Include zero cells correctly and prove batch/sequential validity.
- **Acceptance criterion:** A valid nontrivial likelihood e-value with a rigorous denominator bound and held-out power estimate.
- **Failure/refutation criterion:** A lower/local likelihood replaces an upper bound, or a claimed gain vanishes under the necessary relaxation.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Adversarial null tables, sparse observed counts and a second alias family.
- **Stop rule:** Bound global-search effort; a valid but powerless unconstrained bound is a stop result, not a practical solution.
- **Computational class:** L
- **What may be concluded:** Conditional finite-sample testing without assuming exact MLE availability.
- **What must NOT be concluded:** Global EM optimality, arbitrary online data reuse or numerical e-values as proofs.
- **Repo de-dup / prior:** [C13:HANKEL03B_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b76d84bfbbe9f7368516c053c6ab6105f2d5ae5e5ffe60fe19b4090f52f5ea96); [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Computes safe data-dependent likelihood upper bounds for the declared class; no extra EM optimism.
- **Artifacts:** `fin_post_phy_campaign/MEM-005/REPORT.md`, `fin_post_phy_campaign/MEM-005/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-005/results.json`, `fin_post_phy_campaign/MEM-005/replay.py`, `fin_post_phy_campaign/MEM-005/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 4, 3, 1, 5, 5].

## MEM-006 — Positive process extension of finite block models

Priority: **P1**, score 38; wave 2.

- **Dependencies:** MEM-001, MEM-002.
- **Central question:** Which stored block-null tables extend to the intended positive all-horizon process with a bounded state count?
- **Why it matters:** Links statistical rank witnesses to a legitimate dynamical realization class.
- **Input status:** Finite block rank and all-horizon realization are explicitly separated in C13.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** Positive state-/edge-emitting automata extending prescribed four-jump probabilities.
- **Exact method:** Form exact normalization and consistency equations across cuts; construct an extension or a separating witness. Distinguish time-inhomogeneous reset models from stationary homogeneous HMMs. Verify more than the central 2|2 factorization.
- **Acceptance criterion:** An explicit positive realization with a proved state bound, or infeasibility in a specified class.
- **Failure/refutation criterion:** The central rank certificate cannot be extended with the claimed hidden-state semantics.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Unused longer words and a second reset distribution.
- **Stop rule:** No global ontology inference from a successful finite extension; stop at the declared state-budget cap.
- **Computational class:** L
- **What may be concluded:** Legitimate class membership or a positive-realization obstruction.
- **What must NOT be concluded:** Equality of ordinary rank and nonnegative hidden-state dimension in general.
- **Repo de-dup / prior:** [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010); [C09:HANKEL_MEMORY_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/1f/1f0fa673d18c02f309111b529b0ed0722b8da1b519bb590b1228e11e714e412b); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Pays positivity and cross-horizon consistency, beyond ordinary linear realization or one factorized table.
- **Artifacts:** `fin_post_phy_campaign/MEM-006/REPORT.md`, `fin_post_phy_campaign/MEM-006/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-006/results.json`, `fin_post_phy_campaign/MEM-006/replay.py`, `fin_post_phy_campaign/MEM-006/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 4, 3, 1, 5, 5].

## MEM-007 — Tight horizon obstruction under explicit reset and stationarity rules

Priority: **P1**, score 39; wave 1.

- **Dependencies:** MEM-001.
- **Central question:** How much of the log(H)-to-H gap survives after the hidden-process assumptions are fixed?
- **Why it matters:** Avoids paying for the wrong adversarial memory class.
- **Input status:** C10 proves a logarithmic family lower bound and linear generic witness upper bound, not a query minimax theorem.
- **Claim type:** NO_GO_TARGET
- **Exact mathematical object:** H-state hidden processes constrained by a declared initialization/stationarity class.
- **Exact method:** Test a delayed H-state reset clock for the unrestricted class; it must not be used as a stationary counterexample. For the stationary class seek a positive adversarial construction or sharper upper bound. State any mixing/gap restrictions.
- **Acceptance criterion:** A corrected class-dependent horizon lower/upper bound with explicit process witnesses.
- **Failure/refutation criterion:** The previous claimed lower bound relies on an initialization forbidden by the proposed experiment.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Unknown versus controlled hidden reset and an unseen H.
- **Stop rule:** Stop when assumptions change the class; no exponential word enumeration presented as a universal query lower bound.
- **Computational class:** M
- **What may be concluded:** Scoped horizon complexity, not a physical memory capacity.
- **What must NOT be concluded:** Uniform exclusion of unbounded memory by a finite record.
- **Repo de-dup / prior:** [C10:MEM03_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/15/159582e5baf97850e1afde1f484f2cf5fc9845e44b69a3f32986d3e7670b8222); [C08:MEMORY_ALIAS_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b2/b28c3690e427088b42765e32a6385031db8c616cd27841d6e4582e77d310a287); [C09:HANKEL_MEMORY_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/1f/1f0fa673d18c02f309111b529b0ed0722b8da1b519bb590b1228e11e714e412b).
- **Genuinely new atom:** Separates resettable nonstationary clocks from stationary ergodic aliases and seeks a sharper worst-case theorem.
- **Artifacts:** `fin_post_phy_campaign/MEM-007/REPORT.md`, `fin_post_phy_campaign/MEM-007/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-007/results.json`, `fin_post_phy_campaign/MEM-007/replay.py`, `fin_post_phy_campaign/MEM-007/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 5, 2, 1, 5, 4].

## MEM-008 — Cost-aware adaptive path experiment design

Priority: **P1**, score 33; wave 3.

- **Dependencies:** MEM-004, MEM-005, MEM-006, MEM-007.
- **Central question:** Which resets, path lengths and contrasts maximize worst-case information per admitted cost?
- **Why it matters:** Memory detection and full dimension certification need different experiments.
- **Input status:** C09/C10 already separate easy lag detection from expensive full-rank certification.
- **Claim type:** OPERATIONAL_PROTOCOL
- **Exact mathematical object:** Adaptive controls over reset distributions, horizons and observable contrasts with fixed reset/jump costs.
- **Exact method:** Specify whether the adversary shares parameters across actions. Derive a minimax information criterion and valid predictable allocation. Check loss from three-category compression against hidden alternatives inside discarded categories.
- **Acceptance criterion:** A frozen policy with valid error control and a certified or bounded information advantage over a fixed four-jump design.
- **Failure/refutation criterion:** Apparent advantage relies on known alias structure or unpriced long paths/resets.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** A non-periodic HMM, detector-memory control and a withheld within-category perturbation.
- **Stop rule:** No adaptive feature search on the scoring data without a proved selection correction.
- **Computational class:** L
- **What may be concluded:** A conditional experimental design and its cost/error tradeoff.
- **What must NOT be concluded:** A universal optimal protocol or laboratory feasibility without apparatus costs.
- **Repo de-dup / prior:** [C09:ADAPTIVE_MEMORY_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/a8/a83076879d8a06ccf7c15dfe738b90e8fa4e54cec592502085e4d1d60eb7c6f5); [C10:MEM03_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/15/159582e5baf97850e1afde1f484f2cf5fc9845e44b69a3f32986d3e7670b8222); [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010).
- **Genuinely new atom:** Optimizes a precise action/query model against shared composite alternatives, not more exhaustive words.
- **Artifacts:** `fin_post_phy_campaign/MEM-008/REPORT.md`, `fin_post_phy_campaign/MEM-008/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-008/results.json`, `fin_post_phy_campaign/MEM-008/replay.py`, `fin_post_phy_campaign/MEM-008/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 4, 4, 3, 3, 2, 5, 4].

## MEM-009 — Honest lower/upper power frontier

Priority: **P1**, score 39; wave 3.

- **Dependencies:** MEM-002, MEM-004, MEM-005.
- **Central question:** What sample-cost bracket is actually certified for broad rank exclusion?
- **Why it matters:** Closes the semantic gap between necessary resets, population design points and sufficient power.
- **Input status:** C13 provides a close numerical null and Type-I-valid checkpoints but no minimax-optimal test.
- **Claim type:** IDENTIFIABILITY_TARGET
- **Exact mathematical object:** A testing risk surface indexed by separation, reset cost, nuisance and null class.
- **Exact method:** Use MEM-002 lower bounds and proved power of MEM-004/005. For a matrix test pay both the null threshold and alternative deviation. Prove no uniform finite power as separation tends to zero; distinguish Bayes from two-sided error criteria.
- **Acceptance criterion:** A reproducible valid bracket or an explicit gap too wide for operational use.
- **Failure/refutation criterion:** Claims of practicality rest only on r<population sigma or expected log evidence crossing a threshold.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** At least two deltas and a non-frozen rank-four alternative.
- **Stop rule:** If the required cost exceeds the admitted platform budget, stop the rank assay rather than weaken its null silently.
- **Computational class:** M
- **What may be concluded:** An error-criterion-matched cost frontier.
- **What must NOT be concluded:** Measured power or a globally optimal minimax constant without proof.
- **Repo de-dup / prior:** [C13:HANKEL03B_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b76d84bfbbe9f7368516c053c6ab6105f2d5ae5e5ffe60fe19b4090f52f5ea96); [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010).
- **Genuinely new atom:** Combines certified necessary and sufficient counts for the same hypotheses/errors; does not compare incompatible nulls.
- **Artifacts:** `fin_post_phy_campaign/MEM-009/REPORT.md`, `fin_post_phy_campaign/MEM-009/theorem_or_counterexample.md`, `fin_post_phy_campaign/MEM-009/results.json`, `fin_post_phy_campaign/MEM-009/replay.py`, `fin_post_phy_campaign/MEM-009/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 5, 2, 1, 5, 4].

## OP-001 — Wave response as a complete probability experiment

Priority: **P0**, score 42; wave 1.

- **Dependencies:** AUD-001.
- **Central question:** Which complete observation model licenses treating squared wave displacement as event probabilities?
- **Why it matters:** The existing Chernoff attempt bound needs a probability/readout premise, not only a response amplitude.
- **Input status:** CHANNEL-03/05 use zero-velocity squared wave response and retain every attempt; apparatus interpretation remains conditional.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** C(t)=cos(t sqrt(A)), its off-diagonal squares, return channel and missing-norm outcome.
- **Exact method:** Check subprobability normalization for PSD A. Construct the complementary outcome or conservative doubled-state dilation and distinguish raw intensity, normalized displacement and a calibrated Bernoulli readout. Re-derive the frozen bounds for the chosen record convention.
- **Acceptance criterion:** One complete probability model with all outcomes and readout premises explicit, or a protocol-scope correction.
- **Failure/refutation criterion:** The proposed event probabilities require undocumented renormalization or are not the outcomes used by the concentration bound.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** A two-eigenvalue control and a second initial preparation.
- **Stop rule:** Do not claim a laboratory probability rule from mathematical norm completion.
- **Computational class:** S
- **What may be concluded:** A mathematically consistent conditional wave-readout experiment.
- **What must NOT be concluded:** Born-rule derivation, quantum dynamics or empirical detector behavior.
- **Repo de-dup / prior:** [C10:CHANNEL03_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/db/db2f9e28531ea8ed04c32e4234d40a17158cd95194437635598b499ad81849a7); [C11:CHANNEL04_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/d2/d2289a8c5519c6ebeadc6d163612be4b15ebb441d5421fbb8eb15f27a631ee35); [C13:CHANNEL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/cf/cfab131c8cca3a29ca26cb45b1451974d08a8e8e1fa7dd799eb20719aab7af83).
- **Genuinely new atom:** Pays the complete-outcome/dilation contract; does not repeat leading time powers.
- **Artifacts:** `fin_post_phy_campaign/OP-001/REPORT.md`, `fin_post_phy_campaign/OP-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/OP-001/results.json`, `fin_post_phy_campaign/OP-001/replay.py`, `fin_post_phy_campaign/OP-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 5, 1, 1, 5, 4].

## OP-002 — Clock-map-free escape–profile geometry

Priority: **P0**, score 42; wave 2.

- **Dependencies:** OP-001.
- **Central question:** Can destination shape as a function of escape distinguish channels without a locally linear control clock?
- **Why it matters:** Potentially removes one strong calibration premise of CHANNEL-04/05.
- **Input status:** Leading destination profiles and escape powers are known; existing finite protocols assume bounded control-ratio distortion.
- **Claim type:** THEOREM_TARGET
- **Exact mathematical object:** S=sum_(j!=i) p_ji and r_j=p_ji/S, viewed as an unparameterized curve.
- **Exact method:** Derive strict-operator expansions. Candidate distinction: unitary r-r0=O(S), while zero-velocity wave r-r0=B sqrt(S)+O(S), with B nonzero unless an exceptional spectral relation holds. Prove the exception classification, finite-S remainder, calibration sensitivity and sampling cost. Heat has a different leading profile for unequal weights.
- **Acceptance criterion:** A finite observable-window separation invariant under admissible monotone time reparameterization, or an exact curve-equivalence counterexample.
- **Failure/refutation criterion:** The coefficient vanishes in the claimed class, curves are observationally equivalent, or calibration overwhelms the effect.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Nonlinear control clocks, unseen spectral perturbations and the complete-graph exceptional control.
- **Stop rule:** If asymptotic separation has no useful finite calibrated window, report a theorem-only outcome; no practicality claim.
- **Computational class:** M
- **What may be concluded:** Conditional clock-map-free channel discrimination from joint records.
- **What must NOT be concluded:** A physical clock, physical wave ontology, zero-calibration inference or a smaller sample count without proof.
- **Repo de-dup / prior:** [C10:CHANNEL03_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/db/db2f9e28531ea8ed04c32e4234d40a17158cd95194437635598b499ad81849a7); [C11:CHANNEL04_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/d2/d2289a8c5519c6ebeadc6d163612be4b15ebb441d5421fbb8eb15f27a631ee35); [C13:CHANNEL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/cf/cfab131c8cca3a29ca26cb45b1451974d08a8e8e1fa7dd799eb20719aab7af83).
- **Genuinely new atom:** Eliminates the latent time parameter between two observables; no repeat of slopes at proportional controls.
- **Artifacts:** `fin_post_phy_campaign/OP-002/REPORT.md`, `fin_post_phy_campaign/OP-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/OP-002/results.json`, `fin_post_phy_campaign/OP-002/replay.py`, `fin_post_phy_campaign/OP-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 5, 2, 1, 5, 5].

## OP-003 — Unknown carrier versus structured detector support

Priority: **P1**, score 39; wave 2.

- **Dependencies:** SRC-003, AUD-001.
- **Central question:** Can carrier topology and calibration structure be identified jointly without assuming the cycle being tested?
- **Why it matters:** The CAL-05 neighbor-confusion premise can become circular if inferred from the same carrier claim.
- **Input status:** Known inverse calibration and support-preserving error factors are already proved.
- **Claim type:** IDENTIFIABILITY_TARGET
- **Exact mathematical object:** Unlabeled transition laws, preparation/effect maps and unknown support-preserving confusion patterns.
- **Exact method:** Separate independently calibrated support from graph hypotheses. Search exact alternative graph/detector factorizations. Determine the minimum intervention or calibration record that breaks them, keeping label permutations as gauge.
- **Acceptance criterion:** A joint identifiability theorem or an explicit calibration-only topological alias.
- **Failure/refutation criterion:** The support assumption presupposes the carrier conclusion or independent calibration cannot distinguish alternative graphs.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** A noncyclic carrier and relabeled apparatus channels.
- **Stop rule:** No same-data calibration that selects the desired cycle after inspection.
- **Computational class:** M
- **What may be concluded:** A noncircular conditional carrier protocol or a necessary apparatus obligation.
- **What must NOT be concluded:** Physical adjacency solely from largest observed rates.
- **Repo de-dup / prior:** [C08:CALIBRATION_CHANNEL_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c8/c859dfa8b8f8a35a7034ad1db0e991b2a897aadeda26e5d48255c95755d4ea45); [C13:CAL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/2e/2e8abe1af9ec266162a985ddcfe5594b864d322a8425115f86de05709e24f9f4); [C07:OPERATIONAL_CARRIER_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/09/0980bfba87d5ed26304dcdc5b0c0a5c4c7a6dfbe954a609d70cb056dcf04437f).
- **Genuinely new atom:** Tests joint graph/calibration equivalence classes, not a further numerical radius improvement with a known cycle.
- **Artifacts:** `fin_post_phy_campaign/OP-003/REPORT.md`, `fin_post_phy_campaign/OP-003/theorem_or_counterexample.md`, `fin_post_phy_campaign/OP-003/results.json`, `fin_post_phy_campaign/OP-003/replay.py`, `fin_post_phy_campaign/OP-003/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 4, 4, 2, 1, 5, 5].

## OP-004 — Detector memory and postselection rank aliases

Priority: **P0**, score 41; wave 3.

- **Dependencies:** MEM-001, MEM-006, OP-003.
- **Central question:** How much observed predictive rank can be attributed to the detector rather than FIN dynamics?
- **Why it matters:** A system-memory claim is invalid if readout memory or conditioning can supply the same rank.
- **Input status:** C08/C13 calibration is primarily memoryless; memory aliases of the visible process are known.
- **Claim type:** IDENTIFIABILITY_TARGET
- **Exact mathematical object:** A system process composed with a bounded detector transducer, including no-click records.
- **Exact method:** Prove rank nonincrease/invariance under specified memoryless linear readout maps. Track what changes under finite-state detector memory and length-dependent postselection normalization. Bound system/detector ambiguity and test reset failures; do not assume iid blocks when the detector state persists.
- **Acceptance criterion:** A class-dependent rank attribution bound or explicit indistinguishable system/detector allocations.
- **Failure/refutation criterion:** The observed extra rank is fully explained by admitted calibration dynamics or conditioning.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** An iid system with a two-state detector and a memoryful system with memoryless readout.
- **Stop rule:** If detector memory is unbounded, stop unique system-rank attribution and export the required independent reset/calibration gate.
- **Computational class:** M
- **What may be concluded:** Predictive-dimension statements robust to a declared detector class.
- **What must NOT be concluded:** Microscopic system ontology from visible Hankel rank.
- **Repo de-dup / prior:** [C08:CALIBRATION_CHANNEL_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c8/c859dfa8b8f8a35a7034ad1db0e991b2a897aadeda26e5d48255c95755d4ea45); [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010); [C13:CAL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/2e/2e8abe1af9ec266162a985ddcfe5594b864d322a8425115f86de05709e24f9f4); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Joins temporal calibration with positive realization and raw-attempt bookkeeping.
- **Artifacts:** `fin_post_phy_campaign/OP-004/REPORT.md`, `fin_post_phy_campaign/OP-004/theorem_or_counterexample.md`, `fin_post_phy_campaign/OP-004/results.json`, `fin_post_phy_campaign/OP-004/replay.py`, `fin_post_phy_campaign/OP-004/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 4, 2, 1, 5, 5].

## OP-005 — One joint error ledger for carrier, channel and memory

Priority: **P1**, score 33; wave 3.

- **Dependencies:** MEM-008, MEM-009, OP-002, OP-004.
- **Central question:** Can the successive operational claims share data without invalidating their error guarantees?
- **Why it matters:** Separate valid protocols do not automatically compose into a valid adaptive discovery pipeline.
- **Input status:** C13 has separate rank, calibration and channel designs with different assumptions.
- **Claim type:** OPERATIONAL_PROTOCOL
- **Exact mathematical object:** A frozen protocol graph with shared nuisance variables, filtration and error budget.
- **Exact method:** Specify training/calibration/scoring splits or valid predictable reuse. Couple the nuisance sets for topology, readout and rank. Carry conditional failure probabilities and abstention through the graph; reject multiplication of dependent e-values without a theorem.
- **Acceptance criterion:** An end-to-end error guarantee with every attempted/no-click/reset record retained.
- **Failure/refutation criterion:** Data-dependent selection or reused calibration makes a claimed conditional guarantee false.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Adversarial ordering, stopping and calibration drift within the declared model.
- **Stop rule:** No empirical execution until the joint ledger and external calibration pass.
- **Computational class:** M
- **What may be concluded:** A coherent conditional operational workflow.
- **What must NOT be concluded:** Independent custody, actual calibration or experimental success from a software schema.
- **Repo de-dup / prior:** [C13:HANKEL03B_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/b7/b76d84bfbbe9f7368516c053c6ab6105f2d5ae5e5ffe60fe19b4090f52f5ea96); [C13:CAL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/2e/2e8abe1af9ec266162a985ddcfe5594b864d322a8425115f86de05709e24f9f4); [C13:CHANNEL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/cf/cfab131c8cca3a29ca26cb45b1451974d08a8e8e1fa7dd799eb20719aab7af83); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Proves sequential composition across calibration, feature selection and claims, beyond a two-budget union bound.
- **Artifacts:** `fin_post_phy_campaign/OP-005/REPORT.md`, `fin_post_phy_campaign/OP-005/theorem_or_counterexample.md`, `fin_post_phy_campaign/OP-005/results.json`, `fin_post_phy_campaign/OP-005/replay.py`, `fin_post_phy_campaign/OP-005/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 4, 4, 4, 2, 2, 5, 4].

## EMP-001 — Platform readiness and information-budget gate

Priority: **P1**, score 34; wave 3.

- **Dependencies:** OP-005.
- **Central question:** Is any supplied real platform capable of the required preparations, resets and calibrated records?
- **Why it matters:** Prevents billion-reset mathematical designs becoming unsupported experimental recommendations.
- **Input status:** EXTERNAL-OCB-01 is already gated; no apparatus or raw physical records are supplied here.
- **Claim type:** EMPIRICAL_GATE
- **Exact mathematical object:** A platform capability/calibration/custody record and an admissible experiment budget.
- **Exact method:** Request platform facts only when external execution is authorized. Match reset fidelity, detector memory, efficiencies, control ranges and costs to OP-005/MEM-009. Without real inputs emit NOT_READY with exact missing records.
- **Acceptance criterion:** Independent records satisfy all required gates, or a precise blocked/readiness verdict is issued.
- **Failure/refutation criterion:** Calibration is circular, resets fail the model, or guaranteed information exceeds the actual budget.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** An independent calibration session and held-out negative-control run.
- **Stop rule:** No purchases, contacts, uploads or simulated replacement for missing nature data.
- **Computational class:** S
- **What may be concluded:** Readiness or nonreadiness of a named conditional protocol.
- **What must NOT be concluded:** Empirical validation or FIN ontology from apparatus capability.
- **Repo de-dup / prior:** [C13:CHANNEL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/cf/cfab131c8cca3a29ca26cb45b1451974d08a8e8e1fa7dd799eb20719aab7af83); [C13:NEXT_ATOMS.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/ad/ad3e78c61ca4ce5d4540bf9f3c922968cdeaa7aef8363b0b07d1a870e8590322); [fin_physics_review/MASTER_INTAKE_20260923.md](fin_physics_review/MASTER_INTAKE_20260923.md).
- **Genuinely new atom:** Evaluates readiness for the new joint protocol and its certified cost frontier; does not reannounce an experiment.
- **Artifacts:** `fin_post_phy_campaign/EMP-001/REPORT.md`, `fin_post_phy_campaign/EMP-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/EMP-001/results.json`, `fin_post_phy_campaign/EMP-001/replay.py`, `fin_post_phy_campaign/EMP-001/MANIFEST.sha256`.
- **Activation:** External data/coordination not authorized; readiness may be assessed symbolically.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 5, 5, 2, 2, 2, 4, 4].

## EMP-002 — Preregister a discriminating pilot, not a full rank campaign by default

Priority: **P1**, score 38; wave 3.

- **Dependencies:** EMP-001.
- **Central question:** Which lowest-cost real observation would eliminate the largest admitted model class?
- **Why it matters:** Carrier confirmation, channel discrimination and full memory rank have very different value/cost.
- **Input status:** A platform may be admitted only through EMP-001; otherwise this remains a parameterized design.
- **Claim type:** OPERATIONAL_PROTOCOL
- **Exact mathematical object:** One preregistered observable/control allocation and competing system-plus-detector models.
- **Exact method:** Choose carrier, escape–profile or memory assay by the proved cost frontier. Freeze effect size, nuisance bounds, abstention and maximum cost. Include the closest admitted non-FIN or detector-only alternative and a held-out prediction.
- **Acceptance criterion:** A no-fit pilot with feasible error/power bounds for the admitted platform; else an explicit no-feasible-pilot verdict.
- **Failure/refutation criterion:** Every proposed observation is observationally aliased or requires post-hoc fit freedom.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Independent preparations and scoring data not used for calibration or protocol selection.
- **Stop rule:** Stop or redesign under a new protocol ID after failure; do not weaken the null while keeping the old success label.
- **Computational class:** M
- **What may be concluded:** A preregistered future experiment; EV only after separately authorized physical data.
- **What must NOT be concluded:** A foundation-of-physics confirmation from a synthetic pilot.
- **Repo de-dup / prior:** [C07:OPERATIONAL_CARRIER_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/09/0980bfba87d5ed26304dcdc5b0c0a5c4c7a6dfbe954a609d70cb056dcf04437f); [C13:CHANNEL05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/cf/cfab131c8cca3a29ca26cb45b1451974d08a8e8e1fa7dd799eb20719aab7af83); [C13:MEM05_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c2/c22f4f0dfd0f67c8890a86b869ad64d266193a6154663df5b30886c11151f010).
- **Genuinely new atom:** Chooses among the new joint observable options using certified information, rather than automatically executing the old frozen count.
- **Artifacts:** `fin_post_phy_campaign/EMP-002/REPORT.md`, `fin_post_phy_campaign/EMP-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/EMP-002/results.json`, `fin_post_phy_campaign/EMP-002/replay.py`, `fin_post_phy_campaign/EMP-002/MANIFEST.sha256`.
- **Activation:** External data/coordination not authorized; readiness may be assessed symbolically.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [5, 5, 5, 2, 2, 1, 5, 4].

## GATE-001 — New operational structure before reopening quantum claims

Priority: **P3**, score 22; wave 5.

- **Dependencies:** SRC-003, OP-004.
- **Central question:** Has a new observation/composition object actually escaped the previously refuted classical or ensemble-inconsistent class?
- **Why it matters:** Prevents memory complexity and complex notation being promoted to QM.
- **Input status:** PHY-025/026 and older ensemble-affinity/feedback obstructions remain in force.
- **Claim type:** SPECULATIVE_ONLY
- **Exact mathematical object:** A newly supplied state/effect/composition algebra with operational independence assumptions.
- **Exact method:** Only if SRC-003/OP-004 exports new structure, test affine preparation consistency, classical simulation with the admitted memory class, noncommuting instruments and a sourced composite rule. Unlimited classical memory cannot be excluded by a finite transcript without restrictions.
- **Acceptance criterion:** A precisely defined new quantum-bridge problem passing all prerequisite gates, or STOP.
- **Failure/refutation criterion:** The result is imported QM, classical memory or an ensemble-dependent law violating its claimed operational assumptions.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Equivalent mixed preparations and independent-setting controls.
- **Stop rule:** Default is STOP—DO NOT PROMOTE; no Standard Model or quantum derivation task is authorized.
- **Computational class:** S
- **What may be concluded:** Admission or refusal of a later quantum research lane.
- **What must NOT be concluded:** QM, Born rule or physical nonlocality from the current core.
- **Repo de-dup / prior:** [C00:PHY_CAMPAIGN_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/e6/e619e95ca6274dc84872d358dbf5766f3c6b27e099fb1a819043c515351e8149); [C05:STATE_DYN_PROX_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/c5/c588238d4f51da564e091f141c8119d52240c3e58b08041c5cceccb88f7eb398); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Dormant admission gate for genuinely new outcomes, not another execution of the old CHSH/simplex test.
- **Artifacts:** `fin_post_phy_campaign/GATE-001/REPORT.md`, `fin_post_phy_campaign/GATE-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/GATE-001/results.json`, `fin_post_phy_campaign/GATE-001/replay.py`, `fin_post_phy_campaign/GATE-001/MANIFEST.sha256`.
- **Activation:** DORMANT; require a genuinely new typed object; otherwise STOP—DO NOT PROMOTE.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [3, 4, 4, 1, 2, 4, 3, 2].

## GATE-002 — Common action and universal-probe gate before gravity

Priority: **P3**, score 24; wave 5.

- **Dependencies:** REF-006, GEO-003, DYN-004.
- **Central question:** Does the new dynamic graph support reciprocal backreaction and a common operational geometry?
- **Why it matters:** State-dependent distance alone is not gravity.
- **Input status:** The prior gravity gate is closed; new DYN/REF results would be additional conditional objects.
- **Claim type:** SPECULATIVE_ONLY
- **Exact mathematical object:** A joint field/graph action, its conservation identities and at least two probe responses.
- **Exact method:** Require one variational/balance source for both matter and graph response. Test reciprocity, conservation, causal propagation and probe universality under the same calibration. Expose independent kinetic ratios rather than identifying all speeds with c.
- **Acceptance criterion:** All necessary prerequisites are explicitly met, or a typed missing-atom/STOP result.
- **Failure/refutation criterion:** No backreaction, arbitrary probe-specific metrics or inconsistent energy exchange.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Two distinct probe species and one withheld refinement.
- **Stop rule:** No gravity language unless all gates pass; no GR coefficient fitting or dimensional G target.
- **Computational class:** S
- **What may be concluded:** Admission or refusal of a later gravity bridge.
- **What must NOT be concluded:** GR, equivalence principle or spacetime from graph deformation.
- **Repo de-dup / prior:** [C00:PHY_CAMPAIGN_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/e6/e619e95ca6274dc84872d358dbf5766f3c6b27e099fb1a819043c515351e8149); [C01:SOURCE_FOLLOWUP_REPORT.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/92/9297623bb4fc7c1ce14e2a3596050506e4b8c142b52de0505242cb4b0d3465db); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Dormant consistency test for a genuinely new coupled action, not re-fitting Einstein equations.
- **Artifacts:** `fin_post_phy_campaign/GATE-002/REPORT.md`, `fin_post_phy_campaign/GATE-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/GATE-002/results.json`, `fin_post_phy_campaign/GATE-002/replay.py`, `fin_post_phy_campaign/GATE-002/MANIFEST.sha256`.
- **Activation:** DORMANT; require a genuinely new typed object; otherwise STOP—DO NOT PROMOTE.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [3, 4, 5, 1, 2, 4, 3, 2].

## SYN-001 — Adversarial proof and scope consolidation

Priority: **P1**, score 35; wave 6.

- **Dependencies:** AUD-001, SRC-001, SRC-002, SRC-003, SRC-004, SRC-005, REF-001, REF-002, REF-003, REF-004, REF-005, REF-006, PHA-001, PHA-002, PHA-003, GEO-001, GEO-002, GEO-003, DYN-001, DYN-002, DYN-003, DYN-004, MEM-001, MEM-002, MEM-003, MEM-004, MEM-005, MEM-006, MEM-007, MEM-008, MEM-009, OP-001, OP-002, OP-003, OP-004, OP-005, EMP-001, EMP-002, GATE-001, GATE-002.
- **Central question:** Which newly completed claims survive independent challenge in their actual domains?
- **Why it matters:** A large task queue needs a fail-closed integration boundary.
- **Input status:** All attempted task artifacts, including failures and gated skips; existing master sources remain immutable.
- **Claim type:** NUMERICAL_CERTIFICATE_TARGET
- **Exact mathematical object:** New proof DAG, exact inputs, numerical enclosures and claim ladders.
- **Exact method:** Check hashes, signs, null-class inclusion, lower-versus-upper bound direction and no-go scope. Challenge new results with their declared countermodels. Label same-code versus independent implementation replay; do not launch every old producer again.
- **Acceptance criterion:** Each promoted claim has a complete reproducible dependency path and no unrecorded PM or numerical substitution.
- **Failure/refutation criterion:** A missing premise or invalid certificate remains on the promotion path.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** Corrupted boundaries, changed detector class, stale inputs and deliberate physical-label mutations.
- **Stop rule:** Quarantine only dependent claims and identify the minimal repair; no fabricated global PASS.
- **Computational class:** M
- **What may be concluded:** A review-ready set of scoped new results.
- **What must NOT be concluded:** Physical evidence, independent custody or blanket acceptance of execution-complete tasks.
- **Repo de-dup / prior:** [fin_physics_review/MASTER_INTAKE_20260923.md](fin_physics_review/MASTER_INTAKE_20260923.md); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Audits new certificates and their changed premises, not another master archive deduplication.
- **Artifacts:** `fin_post_phy_campaign/SYN-001/REPORT.md`, `fin_post_phy_campaign/SYN-001/theorem_or_counterexample.md`, `fin_post_phy_campaign/SYN-001/results.json`, `fin_post_phy_campaign/SYN-001/replay.py`, `fin_post_phy_campaign/SYN-001/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 5, 4, 5, 2, 2, 5, 3].

## SYN-002 — Re-optimize the programme after outcomes, then hand off

Priority: **P1**, score 39; wave 6.

- **Dependencies:** SYN-001.
- **Central question:** Which next atom now has highest information gain, including the possibility of stopping?
- **Why it matters:** Avoids sunk-cost continuation and endless boundary polishing.
- **Input status:** SYN-001 and all attempted/gated task dispositions; score sensitivity and kill ledger.
- **Claim type:** OPERATIONAL_PROTOCOL
- **Exact mathematical object:** A versioned source-obligation graph and decision ledger.
- **Exact method:** Update the eight priority dimensions, document demotions and killed branches, and recommend at most three new atoms. Export portable inputs, scripts, theorem scopes and an AGENTS proposal only; preserve physical gates.
- **Acceptance criterion:** A supervisor can distinguish new theorems, conditional models, numerical evidence and external blockers without the conversation.
- **Failure/refutation criterion:** The handoff hides failed directions or recommends the same closed search without a new object.
- **Pretarget information:** Not a source-selection task; premises remain explicitly supplied.
- **Held-out test:** A reviewer reconstructs the first next decision using only the handoff.
- **Stop rule:** Stop after handoff; leave no background workers or automatic downstream campaign.
- **Computational class:** S
- **What may be concluded:** A justified next decision or a scoped research stop.
- **What must NOT be concluded:** Completion of a physical theory from terminal task dispositions.
- **Repo de-dup / prior:** [C13:NEXT_ATOMS.md](FIN_PHY_MASTER_HANDOFF_20260923/store/sha256/ad/ad3e78c61ca4ce5d4540bf9f3c922968cdeaa7aef8363b0b07d1a870e8590322); [fin_physics_review/MASTER_INTAKE_20260923.md](fin_physics_review/MASTER_INTAKE_20260923.md); [AGENTS.md](AGENTS.md).
- **Genuinely new atom:** Uses genuinely new outcomes to change the research order; the present source and OCB queues are not reopened by name.
- **Artifacts:** `fin_post_phy_campaign/SYN-002/REPORT.md`, `fin_post_phy_campaign/SYN-002/theorem_or_counterexample.md`, `fin_post_phy_campaign/SYN-002/results.json`, `fin_post_phy_campaign/SYN-002/replay.py`, `fin_post_phy_campaign/SYN-002/MANIFEST.sha256`.
- **Activation:** PROPOSED_ONLY; authorize execution separately.
- **Scores (information, falsification, foundation, independence, cost penalty, duplication penalty, transfer, premise reduction):** [4, 5, 5, 5, 1, 2, 5, 4].

