# FIN: bounded generative continuation after the methodology reconciliation

Date: 2026-10-14
Repository baseline: c24773ec, plus the current working-tree additions to AGENTS.md and fin_post_326_methodology_review/.
Status: RESEARCH_PLAN. The tasks below have not been executed by this document.
Main campaign: RG-C01–RG-C06. Architectural decision: start with the available finite FIN model and controlled reduction.
Conditional branch: RG-M00, evidence recovery and model freezing for FIN-M0.
Language: English for the execution agent; existing source reports retain their original language.

## 1. Diagnosis and the effect of the new instructions

The 2026-10-14 additions to AGENTS.md are binding. In particular:

- FIN-M0 execution, Kaggle sweeps, saturation time, SOC and soliton claims have no available source/trajectory package in this checkout. They are reported hypotheses or experiments awaiting evidence.
- The triangle-reward derivative is nonpositive, not strictly negative in every configuration.
- The relational-generative chain RG-P01–RG-P08 is a proposed dependency chain.
- The post-326 reconciliation qualifies physical-bridge PASS language and retains the unpaid barrier proof gates.
- A many-body effect produced by marginalization is accepted as a bounded numerical witness, not a completed universality or field-theory result.
- The preparation lane remains useful, with task337 WIP and task335 unopened.

The current finite-N FIN core is a generalized q=12 matrix Curie–Weiss/discrete vector-spin model with a supplied A7 and a declared microscopic kinetic law. The strongest available bridge to known physics is this model-class identity and its static fingerprints. Physical hardware calibration and an independently sourced A7 remain open.

For the generative research question, the most useful available object is the exact marginal law obtained by eliminating a copy. It creates a concrete test of whether several levels of description can be derived from one microscopic law.

The current scientific objective is:

> Derive an effective interaction and an observed process from the same microscopic FIN, measure the approximation error, and establish which effects distinguish FIN from normalized comparison models.

This objective uses existing sources and needs no assumed physical space. Success establishes controlled generative mathematics. A physical interpretation additionally requires an independently specified preparation, instrument and system mapping.

### Status table

| Object | Current scope |
|---|---|
| Finite-N Gibbs/leave-one-out core | Defined mathematical model |
| Matrix Curie–Weiss identity | Exact within that model |
| Static Fourier slopes and N=2 pair law | Admitted model results |
| Exact copy-elimination formula | Algebraic identity |
| Third alternating log-weight contrast, magnitude about 0.27127058458855 | Reviewed finite numerical witness |
| PHYS-013 TV about 0.0004127833 | Nominal CTMC replay using stored calibration parameters |
| Specific A7 source | NO_NEW_SOURCE |
| Global Gamma=B4 | Unpromoted pending input, LP, minimum and coverage payments |
| Cross-N/cross-g response task337 | WIP |
| FIN-M0 dynamics and trajectories | Missing evidence package |
| Quantum/gravity/physical-particle claims | Downstream bridge gates unfulfilled |

The 70,865-line GEN001-gen956.md is an informal research narrative. Selected definitions can be proposed as inputs, but the numbering and prose do not constitute 956 independently verified results. ROADMAP_TOE_FIN.md contains downstream proposals whose acceptance criteria do not by themselves establish GR, gauge dynamics or particle physics. For example, computing nonzero curvature does not prove an Einstein–Hilbert limit, and inserting a connection to restore a chosen gauge invariance does not derive that connection from the admitted FIN core.

## 2. Three distinct questions

1. **Effective generation:** does eliminating variables produce calculable new operators, memory and predictions?
2. **Relational-model generation:** does a separately frozen model of relations produce reproducible nontrivial configurations?
3. **Physical generation:** do independently justified microscopic laws yield measured structures of the physical world?

The first question is immediately executable. The second requires the FIN-M0 model and evidence contract. The third requires additional source and empirical bridges. Results must name the question they actually address.

A valid minimal theory can start from explicit axioms. Its claim of deriving a structure must still identify which features were already contained in those axioms, the allowed operations, controller, cost or initial ensemble.

The nadsoliton ontology is retained as the project's hypothesis: the nadsoliton itself is primordial information in a solitonic state, with no lower information substrate introduced here.

## 3. The mathematical objects for the first campaign

### 3.1 Freeze one parent population

Let the parent contain M labelled copies, sigma_a in Z12. With zero control field,

P_M(sigma) is proportional to exp[g/(2M) sum_(a,b=1)^M A7[sigma_a,sigma_b]].

The first parent is M=4. Use the existing working gain g=5.145228719489142 and g=0 as a negative control. These are already-used inputs, not a new-g test.

The parent denominator M remains fixed under every marginalization. Eliminating one copy does not authorize replacing it with a new population-normalization convention.

For n counting the M-1 surviving labelled copies, the exact log-weight, up to a constant, is

F_(M-1)(n) = g/(2M) n^T A7 n + log sum_j exp[(g/M)(A7 n)_j].

Constant diagonal entries make the removed self-term irrelevant to normalized probabilities. Count probabilities require the additional multinomial degeneracy of the surviving copies.

### 3.2 Use an unambiguous operator decomposition

On the surviving labelled state space, define uniform-reference conditional-averaging operators E_i over label i. The canonical component supported on a subset S is

F_S = [product_(i in S)(I-E_i)] [product_(i notin S) E_i] F.

This fixes the decomposition gauge:

- F_empty: constant;
- |S|=1: one-body terms;
- |S|=2: pair terms;
- |S|>=3: generated higher-body terms.

Every nonconstant component has zero uniform mean in each of its active coordinates. The decomposition sums exactly to F. It is a representation of the static effective log-weight, not a choice of kinetics.

First compare the exact law with its canonical pair truncation. Add the generated third component only if the pair truncation fails the declared error objective. No parameters are fitted to repair a failed truncation.

### 3.3 A usable probability-error certificate

For exact and approximate unnormalized log-weights F and F_tilde, let

r = F - F_tilde,
omega = max r - min r.

A finite-alphabet bound is

TV(P_F,P_F_tilde) <= tanh(omega/4).

The execution agent must include its proof and compute certified extrema, not rely only on floating sample residuals. A constant residual has omega=0 and cannot change the normalized law.

For the totality of any later fixed observation channel K,

TV(P_F K,P_F_tilde K) <= TV(P_F,P_F_tilde).

Conditional postselection is a separate operation. If used, its normalization and probability of selection must enter the error budget. Prefer retaining an explicit unclassified outcome.

### 3.4 Static reduction and dynamic reduction are different operations

A microscopic model is the pair (P_M,L_M), not only its equilibrium law. The projected process must be derived from L_M.

For a three-copy parent, retain the ordered labels (sigma_1,sigma_2) and hide sigma_3. Two microstates (0,0,0) and (0,0,1) have the same observed state. Their rates for refreshing copy1 can differ through

softmax[(g/3)(A7[:,0]+A7[:,sigma_3])].

Thus a strong Markov closure needs a separate test; it does not follow from the known effective stationary distribution.

Construct and compare:

- the exact observed process from the parent;
- the conditional stationary average of the microscopic generator;
- a separately declared Gibbs sampler constructed from the exact marginal stationary law.

The last two can share an equilibrium law without sharing the observed process. Their roles must not be merged in the report.

For equilibrium orthogonal projection P=conditional expectation onto observed variables, Q=I-P, and the reversible generator L, the coupled memory operator is

K(t) = P L Q exp(t Q L Q) Q L P.

Its nonzero source, initial-state contribution and observation window must be recorded. Report scalar spectral positivity only for properly defined quadratic channels of the reversible representation.

### 3.5 “Copy deletion” is not spatial coarse-graining

This campaign marginalizes exchangeable mean-field copies. It does not introduce a length scale, lattice neighborhoods, physical dimension or a spatial renormalization group. These require separate composition and locality inputs.

## 4. Dependency map and first decision

Available model and audited inputs
→ exact generated static operators
→ probability-error control
→ actual projected dynamics and memory
→ normalized countermodel comparison
→ frozen finite-window prediction
→ checkpoint.

Only after that checkpoint may a new source object, relational candidate or operational composition problem be selected.

The independent proof-support lane for 327 may repair its arithmetic and replay the cover. It does not block the small static/dynamic campaign here. Task337 can continue under its existing protocol; its library reduction is not a substitute for the generator/source questions in this campaign.

## 5. Exact input registry

Paths are relative to the repository root.

| Alias | Input |
|---|---|
| I01 | AGENTS.md, especially its 2026-10-14 sections |
| I02 | fin_post_326_methodology_review/REVIEW_20261014.md |
| I03 | fin_post_326_methodology_review/verification.json |
| I04 | fin_post_326_methodology_review/replay_checks.py |
| I05 | fin_physical_bridge_campaign_01/PHYS-001/DICTIONARY.md |
| I06 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py |
| I07 | fin_physical_bridge_campaign_01/PHYS-003/PROOF.md |
| I08 | fin_physical_bridge_campaign_01/PHYS-004/COUNTERMODELS.json |
| I09 | fin_physical_bridge_campaign_01/PHYS-004/CALIBRATION_CONTRACT.json |
| I10 | fin_physical_bridge_campaign_01/PHYS-005/SOURCE_GRAPH.md |
| I11 | FIN_REALITY_GENERATION_ROADMAP_20260929.md |
| I12 | FIN_PHYSICS_BRIDGE_MASTER_ROADMAP_20260929.md |
| I13 | fin_research_327_current_review/INTAKE_20260928.md |
| I14 | FIN_FULL_HANDOFF_327_CURRENT_20260928/02_RESEARCH_331_336/fin333/FROZEN_CONTROL_PROTOCOL_333.json |
| I15 | GEN001-gen956.md, only the explicitly cited candidate-definition sections near GEN-929 and GEN-944 |
| I16 | fin_physical_bridge_campaign_03/PHYS-014/REVISED_PREREGISTRATION.md |

No missing GEN001–956 archive or absent simulator is required for RG-C02–RG-C05. If an input is absent, report its exact role; do not reconstruct scientific outputs from a prose PASS label.

Freeze the exact mathematical kernel convention and its enclosures before claiming validated bounds. The double implementation I06 provides numerical proposals. It does not by itself certify transcendental feature coefficients. A result for a deliberately rounded exact input matrix must explicitly identify that different model.

## 6. Execution budget and deliverable rules

The architect selects the questions; the execution agent implements and reports them.

- Static parent: at most four copies, 12^4=20,736 labelled states.
- Dynamic parent: at most three copies, 12^3=1,728 labelled states, sparse single-copy updates.
- Start with FIN, trace-matched full Potts and flat rank-7 projector. Add only the already-defined perturbation needed to interpret a particular result.
- Do not run N=13, task335, a broad parameter sweep or a new physical experiment.
- Use the existing gain and g=0 controls. An eventual withheld time/observable is a computational validation, not a new physical-gain holdout.
- Use an independent small verifier for the arithmetic/decomposition certificate. A producer repeating its own PASS flag is insufficient.
- Existing code is input; implement new work under fin_generative_campaign_20261014/.
- Do not alter AGENTS.md or historical source packages as part of execution.

Each task exports REPORT.md, RESULTS.json, NONCONCLUSIONS.md, REPLAY.md, INPUTS.sha256 and certificate/input objects sufficient for its claimed scope.

RESULTS.json must include:
execution_status, scientific_verdict, model_id, input_hashes, assumptions,
derived_objects, supplied_objects, fitted_objects, arithmetic_status,
error_budget, discrepancies, dependencies and next_action.

Use numerical evidence, conditional theorem, validated bound, exact theorem, no-go, open-input and design-only as separate status types. “Completed” never means “physics proved”.

## 7. Task cards

### RG-C01 — Bind the current scientific model and evidence

Class: [A]. Priority: P0. Computational class: S.
Question: which available objects can support this specific bounded campaign?
Purpose: RG-P01, model and evidence contract.
Inputs: I01–I14 and the present plan.
Dependencies: none.

Steps:
1. Hash actual inputs, including the uncommitted methodology review used in this workspace.
2. Reuse the review-time verifier results with their source hashes; rerun only a relevant check whose source has changed.
3. Freeze q, parent M, g, A7 construction, labelled-copy convention, stationary law and actual kinetic law.
4. Freeze the observation map and distinguish labelled and count representations.
5. Specify the enclosure convention and proposal/verifier separation.
6. Export the dependency ledger for every numerical or theorem-level claim.

Deliverables: MODEL_CONTRACT.json, EVIDENCE_LEDGER.md and INPUTS.sha256.
Acceptance: every task uses one explicitly identified parent model, gain and observation map; no missing artifact is silently promoted.
Failure/kill: a dependency is unavailable or model definitions disagree. Mark the affected claim OPEN_INPUT and continue unaffected tasks.
Forbidden shortcuts: trusting report counts as evidence; changing the M denominator under deletion; using an unvalidated high-precision point as an exact-kernel interval.
After PASS: RG-C02 and bounded RG-C03.
After FAIL: repair the contract; a new physical premise requires architectural review.
Must not conclude: A7 is sourced, a simulator exists, or the physics programme is complete.

### RG-C02 — Generated operators and a controlled static truncation

Class: [B]. Priority: P0. Computational class: S/M.
Question: what changes in the effective interaction after 4→3→2 copy deletion?
Purpose: static part of RG-P04.
Inputs: RG-C01, I04–I09 and section3.
Dependencies: RG-C01 PASS.

Steps:
1. Verify the exact 4→3 marginalization formula against direct labelled summation.
2. Fix the alternating-contrast convention: plus at (0,0,0) for label choices {0,1}, {0,2}, {0,4}.
3. Reproduce its nonzero contrast; certify an interval excluding zero if claiming a validated witness.
4. Compute the uniform-reference operator decomposition over all 12^3 states.
5. Construct pair-only and, if required, third-order-inclusive log weights.
6. Report residual oscillation, its TV bound and actual finite-enumeration TV.
7. Compare sequential 4→3→2 with direct 4→2 while retaining the parent denominator.
8. At g=0 require all generated nonconstant corrections to vanish.

Deliverables: GENERATED_OPERATORS.json, certificate data, TRUNCATION_REPORT.md, independent verifier.
Acceptance: exact representation is reconstructed; the numerical/certified residual status is explicit; each truncated law has a valid probability/error statement.
Failure/kill: pair-only closure fails its declared target. Retain generated operators; do not replace them by a refitted gain.
Forbidden shortcuts: labelling generated operators as new fundamental forces; treating count multiplicity as an interaction; claiming spatial RG.
After PASS: RG-C04 and the static part of RG-C05.
After FAIL: report the minimal required operator family and stop unsupported compact-closure claims.
Must not conclude: new space, rank-seven necessity or a renormalization fixed point.

### RG-C03 — Derive the projected process from the same parent

Class: [B]. Priority: P0. Computational class: M, capped at 1,728 states.
Question: does the static effective law suffice to predict the actual observed dynamics?
Purpose: dynamic/memory part of RG-P04.
Inputs: RG-C01, I05–I06, section3.4.
Dependencies: RG-C01 PASS; use the same conventions as RG-C02.

Steps:
1. Build the three-labelled-copy leave-one-out generator from the existing formulas.
2. Verify its row sums, positivity, stationarity and detailed balance.
3. Define the observation as the first two ordered labels.
4. Exhibit or refute the two-microstate strong-lumpability defect analytically.
5. Construct conditional-equilibrium projection using the parent Gibbs law.
6. Compare the averaged coarse generator with a Gibbs sampler constructed from the exact two-label marginal.
7. Compute a small declared set of projected correlations and transition rows.
8. Quantify memory and initial/preparation contributions for declared observable channels.
9. Test finite-window predictions; keep hidden variables or memory if the autonomous approximation fails.
10. Include equilibrium-prepared and one explicitly nonstationary hidden preparation; separate their implications.

Deliverables: PROCESS_MAP.md, LUMPABILITY_WITNESS.json, MEMORY_REPORT.json and finite-window diagnostics.
Acceptance: equilibrium identity and process equality are independently assessed; any Markov approximation has a declared preparation, window and error.
Failure/kill: only equal stationary probabilities or projected first derivatives have been shown. That is insufficient for dynamic closure.
Forbidden shortcuts: selecting a new Gibbs sampler and calling it the actual marginal process; discarding initial slip; symmetrizing a failed detailed-balance check.
After PASS or a scoped memory no-go: RG-C04.
After implementation FAIL: repair the generator before computing fingerprints.
Must not conclude: fundamental time, irreversible information loss or physical causal cones.

### RG-C04 — Countermodels and the value of the generated structure

Class: [A]. Priority: P0. Computational class: S/M.
Question: which effects are specific to the supplied FIN interaction, and which occur throughout the comparison class?
Purpose: RG-P07 restricted to the available finite effective models.
Inputs: RG-C02/03, I08–I10.
Dependencies: C02 and C03 completed with scoped evidence.

Steps:
1. Start with FIN, FULL_POTTS_TRACE and FLAT_P7_TRACE.
2. Use the existing trace/gain calibration rule for all of them.
3. Apply exactly the same marginalization, operator decomposition and observation map.
4. Compare generated higher-order structure, truncation errors, static observables and projected memory.
5. Compare alternative reversible kinetics only as a separately marked sensitivity axis.
6. Classify common mechanisms, kernel-dependent numbers and unresolved effects.
7. Identify at most one contrast worth a withheld computational test.

Deliverables: COUNTERMODEL_COMPARISON.csv and SPECIFICITY_VERDICT.md.
Acceptance: no advantage arises from changing normalization, observables or errors for only one model; positive and negative results are retained.
Failure/kill: a supposedly FIN-specific phenomenon also appears in the controls. Report class-level behaviour rather than FIN necessity.
Forbidden shortcuts: inferring universality from a few examples; accepting rank/symmetry as a source theorem; optimizing for dimension3.
After terminal result: RG-C05, then the checkpoint.
Must not conclude: physical selection of A7 or a unique generative ontology.

### RG-C05 — Freeze one finite-window computational prediction

Class: [A]. Priority: P0. Computational class: S.
Question: can the chosen reduced description predict an observable/time-window not used to select its truncation?
Purpose: bounded RG-P08.
Inputs: RG-C02–C04, I07 and the full observation/preparation contract.
Dependencies: a candidate with a stated error budget.

Steps:
1. Before generating the test output, freeze the retained operators, model, preparation, observation and candidate window.
2. Reserve one previously unused observable or finite time window at the same declared gain.
3. Sum input, truncation, numerical propagation, memory/preparation and readout errors.
4. Compare the frozen prediction with a separately regenerated parent-model result.
5. Use controls from C04 under the same protocol.
6. Record whether this was an exhaustive consistency check, a computational holdout or empirical data. The first two cannot become laboratory validation.

Deliverables: PREREGISTRATION.json, PREDICTIONS.json, VALIDATION.json and complete errors.
Acceptance: the actual test discrepancy respects the predeclared bound, or a scoped failure is explained without changing the original task identifier.
Failure/kill: an extra parameter or a revised preparation is needed after opening the test.
Forbidden shortcuts: retuning the error bound after the test; using task335's unopened gain; calling synthetic evidence physical-source validation.
After result: stop and return RG-C06 to the architect.
Must not conclude: a physical world has been generated.

### RG-C06 — Architectural checkpoint

Class: [S]. Priority: P0. Computational class: S.
Question: what is the strongest justified next physics question?
Inputs: C01–C05 handoff, claim register and failures; current AGENTS.md.
Dependencies: preceding tasks have terminal, scoped dispositions.

The architect answers:
1. Are static and dynamic reductions mutually consistent?
2. Is a small effective description accurate enough for its declared observables?
3. Does any contrast depend materially on A7 rather than generic mean-field behaviour?
4. Is a new source object needed, or is the next useful action an engineered measurement?
5. Does the relational model have a reproducible contract capable of a separate degeneracy test?

Deliverable: CHECKPOINT_VERDICT.md with one selected branch:
- CONTROLLED_FINITE_EFFECTIVE_THEORY;
- RETAIN_EXPLICIT_MEMORY;
- FIN_M0_MODEL_GATE;
- ENGINEERED_MEASUREMENT_GATE;
- NO_NEW_FUNDAMENTAL_SOURCE.

Acceptance: one bounded next step, explicit prerequisites and a stop condition.
Kill: proposing hundreds of downstream tasks to avoid a negative source or closure result.
Forbidden shortcuts: automatic particle/QM/GR promotion; adding a graph or clock at each level and calling it derived.
After verdict: the next campaign is separately scoped. The execution agent cannot make this architectural choice.

## 8. Conditional branch RG-M00: restore or freeze FIN-M0

This branch implements RG-P01–RG-P03 only after an executable model is available. It is independent of C02–C05.

**Inputs required:** actual producer, all state variables, cost, transition rule, controller/update schedule, initialization, boundary/projection rule, memory equation, numerical settings and recorded trajectories or a replay recipe.

If the historical simulator cannot be recovered, a new candidate may be defined from explicitly selected hypotheses in GEN001-gen956.md. Give it a new model ID; do not count it as a replay or as evidence that the original experiment happened. The selection of a new cost/dynamics remains a modelling premise requiring architectural review.

### Triangle and degree negative controls

For the stated undirected nonnegative W:

d[-alpha Tr(W^3)/6]/d w_ij = -alpha (W^2)_ij.

This derivative may vanish. With bounded edges, the triangle term favours dense triangular connections. With unbounded edges, boundedness of the cost is itself a required check.

Adding a degree constraint does not automatically solve the locality problem. For a simple d-regular unweighted graph on n vertices, the triangle count satisfies

T <= n d(d-1)/6.

Each triangle is counted at three vertices, and every vertex contributes at most choose(d,2). Equality requires each neighbour set to be complete, so every connected component is a clique K_(d+1). Equality is achievable by a union of such cliques when n is divisible by d+1.

This is a conditional combinatorial control, not a statement about the missing weighted/noisy/memory simulator. It shows why “sparsity plus triangles” can favour fragmentation rather than extended geometry. Connectivity, degree or transport constraints must remain explicitly listed as supplied premises when used.

### Initial execution scope

1. Test empty, complete, uniform, disconnected-clique and random configurations.
2. Use the same permutation-invariant definition of observables across seeds.
3. Verify analytic gradients and stationary candidates before trajectories.
4. Compare at least three time steps or demonstrate exact discrete updating.
5. Include memory-off and specified memory-on dynamics as distinct contracts.
6. Use small n=4,6,8 before any larger run.
7. Record whether an apparent result is sensitive to clipping or update order.
8. Define stochastic entropy on the actual ensemble; the Frobenius pattern of W alone has no automatic Shannon entropy.
9. Before adding a competing cost term, state the independent motivation and freeze one new model.
10. If only trivial states survive, close that candidate at Checkpoint A.

Deliverables: MODEL_CONTRACT.json, baseline proofs, finite trajectories, replay and NEGATIVE_RESULTS.md.
Success: a reproducible finite model with honest outcomes. It is not a particle or space certificate.
Missing package: OPEN_INPUT. Work on the available finite-N branch continues.

## 9. Later bridges and activation conditions

| Next layer | Minimum activation condition |
|---|---|
| Persistent pattern | Nontrivial reproducible states after model/degeneracy checks; identity and lifetime modulo permutation defined |
| Pattern composition | Joint protocol derived from the same admitted dynamics and compared with independent composition |
| Operational locality | Intervention influence and propagation derived with a size-controlled support/rate law |
| Effective geometry | A defined distance or directed cost with its actual axioms; stable under the declared coarse map |
| Scaling/continuum | A family and order of limits; multiple independent diagnostics and countermodels |
| Effective fields | Explicitly retained slow/conserved variables; gradient/source terms derived from the composition law |
| Stable effective matter | Persistent propagating patterns with stability and interactions, beyond a stationary minimum |
| Physical prediction | Independent preparation/readout/units mapping and data outside calibration |
| Quantum/gravity | Separate bridge premises and tests that earlier layers do not supply |

Nothing here selects a physical dimension, fundamental velocity or particle number by target fitting.

The immediate missing laws remain: a physical source for A7, a kinetic/controller premise, simultaneous composition/incidence, protection of slow conserved quantities, and an operational mapping to independently collected data. Their dependencies should be resolved in this order only when the relevant lane has its required inputs.

## 10. Existing parallel lanes

### Barrier 327

The required repair is narrower than “redo everything at higher precision”:

- bind exact mathematical kernel/features and their enclosures;
- certify LP extrema, for example with a dual bound rather than a rounded allocation;
- use a certified analytic scalar minimum rather than point evaluation at a rounded argmin;
- interval-evaluate the Krawczyk centre residual and derivatives;
- bind local cube, path and global cover to the same inputs;
- retain complete, independently checkable coverage records.

Gamma=B4 and an unconditional capacity exponent remain conditional until this chain passes. No extra metastable prefactor claim may borrow an unpaid global premise.

### Preparation 337/335

The current preparation controller is finite-time, while some response weights use biased equilibrium. Carry controller mismatch explicitly or use finite-time weights.

Use the symmetry stabilizer of the chosen seed for a single preparation. Full D12 covariance relates different seeds; it does not automatically identify them within the same controlled experience.

The 78-state low-defect set has a finite-N domain. Freeze a complete response/dynamics predictor with tail and propagation errors before task335. Existing N-specific tables do not supply a cross-g law.

### Hardware PHYS-014

Real q12 precalibration is an independent route. Nominal CTMC replay, programme accuracy and total record count do not prove a hardware TV gate or an independent physical source for A7. Keep calibration and validation separate, including dependence between trajectory samples.

The present continuation authorizes planning and local finite checks only; no external experiment, spending or messaging is inferred.

## 11. Handoff instruction for the execution agent

~~~text
Read AGENTS.md, especially the 2026-10-14 additions, and
FIN_GENERATIVE_NEXT_CAMPAIGN_20261014.md.

Execute RG-C01–RG-C05 only in fin_generative_campaign_20261014/.
Use at most four copies for static elimination and three labelled
copies for sparse dynamics. Use the already-declared gain and g=0.
Do not open task335 or introduce a physical graph by interpretation.

Start from the versioned safe-seed implementation and the methodology
review. Do not treat GEN001-gen956.md as an accepted claim register.
Reuse completed reviews where inputs have not changed.

Track the exact parent normalization under every elimination.
Give both generated static operators and the actual projected dynamics.
Retain memory or hidden preparation variables when closure fails.

Use trace-matched full Potts and flat-P7 controls.
Before a withheld computational observation, freeze the whole
prediction and all error components.

For every task deliver REPORT.md, RESULTS.json, NONCONCLUSIONS.md,
REPLAY.md and input/certificate hashes. Save failed controls.
Distinguish execution_status from scientific_verdict.

After RG-C05 export HANDOFF.md, CLAIM_REGISTER.json and MANIFEST.sha256,
then stop. RG-C06 is the architect's decision.

Do not run FIN-M0 sweeps unless a separately frozen executable model
and its evidence contract are available. If absent, mark OPEN_INPUT.
Do not modify AGENTS.md, acquire external hardware or claim physical
source/QM/GR/ToE closure.
~~~

## 12. The recommended first three tasks

1. RG-C01: bind the available model and evidence.
2. RG-C02: identify generated operators and quantify static truncation.
3. RG-C03: determine what extra memory/dynamic information is required by the same reduction.

The first checkpoint follows the controls and frozen observation in RG-C04/05. A missing FIN-M0 package does not prevent this work.

The central test is whether one microscopic FIN supplies accurate effective laws with a bounded operator and memory ledger. A positive result supports a controlled theory of a model; independent source and empirical bridges determine how far its physical interpretation can go.
