# FIN after RI: microscopic parity likelihood and bounded predictive sufficiency

Date: 2026-10-15. Baseline: `21a875c5` plus preserved UF/MT/OBS/RI intakes.
Status: **RESEARCH_PLAN; LI-001–LI-005 have not been executed.**
Execution directory: `fin_parity_record_sufficiency_campaign_20261015/`.
Language: English for the execution agent.
Decision: [LI_DECISION.md](fin_record_information_review_20261015/LI_DECISION.md).

## 1. One question, not another precision campaign

Does the **entire microscopic one-tag parity history** permit a nonvanishing
extra improvement in predicting the later collective sign, beyond the
already derived continuous classical observation experiment?

Target only the declared lambda0, k1, observation horizon 1/2 and future
target time 1. Desired theorem:

    r_N(full parity history to 1/2) -> r_*(continuous Y to 1/2),
    target C_N=sign(Z_N(1)), retaining ZERO at finite N.

This is a bounded 0-1 decision-risk bridge. It is **not** universal Le Cam
equivalence, convergence of every microscopic filter, entropy preservation,
a physical observer or a fundamental measurement postulate.

RI already proved the fixed nine-record theorem and qualitative advantage.
Do not recompute its risk to more digits. The limiting filter is available;
do not re-derive it as the new result. A new graph, controller, medium,
kinetic law, noise source or external experiment is neither needed nor allowed.

## 2. Exact input registry

| Alias | Repository input |
|---|---|
| L01 | AGENTS.md, current RI/OBS/MT/UF and unchanged physical guardrails |
| L02 | fin_record_information_review_20261015/INTAKE_20261015.md, CLAIM_REGISTER.json and verification.json |
| L03 | fin_record_information_review_20261015/PROOF_COMPLETIONS.md |
| L04 | fin_record_information_review_20261015/TAG_FEEDBACK_REMOVAL_LEMMA.md |
| L05 | fin_record_information_campaign_20261015/RI-001_RECORD_FIELDS_AND_TAG_GENERATOR.md and RI-004_FULL_RECORD_GATE.md |
| L06 | fin_record_information_campaign_20261015/RI-002_FINITE_RESOLUTION_THEOREM.md, RI-003_CAUSAL_FILTER.md and FIXED_EXPERIMENT.json |
| L07 | fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md and REVERSIBILITY_AND_RECORD_SCOPE.md |
| L08 | fin_critical_metastability_review_20261015/INTAKE_20261015.md and PROOF_COMPLETIONS.md |
| L09 | fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md and fin_uniform_critical_review_20261014/INTAKE_20261014.md |
| L10 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; fin_record_information_review_20261015/verify_intake.py |

Hash actual inputs. If an uploaded wrapper lacks a required review proof,
mark that missing dependency; do not invent its hash or historical access.
Obtain the actual files from the supplied workspace/handoff or retain an
explicit conditional statement. Do not substitute an old plan for an intake.

## 3. Frozen primary process and new analysis references

The primary FIN remains the q12 labelled unit-clock LOO heat bath with exact
`g_N=g_c`, iid-uniform initialization and deterministic post-MT observation
origin t_N. The complete horizon through target time must fit the MT window.
Use the inherited local-Gibbs lift only by the paid TV path comparison;
never discard outside/escaped trajectories from the actual ensemble.

Distinguish three laws, never silently calling them the same microscopic FIN:

1. P_N: unchanged interacting labelled FIN.
2. R_N: an **auxiliary** law omitting the one tagged value from background
   refresh fields, with original denominator N and clocks unchanged. The
   tag still responds to background. L04 already bounds full-path TV to P_N
   by O(N^(-1/4)) over sqrt(N) times each fixed horizon, for equal starts.
3. Q_N: autonomous background as in R_N, but tag parity is a fair two-state
   flip chain of rate 1/2. After the separate initialization payment, its
   initial sign is fair and independent of background. This is a reference
   probability measure, not new microscopic physics.

All tagged q12 degrees may be projected away **after** feedback deletion:
under R_N the background/tag-parity pair is Markov. Full q12 label paths and
refresh-attempt marks are not observed. Do not use their richer information
to solve a parity-only question.

Write n_bg for N-1 background counts, but keep

    p_bg=n_bg/N, alpha_N=(N-1)/N,
    m_bg=sigma^T p_bg,
    w_bg=p_bg-alpha_N u-m_bg sigma/12,
    Z_bg,N=N^(1/4)m_bg,
    q_N=sigma^T softmax[(g_c/N) A n_bg].

Missing population mass is in alpha_N u; A u=0. Do not silently normalize
the background to a new N-1 parent and borrow a different gain/clock.
The initial background marginal is the **actual** local-Gibbs marginal after
tag elimination, not a newly supplied Gibbs measure without that tag.

## 4. A concrete likelihood route

Under Q_N let S_t be the fair tagged sign, J_t its flip count, and

    M_S(t)=S_t-S_0+integral_0^t S_r dr,
    <M_S>_t=2t,  M_J=J-t/2.

These clocks follow from the existing unit-attempt convention. The reference
integrated record is `Y_N(s)=N^(-1/4) integral_0^(sqrt(N)s) S_t dt`.

Conditionally on autonomous background, the biased-parity likelihood to
model time T_obs=sqrt(N)/2 is exactly

    log L_N=integral log(1-S_(t-) q_N(t)) dJ_t
             +(1/2) integral S_t q_N(t) dt.

All integrands are predictable with left limits where needed. In particular
the leading score is the **exact** identity

    -integral S_(t-) q_N(t) (dJ_t-dt/2)
       =(1/2) integral q_N(t) dM_S(t).

This avoids imposing an independent measurement-noise law or treating a
finite-variation Y_N as Brownian before proving its martingale replacement.
The proposed joint reference limit is

    (Z_bg,N,Y_N,log L_N)
       => (Z,sqrt(2)W,
           (1/2) integral_0^(1/2) Z_s dY_s
             -(1/4) integral_0^(1/2) Z_s^2 ds),

with Z the admitted stationary quartic diffusion and independent reference
W. These are **targets** except the exact finite likelihood/score identity.
Likelihood tilting would then produce the already admitted observation drift.

A joint likelihood limit is not the final decision theorem. Conditional
densities given full parity histories must be approximated in L1 by functions
of finite observed cylinders. This is the separate LI-004 acceptance gate.

## 5. Task cards

### LI-001 — Pay the autonomous-background and initial-law reduction

**Title:** One tag may be negligible in feedback, not automatically in information.
**Scientific question:** can the actual prediction experiment be compared
with an autonomous-background experiment without changing FIN's premises?
**Why it matters:** a reference measure must be justified before importing
independent-signal filtering tools.
**Exact inputs:** L01–L10, especially L04 and MT's local moment/boundary law.
**Dependencies:** admitted RI/OBS/MT; no new kinetic/source object.
**Method:**

1. Bind P_N, t_N, tag, target, exact gain, denominator, clock and horizon.
   Reuse the proved path-TV feedback bound through **time 1**, not just 1/2.
2. Start the comparison at the actual lifted pi_D,N using the inherited
   coupling, with its initial TV error and future exit error retained.
3. Derive the conditional tagged-label law given n_bg from the constant-
   diagonal Gibbs weight. Where all twelve completions n_bg+e_j lie in D_N,
   it is exactly softmax[(g/N)A n_bg]. Where completions cross the domain,
   retain the truncation instead of dropping it.
4. Pay the bad-completion mass using MT's true departure boundary layer:
   an allowed completion and a forbidden completion differ by one legal
   tagged count move. The corresponding actual state is a boundary state.
5. Use MT moments to bound the good-region conditional parity bias by
   C E(|m_bg|+||w_bg||+1/N)=O(N^(-1/4)). Prove the product-start TV payment
   against the unchanged background marginal times a fair initial sign.
   Background autonomy by itself would not give initial independence.
6. Retain the actual eliminated-tag background marginal; do not replace it
   by a re-normalized N-1 Gibbs model. Bind soft/fast moments and the reference
   critical path law. Reuse path-TV transfer for weak convergence; pay any
   needed unbounded moment estimates natively via the background generator.
7. Replace the full count target by its background sign only after paying
   the event `|Z_bg,N(1)|<=N^(-3/4)`. Weak convergence to a density gives its
   probability ->0; it is not a finite-N rate or a discarded ZERO outcome.

**Deliverables:** REFERENCE_PROCESS_CONTRACT.json, INITIAL_PRODUCT_PAYMENT.md,
TARGET_TRANSFER.md, BATH_MOMENT_AND_PATH_SCOPE.md, INPUTS.sha256.
**Acceptance criterion:** a TV comparison of the bounded prediction problem,
including correlated-start/domain errors, or the exact missing initial lemma.
**Failure/refutation criterion:** feedback is removed without a law-level
bound, initial correlations are assumed absent, or the denominator/target changes.
**What becomes possible:** an actual fair-parity reference likelihood for
the same bounded target, not an external-signal fiction.
**Must NOT conclude:** a physically decoupled subsystem, locality, external
detector, new kinetic law or microscopic information independence.
**Estimated computational class:** S/M, analytical comparison.
**Priority:** P0.

### LI-002 — Exact parity likelihood and paid score remainders

**Title:** Replace arbitrary microscopic timestamps by a source-derived score.
**Scientific question:** which terms in the complete parity-history likelihood
survive critical scaling?
**Why it matters:** this tests informational content, not only amplitude paths.
**Exact inputs:** LI-001, L03–L05, section 4's exact fair-flip identities.
**Dependencies:** autonomous background and paid initial-parity convention.
**Method:**

1. Derive the finite conditional point-process likelihood, including the
   initial parity mass/product-law payment. Retain full killed diagonals in
   the alternative hidden-state representation; do not normalize twice.
2. Expand log(1-Sq) to second order on a declared compact-soft/fast stopping
   region. The first-order score is `(1/2) integral q dM_S` exactly.
3. Split the quadratic term into `-(1/4) integral q^2 dt` and
   `-(1/2) integral q^2 dM_J`. Pay the latter bracket and the cubic remainder
   with the actual rate1/2 flip count, not a rate1 attempt count.
4. Derive a sufficient bound on q_N-m_bg from the existing softmax/fast
   structure. A bound `C(||w_bg||+|m_bg|^3+1/N)` suffices at lambda0.
5. Pay `E integral_0^T (q_N-m_bg)^2 dt ->0` on the declared scope. Expected
   fast occupation over model time sqrt(N) is relevant. Do **not** replace
   this time-integrated bound by an unproved uniform scaled-fast estimate.
6. State all stopping, initial tail and background-jump contributions.
   Any clipping/killing is a reference proof device with a paid transfer,
   not an undisclosed controller of the actual experiment.

**Deliverables:** EXACT_PARITY_LIKELIHOOD.md, SCORE_EXPANSION.md,
REMAINDER_LEDGER.json, PROOF_OR_OPEN_LEMMA.md.
**Acceptance criterion:** exact likelihood plus explicit vanishing remainder
estimates, with every reference/original law distinguished.
**Failure/refutation criterion:** a hidden mark carries the claimed information,
an unobserved attempt is treated as a parity flip, or a fast remainder survives.
**What becomes possible:** a joint stochastic-integral likelihood limit.
**Must NOT conclude:** full record sufficiency, filter convergence, entropy
conservation, quantum amplitudes or a physical noise bath.
**Estimated computational class:** M, generator/martingale estimates.
**Priority:** P0.

### LI-003 — Joint reference score limit, mass and uniform integrability

**Title:** Turn the microscopic likelihood into the correct limiting experiment.
**Scientific question:** does its entire likelihood converge without loss of mass?
**Why it matters:** a formal log expansion is not a probability-measure bridge.
**Exact inputs:** LI-001–002, inherited compact-stopped martingale machinery.
**Dependencies:** paid references, time-integrated score remainders and initial law.
**Method:**

1. Under the fair product reference, prove the joint background critical law
   and additive tagged martingale limit. The exact relation between Y_N and
   N^(-1/4)M_S has a uniformly bounded endpoint error. Do not identify their
   path sigma-fields with those of the limit without this next likelihood work.
2. Identify the score jointly through its bracket, cross bracket with the
   record martingale, and zero cross bracket with background jumps. Work in
   the correct joint natural filtration and pay square-martingale UI on stops.
3. Prove joint convergence with the **future** Z_bg,N(1) as well, although
   likelihood weighting/observation stops at 1/2. Future signal is a target,
   not extra observed data.
4. Obtain the limiting log likelihood in section 4 and verify its coefficients
   1/2 and 1/4 against observation variance2. Its conditional mean likelihood
   given the background path must be exactly1 under the reference law.
5. Exploit a simpler possible UI route: L_N>=0, E_QN L_N=1, convergence in
   distribution to L>=0 and E_Q L=1 imply UI. Prove the premises and the
   truncation argument; do not use an unexplained UI flag or infer exponential
   moments from o(1) TV. This result follows by testing min(L_N,K) then K->infinity.
6. Remove localization carefully. For an event depending only on background,
   conditional mean-one likelihood gives `E_QN[L_N 1_bad]=Q_bg,N(bad)`.
   This can pay weighted rare-background tails without a new global moment
   hypothesis. Record any initial-state/error terms not covered by that identity.

**Deliverables:** JOINT_REFERENCE_LIMIT.md, MASS_AND_UI_PROOF.md,
LOCALIZATION_AND_TARGET_LEDGER.json, NONCONCLUSIONS.md.
**Acceptance criterion:** a joint likelihood/target/record theorem with paid
mass/UI, not merely the already known OBS joint state/record limit.
**Failure/refutation criterion:** mass escapes on rare histories, score terms
retain extra noise, or the limit is not measurable in the asserted joint paths.
**What becomes possible:** conditional likelihood comparison under a product
reference, which LI-004 must still prove before any risk-equivalence claim.
**Must NOT conclude:** statistical sufficiency from UI alone, arbitrary-start
universality, physical action, fundamental noise or empirical validation.
**Estimated computational class:** M, analytical proof; no trajectory sweep.
**Priority:** P0.

### LI-004 — Pay conditional finite-cylinder projection and the Bayes-risk bridge

**Title:** Prove absence of extra microscopic predictive information, or expose it.
**Scientific question:** does the full microscopic parity sigma-field improve
the stated bounded forecast by a nonvanishing amount over continuous Y?
**Why it matters:** this is the actual RI-LIKELIHOOD-UI decision gate, not a
consequence of an excellent likelihood plot.
**Exact inputs:** LI-001–003, RI finite-record theorem and limiting filter.
**Dependencies:** reference product structure, exact likelihood, joint convergence
and UI; terminal source/target scopes must match.
**Method:**

1. Write full-record target densities under the fair reference as
   `gamma_c,N(O)=E_bg,N[1{C_bg,N=c} L_N(bg,O)]`, with no hidden independence
   assumption beyond the paid product reference. Keep all three target classes.
2. Express the optimal risk as `1-E_QN max_c gamma_c,N(O)`. Optimality is over
   all measurable parity-history decision rules, not an enlarged q12/mark observer.
3. Approximate the limiting likelihood in L1 by bounded continuous functions
   of finitely many background and observed path coordinates. Dense finite
   grids generate the continuous-path sigma-fields. Explain that this is a
   proof approximation, not a changed microscopic measurement law.
4. Pull the approximation back to finite N using joint convergence and UI.
   Then use L1 contraction of conditional expectation. Product references
   permit integrating out background; pay convergence/equicontinuity of the
   resulting finite-cylinder target densities, including the null sign boundary.
5. This is a possible route, not a free inference: provide the approximation
   and each limit explicitly. If it fails, identify the microscopic statistic
   that survives or the exact missing conditional-projection estimate.
6. First take N->infinity at fixed cylinders, then refine cylinders. Prove
   the full continuous-Y risk is the limit of its nested finite-grid risks
   by conditional-martingale L1 convergence. No uniform-in-history supremum
   or unlicensed exchange of resolution and N limits is required.
7. Transfer the bounded risk back through LI-001 to the actual iid-start FIN.
   State exactly the admitted target, time, initialization and kinetic scope.

**Deliverables:** CONDITIONAL_CYLINDER_COMPARISON.md,
BOUNDED_RISK_SUFFICIENCY_THEOREM_OR_BLOCKER.md, LIMIT_ORDER_LEDGER.json,
NEGATIVE_RESULTS.md.
**Acceptance criterion:** the stated equality of full-parity and continuous-Y
limiting forecast risks, a genuine FIN-specific refutation, or a precise open
conditional estimate. The last disposition does not promote sufficiency.
**Failure/refutation criterion:** a bounded forecast advantage remains in fine
microscopic histories, or only joint weak convergence/UI has been shown without
the conditional comparison. A generic C/N toy is not a FIN refutation.
**What becomes possible:** a controlled predictive observation theory for this
target, or a concrete minimal statistic needed beyond the proposed record.
**Must NOT conclude:** all-filter or entropy/MI equivalence, full-label-history
sufficiency, microscopic reversibility-to-unitarity or physical measurement closure.
**Estimated computational class:** M, hard proof gate; no brute-force escalation.
**Priority:** P1, after the three P0 tasks.

### LI-005 — Scope handoff and architect checkpoint

**Title:** Return one paid bridge or one genuinely remaining obstruction.
**Scientific question:** what observational theory follows from the same FIN law?
**Why it matters:** bounded prediction, full process equivalence and physical
observation are distinct scientific claims.
**Exact inputs:** all LI results/failures and L01–L10.
**Dependencies:** terminal scoped dispositions, not forced proof success.
**Method:**

1. Export the dependency DAG from P_N to initial/background comparison,
   fair reference, score, likelihood mass, conditional cylinders and bounded risk.
2. Separate reused RI/OBS/MT/UF facts, the intake feedback lemma, newly proved
   limits, diagnostics and any unpaid comparison step.
3. Preserve finite ZERO outcomes and actual initial/exit exceptions. State
   that no practical N0, finite-N risk-error rate or physical clock was obtained
   unless independently proved in this bounded campaign.
4. Reuse stationary reversal/nonconservation boundaries. Information in a
   classical observation is not an irreversible physical record or a sourced
   observer. Canonical quartic behavior does not uniquely select A7.
5. Return one exact next atom and stop for architecture review. No automatic
   further precision, entropy, quantum, gravity or external-data programme.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
NEXT_OPEN_ATOM.md, REPLAY.md, INPUTS.sha256, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** a reproducible scoped result or a precise blocker,
with theorem versus proof-submission status distinguished.
**Failure/refutation criterion:** a reference is promoted to a new microscopic
law, a conditional step is silently skipped, or physical interpretation outruns proof.
**What becomes possible:** one evidence-based next architectural decision.
**Must NOT conclude:** physical source/units/locality/conservation/EV,
selector/role transfer, FIN-M0/B4/task335, QM/GR/SM/L_total or ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution and handoff rules

Begin LI-001–003; LI-004 is a separate explicit proof gate. Do not replace an
open estimate by more Monte Carlo, an optimized PDE grid or a new noise model.
Tiny exact combinatorial and actual-rate checks may verify conventions, not
substitute for growing-horizon likelihood arguments.

Work only in the new directory; preserve AGENTS and all historical packages.
Export per-task REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,
NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. RESULTS must distinguish
execution/science status, original/reference laws, supplied/derived/fitted
objects, initial/target/record/time scopes, arithmetic assumptions, errors,
failed approaches and next action. Return a proposed AGENTS patch.

No large-state enumeration, broad N/g sweep, extra force/drive/graph/storage
medium, experiment, purchase or external message. Read any theorem actually
used and pay its hypotheses. Classical change-of-measure references are
analysis tools, not physical interventions or new ontological layers.

Success would close a substantial **predictive statistical bridge** for the
supplied finite-family model. It would not derive its kernel, kinetics,
gain/initial tuning, units, physical composition/locality, conserved energy,
apparatus or natural reality. These independent physical gaps remain open;
no naming change to a likelihood, innovation or conditional law closes them.
