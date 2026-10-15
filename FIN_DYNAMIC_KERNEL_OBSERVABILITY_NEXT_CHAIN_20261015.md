# FIN after KO: dynamic sign spectroscopy and the cost of partial observation

Date: 2026-10-15. Baseline:380e9dc0 plus preserved LI/KO intakes.
Status: **RESEARCH_PLAN; DKO-001–DKO-005 have not been executed.**
Execution directory: `fin_dynamic_kernel_observability_campaign_20261015/`.
Language: English for the execution agent.
Decision: [DKO_DECISION.md](fin_kernel_observability_review_20261015/DKO_DECISION.md).

## 1. A new question, not another static digit campaign

KO sufficiently answered the fixed static comparison: common finite parity
records, different fast half-block laws, and a positive bounded ideal contrast.
No more precision of that event or global saddle work is the default.

The new question is:

> Can temporal **sign-only** observations identify the fast spectral ratios
> independently of amplitude units and a common clock-unit choice, and what
> fraction/time resources does partial readout require?

This activates a separately scoped use of the formerly deferred fast-time
gate. KO-004 remains historically DEFERRED_NOT_REQUIRED. No process theorem
is borrowed from its static covariance. No extra microscopic law, physical
detector, graph, drive, storage medium or external experiment is introduced.

Empirical promotion remains OPEN_INPUT without independent physical source/
preparation/readout evidence. That does not prohibit the theoretical work
below on consequences of the already supplied FIN generator.

## 2. Exact input registry and deduplication

| Alias | Repository input |
|---|---|
| D01 | AGENTS.md, latest KO/LI and unchanged physical guardrails |
| D02 | fin_kernel_observability_review_20261015/INTAKE_20261015.md, PROOF_COMPLETIONS_AND_ARITHMETIC.md, CLAIM_REGISTER.json and verification.json |
| D03 | fin_kernel_observability_campaign_20261015/KO-002/EXACT_FAST_READOUT_KERNEL.md and KO-002/JOINT_SNAPSHOT_LIMIT.md |
| D04 | fin_kernel_observability_campaign_20261015/KO-001/STRUCTURAL_HYPOTHESIS_AUDIT.md; fin_parity_record_sufficiency_review_20261015/KERNEL_FAMILY_SEED.md and KERNEL_FAMILY_SEED.json |
| D05 | fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md and INTAKE_20261015.md |
| D06 | fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md; fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md |
| D07 | fin_intrinsic_observation_campaign_20261015/OBS-C03_RAW_PATH_NO_GO.md and OBS-C04_INTEGRATED_RECORD_THEOREM.md |
| D08 | fin_theory_synthesis_20261014/GEN_EVIDENCE_MAP.md and CORE_DERIVED_BRIDGES.md; fin_intrinsic_theory_review_20261014/INTAKE_20261014.md |
| D09 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/EXCHANGEABLE_GROUP_COMPOSITION_NOGO_148.md |
| D10 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; fin_rank7_followup/inputs/fin_handoff_audit/results.json |

Hash actual inputs. Reuse the already known drift/noise, stable fast
linearization, stationary Gaussian covariance, exchangeable-group bare decay
and parity private-OU arguments. A new proof must address the **joint fast
partial-readout path in the critical local preparation**, not relabel a known
stable-phase OU result as new. If an exact-scope theorem already exists,
bind it and omit its reconstruction.

## 3. Freeze the parent, record and orders of limits

Primary model is unchanged q12 labelled unit-attempt LOO HB at exact g_c,
N in12N, iid-uniform start, deterministic common post-MT origin t_N. Compare
only the already frozen A7 and A7+.01(P3-P4), with separate model IDs. Keep
gain, denominator N and microscopic attempt clocks unchanged.

Let B_N be a predetermined block, k_N/N->alpha. The primary case is alpha=1/2.
For the optional resource benchmark prove the process law for each **fixed**
alpha in(0,1); no uniform alpha_N->0 result is assumed. Use the normalized
Fourier features f_r(j)=sqrt(2)cos(2*pi*r*j/12), r=1,3,4, with

    F_r,N(t)=sqrt(N) sum_j f_r(j)p_j(t_N+t),
    V_r,N(t)=sqrt(N)/k_N sum_(a in B_N) f_r(xi_a(t_N+t)),
    U_r,N(t)=V_r,N(t)-F_r,N(t),
    c_alpha=(1-alpha)/alpha,
    gamma_1=1, gamma_r=1-Lambda_r/Lambda6 for r=3,4.

These use **fast model time t=O(1)**, not critical time sqrt(N)s. The critical
soft amplitude is asymptotically frozen over this window, not set to zero.
The observer is a new q12 block-aggregate readout, not the one-tag parity
observer. Tags are not physical space; internal Fourier harmonics are not
wavevectors or particle species.

The candidate stationary Gaussian process, still unproved here, is

    dF_r=-gamma_r F_r dt+sqrt(2)dB_r,
    dU_r=-U_r dt+sqrt(2c_alpha)dW_r,
    V_r=F_r+U_r,
    C_r(t)=Cov(V_r(t),V_r(0))
          =gamma_r^(-1)exp(-gamma_r|t|)+c_alpha exp(-|t|).

DKO must derive the independent initial laws/drivers and actual brackets.
V_r alone need not be Markov. Only fixed finite observation times/decisions
are used below; weak path convergence is not full microscopic statistical
equivalence of continuous histories.

## 4. The high-value bounded signature

Freeze alpha=1/2 and model-time lag h=1. Define the finite-N sign statistic

    H_r,N=sign(V_r,N(0))*sign(V_r,N(h)), sign(0)=0.

Keep every ZERO finite-N outcome. A nondegenerate Gaussian limit has no
boundary mass, so its expectation is a transferable bounded observable.
For a centered Gaussian pair, derive the angular/orthant identity

    s_r(h)=E H_r=(2/pi)arcsin(rho_r(h)),
    rho_r(h)=C_r(h)/C_r(0).

Now test a *common time-unit reparameterization* t->nu t, not a newly
chosen microscopic rule. Set q=rho_1(h)=exp(-nu h). The proposed invariant is

    rho_r(h)=[q^(gamma_r)+c_alpha gamma_r q]/[1+c_alpha gamma_r],
    rho_r=sin(pi*s_r/2).

For q in(0,1), c_alpha>0, gamma in(0,1), this function is strictly decreasing
in gamma: write it as `q+(q^gamma-q)/(1+c_alpha gamma)` and differentiate.
Thus bounded sign correlations identify gamma uniquely in the limiting
model, using the null mode to remove a common clock gauge.

Time-independent nonzero multiplicative amplitude factors in each channel
cancel from the two-time sign product. Unknown additive offsets, sign errors,
noise, nonlinear instruments or channel-dependent clocks do **not** cancel
by this argument; physical calibration gates remain open. These are ideal
law functionals, not estimates from a single two-time realization.

For the frozen epsilon comparison, lower gamma3 increases rho3 and the
nonzero-sign agreement probability. This gives a prospective bounded dynamic
contrast without transporting unbounded scaled covariance moments through TV.

## 5. Task cards

### DKO-001 — Bind the actual joint fast characteristics

**Title:** One microscopic law, global and private fast variables.
**Scientific question:** which joint drift/brackets follow from the supplied
generator for the new block readout?
**Why it matters:** static Gaussian sampling alone does not determine a process.
**Exact inputs:** D01–D10 and the frozen section3 contract.
**Dependencies:** accepted KO static law and comparator local windows.
**Method:**

1. Freeze models, origin, labels, readout, alpha=1/2, lag h=1 and fast horizon
   at least1. Bind any separately authorized fixed-alpha extension explicitly.
2. Deduplicate existing leading OU/CLT results. Audit their initialization,
   clock, mode and stopping scope against this critical local preparation.
3. Derive L F_r,N from the actual LOO generator. Use exact parity-line
   invariance to eliminate pure soft forcing, then pay O(m) fast-linearization,
   higher-order occupation and finite-N LOO corrections.
4. Generalize `L(m_B-m)=-(m_B-m)+O(1/N)` to the bounded f_r features.
   This supplies private decay, not a general parity-lumpability theorem.
5. Derive every global/global, private/private and global/private bracket
   from actual jumps. For general f the uniform conditional squared jump
   includes `1+f(old)^2`, not parity's constant2. Do not copy that shortcut.
6. Pay mixed brackets by conditional block exchangeability/hypergeometric
   concentration and the local comparison. Initial global/private independence
   is the KO sampling theorem, not an arbitrary Gaussian premise.

**Deliverables:** DYNAMIC_READOUT_CONTRACT.json, DEDUPLICATION_LEDGER.md,
ACTUAL_CHARACTERISTICS.md, BRACKET_AND_INITIAL_PAYMENT.md, INPUTS.sha256.
**Acceptance criterion:** correct actual joint characteristics and paid
initial/error terms, or a precise nonvanishing term invalidating the candidate.
**Failure/refutation criterion:** drift/noise is guessed from one-time covariance,
parity brackets are copied incorrectly, or initialization is silently changed.
**What becomes possible:** a genuinely derived fast observation process.
**Must NOT conclude:** physical field modes, spatial subsystems, new noise,
fundamental clock or changes to primary FIN's microscopic law.
**Estimated computational class:** S/M, symbolic bounds and tiny actual-rate fixtures.
**Priority:** P0.

### DKO-002 — Prove the finite-fast-horizon joint observation law

**Title:** A reusable Gaussian state/readout process, with memory retained.
**Scientific question:** does the proposed global/private system approximate
the actual parent on fixed fast horizons?
**Why it matters:** a positive process-level prediction is stronger than
another matching marginal, while avoiding unproved full-history equivalence.
**Exact inputs:** DKO-001, KO local measure/sampling theorem and MT path coupling.
**Dependencies:** paid drift, brackets, initial law and fixed-alpha scope.
**Method:**

1. Establish tightness and joint martingale identification for all declared
   modes on one fixed fast horizon. Keep the frozen random critical soft
   initial coordinate and pay its vanishing change on that horizon.
2. Use actual jump bounds, moment/localization estimates and square-martingale
   UI in the joint natural filtration. Independence requires joint Gaussian
   limits with paid cross brackets and initial law, not zero covariance alone.
3. Transfer to the actual unconditioned iid-start parent using labelled
   initial/exit coupling. Do not move unbounded moments through an o(1) TV error.
4. Prove or correct the section3 process and covariance formulas under the
   stationary comparison. Identify when the observed V alone is non-Markov;
   retain F,U or the implied memory rather than choosing a new sampler.
5. Export fixed two-time distributions for lag h=1. If proving each fixed
   alpha in(0,1) for DKO-004, state pointwise parameter scope explicitly.
   Do not silently claim alpha_N->0 or increasing-window control.

**Deliverables:** JOINT_FAST_PATH_THEOREM_OR_BLOCKER.md,
OBSERVED_MEMORY_SCOPE.md, TWO_TIME_LAW.json, ERROR_AND_LIMIT_ORDER_LEDGER.json.
**Acceptance criterion:** a paid finite-fast-horizon joint law and bounded
two-time transfer, or a precise failed estimate; not a covariance-only PASS.
**Failure/refutation criterion:** coupling, nonlinearity or bracket leakage
survives, or full microscopic information equivalence is inferred from weak paths.
**What becomes possible:** bounded temporal kernel fingerprints from the same FIN.
**Must NOT conclude:** spatial propagation, physical mass gap, geometry,
energy conservation, irreversible storage or empirical validation.
**Estimated computational class:** M, proof gate; no large matrix enumeration.
**Priority:** P0.

### DKO-003 — Sign-only temporal spectroscopy and unit-gauge tests

**Title:** A bounded kernel signature without amplitude calibration by fitting.
**Scientific question:** can signs of two-time fast readouts identify gamma
after removing a common clock-unit choice?
**Why it matters:** the KO static contrast is amplitude-sensitive; a derived
bounded temporal invariant gives different operational content.
**Exact inputs:** DKO-001–002 and the section4 frozen sign statistic.
**Dependencies:** valid nondegenerate two-time Gaussian law, not just snapshots.
**Method:**

1. Derive the Gaussian sign/orthant identity by angular integration or a
   fully read primary theorem. Prove all finite-N sign-boundary effects vanish.
2. Transfer E H_r,N by bounded weak-law arguments. Keep zero outcomes and
   correlations between different record times; no iid assumption is needed
   to define these ensemble expectations.
3. Use r1 null mode to define q, derive the section4 inversion and prove
   strict monotonicity/uniqueness for the paid gamma range. Derive exactly
   which time/amplitude reparameterizations leave the construction invariant.
4. Freeze the event `V3,N(0)*V3,N(1)>0` as the bounded dynamic comparison.
   Show its strict limiting model contrast without retuning lag, thresholds,
   gain or amplitudes. Numerical values require disclosed arithmetic/error scope.
5. Distinguish identifiable ideal law functionals from a finite-record
   estimator with statistical error. Do not declare a single realization
   an exact gamma measurement or call the internal clock physical seconds.
6. List nuisances that do not cancel: offsets, detector errors, channel-specific
   clocks, disturbance and unknown feature-label mapping. No physical source
   or hardware-readiness inference is permitted.

**Deliverables:** FROZEN_SIGN_CONTRAST.json, BOUNDED_SIGN_LIMIT.md,
CLOCK_AND_AMPLITUDE_GAUGE_THEOREM.md, FINGERPRINT_AND_NUISANCE_SCOPE.md.
**Acceptance criterion:** a strict bounded dynamic contrast and uniquely
defined ideal invariant, or a retained failure of identifiability/gauge removal.
**Failure/refutation criterion:** gamma is confounded even under the declared
unit gauges, the Gaussian pair is singular, or success needs fitted readout noise.
**What becomes possible:** model-specific temporal spectroscopy using derived
classical signals, with a clear physical-interface contract still unpaid.
**Must NOT conclude:** natural constants, physical clock/source, all-instrument
robustness, laboratory discrimination or a new measurement postulate.
**Estimated computational class:** S/M, analytic identities and tiny certificates.
**Priority:** P0.

### DKO-004 — Observation-fraction cost in the limiting Gaussian experiment

**Title:** Can time compensate for a small observed fraction, and at what cost?
**Scientific question:** how does operator information scale with alpha and
the number of dependent fixed-spacing records in the derived limit model?
**Why it matters:** an ideal half-population observation is costly; a claim
that a tiny block solves it must pay an information/time budget.
**Exact inputs:** DKO-002 fixed-alpha process, DKO-003 sign signature, two fixed gammas.
**Dependencies:** process law for each fixed alpha and a declared Gaussian-
experiment scope; no triangular finite-N limit is assumed.
**Method:**

1. Normalize W=sqrt(alpha)V, which preserves signs. At fixed spacing h>0,
   its stationary covariance is

       K_gamma,alpha(i,j)=alpha/gamma exp(-gamma h|i-j|)
                           +(1-alpha)exp(-h|i-j|).

2. For M samples derive the exact Gaussian relative entropy from the two
   covariance matrices. Show uniform eigenvalue bounds for alpha sufficiently
   small at fixed h and positive fixed gammas. Use the private OU component
   for a positive lower bound; include dependence, not an empirical ESS.
3. Since their difference is alpha times an exponentially summable Toeplitz
   sequence, attempt the two-sided bound `c M alpha^2<=KL<=C M alpha^2` for
   distinct gammas. State constants and parameter ranges symbolically if
   sufficient; no broad grid search is needed.
4. The upper KL bound and Pinsker give a no-discrimination result if
   M alpha^2->0, even for full amplitudes and therefore also signs.
5. Do not infer successful discrimination from large KL alone. For an upper
   error bound, derive the nonzero O(alpha) difference of lag-one sign-product
   means and pay its O(1/M) variance by the stationary latent two-OU semigroup
   mixing bound. Alternatively pay a Gaussian Hellinger-affinity bound.
   Establish the actual risk implication when M alpha^2 grows.
6. Take N->infinity first at fixed alpha and finite horizon. Only inside the
   resulting Gaussian family take alpha small or M large. The formal scaling
   M~(N/k)^2 is **not** an actual shrinking-block/growing-time FIN theorem.
   Any such parent-model claim requires a separately scoped uniform limit.

**Deliverables:** GAUSSIAN_RESOURCE_CONTRACT.json, DEPENDENT_RECORD_INFORMATION.md,
LOWER_AND_UPPER_DECISION_BOUNDS.md, LIMIT_ORDER_AND_TRANSFER_GATES.md.
**Acceptance criterion:** a paid resource law for the stated limiting experiment,
or a precise failure; no finite-N/practical sample budget is fabricated.
**Failure/refutation criterion:** record dependence, small-alpha limits or
slow critical modulation are ignored in a claim about the finite parent.
**What becomes possible:** a falsifiable theoretical cost benchmark and a
precise future uniform-limit question, not an experimental shopping list.
**Must NOT conclude:** iid real records, finite-N resource feasibility, physical
time, calibrated sample requirements or kernel selection in nature.
**Estimated computational class:** M, structured Gaussian bounds only.
**Priority:** P1, after the three P0 tasks.

### DKO-005 — Scope handoff and architectural checkpoint

**Title:** A dynamic fingerprint, a cost law, or one exact surviving obstruction.
**Scientific question:** what additional predictive content follows from the
same microscopic FIN without an extra primitive?
**Why it matters:** correct intrinsic spectroscopy must not become an unsupported
claim about apparatus, locality or fundamental law selection.
**Exact inputs:** all DKO results/failures and D01–D10.
**Dependencies:** terminal scoped dispositions, not forced success of DKO-004.
**Method:**

1. Export the DAG from supplied generator to joint characteristics, fast
   process, bounded sign contrast, gauge invariant and separate resource limit.
2. Separate reused static/OU facts, new scope theorems, numerical evidence,
   exact identities and open uniform finite-N estimates.
3. State which effects are class-level and which ratios depend on A7.
   A kernel-sensitive ratio is not a derivation of its value or nature's choice.
4. Preserve actual initial/exit exceptions, all sign-zero outcomes and the
   distinction between fixed horizons and growing observation time.
5. Return one exact next question and stop. Empirical activation still needs
   independent source/preparation/readout evidence and human authority.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
NEXT_OPEN_ATOM.md, REPLAY.md, INPUTS.sha256, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** reproducible scoped results or a precise blocker,
with no automatic extension into another campaign or physical claim.
**Failure/refutation criterion:** covariance is called a field/metric, signs
become quantum outcomes, or a Gaussian-family cost is sold as a finite-N experiment.
**What becomes possible:** an evidence-based decision about a genuinely missing
bridge rather than more ideal static digits.
**Must NOT conclude:** physical A7/kinetic/gain/initial-source, units/locality/
conservation/EV, ontology/selector/legacy-role, FIN-M0/B4/task335 or
QM/GR/SM/L_total/ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution and stopping rules

Start DKO-001–003. DKO-004 is a separate limiting-experiment proof, not a
default actual-parent long-run sweep. Do not launch growing-time/shrinking-
block finite-N simulations to replace the missing uniform theorem.

Work only in the new directory. Preserve AGENTS, sources and reviews;
return a proposed patch. Per task export REPORT.md, RESULTS.json,
PROOF_OR_OPEN_LEMMA.md, NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256.
RESULTS separates execution/scientific status, supplied/reused/derived/fitted
objects, original/comparison laws, readout/time/parameter/initial scopes,
arithmetic assumptions, errors, failures and next action.

No added microscopic clock/noise/graph/drive/controller, source-inventory
loop, B4/FIN-M0/task335 activation, hardware, experiment, spending or messages.
Any literature theorem used must be read and its hypotheses paid. No plot,
large assertion count or precision increase substitutes for the proof gates.

Success would establish a dynamic kernel-sensitive **statistical** observation
theory and a scoped information-cost benchmark. It would not source the
microscopic laws, derive physical spacetime/conservation/units, create a
recording medium or establish QM/GR. Those are independent gaps.
