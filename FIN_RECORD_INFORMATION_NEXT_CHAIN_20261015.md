# FIN after OBS: predictive information in intrinsic classical records

Date: 2026-10-15. Baseline: `21a875c5` plus the preserved UF/MT/OBS intakes.
Status: **RESEARCH_PLAN; RI-001–RI-005 have not been executed.**
Execution directory: `fin_record_information_campaign_20261015/`.
Language: English for the execution agent.
Decision: [RI_DECISION.md](fin_intrinsic_observation_review_20261015/RI_DECISION.md).

## 1. One objective

Determine what the already derived FIN record actually permits one to infer
or predict, and whether the mesoscopic record is statistically sufficient
for that purpose. Derive a bounded observational prediction from the same
microscopic law, without inserting an observer, detector, storage medium,
drive, physical clock or information-conservation principle.

OBS paid a joint weak state/record limit. It did **not** pay convergence of
conditional laws, filters, entropy or full-record statistical experiments.
This is the missing bridge now, not another local-root/gap precision task.

## 2. Fixed parent and inherited results

Use only the supplied q12 A7, labelled unit-clock LOO heat bath, exact
`g_N=g_c+lambda/sqrt(N)`, and the admitted large N in 12N regime. Primary
initialization remains iid uniform, followed by deterministic t_N after the
MT burn-in and before the actual local-sector exit window. Keep every initial
exception and path-coupling/exit error; do not postselect away failed runs.

For each fixed scaled horizon the admitted joint limit is

    dZ=(chi lambda Z-Z^3/3)ds+sqrt(2)dB,
    dY=Z ds+sqrt(2/k)dW,
    Z0~rho_lambda, Y0=0, B,W independent,
    rho_lambda(z) proportional to exp(chi lambda z^2/2-z^4/12).

Here k is fixed. The raw k~sqrt(N) snapshot path is a *different observation*
and is not J1-tight on critical time. A one-time Gaussian channel cannot
replace a history law. All these objects are M-level mathematical statements.

The leading channel has no fitted coefficients. That does not select A7,
derive the microscopic kinetic rule, or turn chi into a physical constant.

## 3. Exact input registry

Paths are relative to the repo root. Hash actual files and access the real
review rather than trusting the name of an archived attachment.

| Alias | Input |
|---|---|
| R01 | AGENTS.md, current OBS/MT/UF and unchanged physical guardrails |
| R02 | fin_intrinsic_observation_review_20261015/INTAKE_20261015.md, CLAIM_REGISTER.json, verification.json |
| R03 | fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md |
| R04 | fin_intrinsic_observation_review_20261015/REVERSIBILITY_AND_RECORD_SCOPE.md and observation_kernels.py |
| R05 | fin_intrinsic_observation_campaign_20261015/OBS-C01_LABELLED_TRANSFER.md and OBS-C04_INTEGRATED_RECORD_THEOREM.md |
| R06 | fin_intrinsic_observation_campaign_20261015/OBS-C02_SNAPSHOT_CHANNEL_THEOREM.md and OBS-C03_RAW_PATH_NO_GO.md |
| R07 | fin_critical_metastability_review_20261015/INTAKE_20261015.md and PROOF_COMPLETIONS.md |
| R08 | fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md and fin_uniform_critical_review_20261014/INTAKE_20261014.md |
| R09 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; FIN_INTRINSIC_OBSERVATION_NEXT_CHAIN_20261015.md |
| R10 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/EXCHANGEABLE_GROUP_COMPOSITION_NOGO_148.md |

The archived file `source_context/MT_INTAKE_USER.md` in OBS is the MT plan,
not the accepted MT review. Do not repeat that provenance error. Earlier
source packages are historical objects; no writer runs against them.

## 4. Observation contracts and two orders of resolution

Distinguish at least:

- the k tagged **parity histories** on [0,S];
- the entire exact integral path Y_N, which can retain its microscopic right
  derivative and is not automatically an information-reducing map;
- a fixed finite set of quantized integral increments;
- the continuous mesoscopic Y and its own observation filtration;
- full q12 tagged-label histories, which are a richer candidate observation,
  not silently interchangeable with parity histories.

The first safe bridge is `N->infinity` at **fixed** finite observation resolution.
Only afterwards may resolution be refined. Allowing resolution to grow with N
or exchanging those limits needs a separate theorem. A continuously readable
exact microscopic integral is not a finite-resolution laboratory instrument.

For a minimal test, freeze lambda=0, k=1, horizon S=1, observations at scaled
times 1/4 and 1/2, and forecast target at time 1:

    C_N=sign(Z_N(1)), with a distinct ZERO outcome,
    V1=Y_N(1/4), V2=Y_N(1/2)-Y_N(1/4),
    Q(v) = LOW if v<=-1/2,
           MIDDLE if -1/2<v<=1/2,
           HIGH if v>1/2,
    R_N=(Q(V1),Q(V2)).

The target has three declared outcomes and the record nine. Retain all 27
joint cells; the limiting ZERO-target cells vanish but are not deleted from
finite-N bookkeeping. This is an analyst-chosen observable, not a new
microscopic axiom or an apparatus. Other choices require a new identifier,
not tuning thresholds after seeing whether the original one is informative.

No useful numeric N0 or finite-N rate is inherited. A convergent asymptotic
table is not a certified practical record requirement.

## 5. Task cards

### RI-001 — Bind record sigma-fields and an exact count/tag state

**Title:** One parent, distinct records, no hidden information upgrade.
**Scientific question:** what data are actually available in each claimed
observation, and which exact Markov state generates them?
**Why it matters:** a weak amplitude limit can lose statistically useful
information even when the underlying paths converge jointly.
**Exact inputs:** R01–R10, especially the actual labelled generator and R03.
**Dependencies:** admitted OBS result; no new microscopic transition.
**Method:**

1. Freeze the records, target, t_N, gain, clock, initial scope and horizon.
   Export separate sigma-field/observation IDs for parity path, integral
   path, quantized increments and any proposed full-label comparator.
2. Derive the exact finite-N Markov quotient `(n,b_1,...,b_k)`, where n is
   the total count and b lists tagged labels. A tagged i->j move changes n
   and that b coordinate with its actual refresh rate. An untagged i->j
   move changes n only, at `(n_i-c_i(b))*R_ij(n)`. State its admissible domain
   `n_i>=c_i(b)`; preserve self-event and clock conventions.
3. Prove lumping/intertwining directly. This is not the false claim that
   parity alone is exactly lumpable or that counts alone determine histories.
   The quotient may be large; derive its generator without enumerating it.
4. Show what the exact integral right derivative recovers, with endpoint/null
   jump conventions. Separate measurability from continuity of that inverse.
5. Reproduce the small logical example `C=+/-1, Y_N=C/N`: weak convergence
   to Y=0 loses a perfectly readable bit. It is a counterexample to an
   inference rule, not a no-go theorem for FIN itself.
6. Test tiny exact support fixtures with the corrected helper, including
   `(N,M,k1,k2,h1,h2)=(12,2,4,3,3,0)` and degenerate parity populations.

**Deliverables:** RECORD_CONTRACT.json, EXACT_TAG_COUNT_GENERATOR.md,
SIGMA_FIELD_MAP.md, INFORMATION_LIMIT_COUNTEREXAMPLE.md, INPUTS.sha256.
**Acceptance criterion:** one reproducible parent and explicit observation
maps; no unsupported sufficiency, independence or recoverability claim.
**Failure/refutation criterion:** a claimed prediction needs hidden tagged
data, invisible swaps, a new sampling law, or an unprovided microscopic lift.
**What becomes possible:** a well-defined predictive statistical experiment.
**Must NOT conclude:** a physical observer/storage medium, spatial factors,
microscopic information conservation or a full-record diffusion equivalence.
**Estimated computational class:** S/M; finite algebra, no large enumeration.
**Priority:** P0.

### RI-002 — A finite-resolution predictive experiment with a real risk limit

**Title:** Convert a weak record limit into a bounded prediction theorem.
**Scientific question:** what can the frozen nine-outcome past record predict
about the later collective sign, and does that result transfer to finite FIN?
**Why it matters:** a joint finite-alphabet probability table supports an
actual risk theorem without illicit entropy or filtration convergence.
**Exact inputs:** RI-001, R02/R03/R07/R08 and section 4's frozen experiment.
**Dependencies:** correct record/target maps and actual iid path transfer.
**Method:**

1. Prove every quantizer boundary and target-zero event has zero limiting
   probability. Nondegenerate independent observation noise supplies continuous
   increment densities; the stationary quartic state has no sign atom.
2. Deduce convergence of all 27 cell probabilities and hence joint TV on
   this fixed finite alphabet. Retain the MT/path-exit o(1) errors.
3. Prove convergence of the finite-experiment Bayes risk
   `r_N=1-sum_r max_c P(C_N=c,R_N=r)` using that finite joint-TV convergence.
   This is optimal only among rules using the frozen finite record.
4. Compute or bound that one limiting table by a justified small diffusion/
   additive-functional method. If using a truncated PDE/quadrature, pay its
   domain, discretization and tail errors. Pure Monte Carlo is diagnostic.
5. Test informativeness against the no-record symmetric-sign baseline 1/2.
   A positive numerical advantage without an error payment is not a certified
   strict inequality. If the coarse statistic is inconclusive, retain it;
   do not shift thresholds or readout times under the old experiment ID.
6. Distinguish prediction of the later sign from recovery of every state,
   permanent sign memory and a physical detector performance claim.

**Deliverables:** FIXED_EXPERIMENT.json, FINITE_ALPHABET_LIMIT_PROOF.md,
JOINT_TABLE_OR_BOUNDS.json, BAYES_RISK_AND_ERROR_SCOPE.md, NONCONCLUSIONS.md.
**Acceptance criterion:** a genuine finite-outcome/risk convergence theorem,
with numerical values separately identified by their arithmetic status.
**Failure/refutation criterion:** a boundary atom is ignored, a conditional
posterior limit is inferred from weak convergence alone, or the record has
no paid predictive advantage. A correct no-advantage result is admissible.
**What becomes possible:** a bounded intrinsic observational prediction of
the same microscopic FIN, with no new law or external experiment.
**Must NOT conclude:** full-record optimality, practical finite-N sample
complexity, detector calibration, natural A7 selection or empirical validation.
**Estimated computational class:** M; one small fixed experiment only.
**Priority:** P0.

### RI-003 — Derive the causal filter of the limiting model, not a new law

**Title:** Conditional state prediction from the admitted continuous record.
**Scientific question:** what causal conditional law is determined by the
already derived joint classical SDE?
**Why it matters:** independent driving noises do not make state and record
independent, and the microscopic noise W is not the observation innovation.
**Exact inputs:** RI-001, OBS joint theorem, stationary quartic moments.
**Dependencies:** correct limiting law and natural observation filtration;
RI-002 may proceed independently after RI-001.
**Method:**

1. Define `pi_s(f)=E[f(Z_s)|sigma(Y_r:0<=r<=s)]`. Derive, rather than merely
   name, the candidate equation for `Gf=(chi lambda z-z^3/3)f'+f''`:

       d pi(f)=pi(Gf)ds
          +(k/2)[pi(zf)-pi(z)pi(f)] [dY-pi(z)ds].

2. Pay the change-of-measure/localization and integrability conditions for
   the unbounded signal z and cubic drift. The quartic stationary tails and
   finite-horizon exponential quadratic moments are relevant inputs, not
   permission to invoke every bounded-coefficient theorem unchanged.
3. Prove that `I=sqrt(k/2)(Y-integral pi(z)ds)` is Brownian in the observation
   filtration. Do not identify I with the hidden W or make it independent of
   the signal merely because it is an innovation.
4. Give a well-defined conditional-law solution/uniqueness statement in the
   declared class. A Gaussian/Kalman closure is not exact for quartic Z.
5. Recover the finite-record probabilities as conditional expectations under
   the same joint law. Keep approximate filter numerics separate from the
   exact equation, with errors if numerical claims are made.
6. Read any external theorem actually used and bind its assumptions/source.
   Do not turn a literature analogy into a physical premise.

**Deliverables:** LIMIT_FILTER_THEOREM.md, INNOVATION_AND_MOMENT_PAYMENT.md,
FILTER_SCOPE.json, REPLAY.md, or an exact missing lemma if proof fails.
**Acceptance criterion:** a paid filter for the *limiting mathematical*
experiment, without declaring the finite-N microscopic filter convergent.
**Failure/refutation criterion:** closure needs fitted noise, discarded higher
moments, future data, an unlicensed filtration, or an unpaid theorem hypothesis.
**What becomes possible:** rigorous causal inference within the mesoscopic
model and a clearly stated target for the microscopic sufficiency test.
**Must NOT conclude:** a new physical measurement postulate, fundamental
observer, Born rule, information conservation or autonomous memory hardware.
**Estimated computational class:** M; analytical filtering, not a model sweep.
**Priority:** P0.

### RI-004 — Full-record predictive sufficiency: prove it or isolate its failure

**Title:** Test information lost by the diffusion limit.
**Scientific question:** is observing the entire finite-N parity history
asymptotically equivalent, for the chosen bounded target, to its continuous
mesoscopic record experiment?
**Why it matters:** this is the bridge from correct coarse trajectories to
a genuinely reusable predictive observation theory.
**Exact inputs:** RI-001–003, exact tagged/count generator and parent couplings.
**Dependencies:** frozen observation scope and a paid limiting filter or
an explicitly narrower finite-resolution comparison target.
**Method:**

1. Specify one precise comparison: risks for the same bounded sign target,
   or convergence of stated bounded posterior functionals. Do not invoke an
   undefined universal Le Cam equivalence over an unspecified model family.
2. Write the exact observation likelihood/filter for the finite-N tagged/count
   model, marginalizing hidden states and preserving their preparation law.
   A prescribed deterministic external signal is not the interacting parent.
3. Attempt a compact-localized likelihood/martingale comparison. Pay the
   influence of tagged copies on the collective state, fast occupation
   corrections, rare exits and UI of likelihood ratios, not just their logs.
4. Identify any residual informative jump/mark/timing statistic omitted by
   the coarse observation. Full q12 label changes and refresh-attempt marks
   are *different* data contracts; do not slip them into a parity observer.
5. Keep `N->infinity` at fixed finite resolution separate from increasing
   temporal/quantization resolution. A lower bound for full-record Bayes risk
   needs its own proof; finite-bin convergence gives only an available rule.
6. If the comparison fails, give one genuine witness of extra predictive
   information or one precise unpaid likelihood lemma. Do not repair the
   microscopic law or assert information equivalence from excellent plots.

**Deliverables:** MICROSCOPIC_LIKELIHOOD_CONTRACT.md,
SUFFICIENCY_THEOREM_OR_OPEN_ATOM.md, LIMIT_ORDER_LEDGER.json, NEGATIVE_RESULTS.md.
**Acceptance criterion:** a proof, a scoped counterexample, or an explicit
open estimate that prevents promotion of full-record sufficiency. Preserve
the narrower RI-002 theorem even if this stronger target remains open.
**Failure/refutation criterion:** finite-N records improve bounded prediction
by a nonvanishing amount over the proposed limit, or the necessary likelihood
control cannot be paid. Only the claimed sufficiency is killed, not OBS.
**What becomes possible:** an honest sufficient observational effective law,
or identification of the exact extra retained statistic needed.
**Must NOT conclude:** entropy/MI convergence without separate bounds,
quantum information, irreversible collapse or a physical detector bridge.
**Estimated computational class:** M, hard proof gate; no brute-force escalation.
**Priority:** P1, after the three P0 tasks.

### RI-005 — Scope synthesis and architect return

**Title:** Which predictive record theory is actually supported?
**Scientific question:** what part of state evolution, observation and
inference now follows from the same FIN parent without a new primitive?
**Why it matters:** accurate records, information retention, irreversible
storage and physical observation are different claims.
**Exact inputs:** all RI results/failures and R01–R10.
**Dependencies:** terminal scoped dispositions, not forced success of RI-004.
**Method:**

1. Export the DAG: supplied law/preparation -> local window -> joint record
   -> finite-resolution probabilities -> causal filter -> sufficiency gate.
2. Distinguish reused OBS/MT/UF results, exact finite identities, analytical
   limits, bounded risk theorems, numerical tables and unresolved estimates.
3. Reuse the admitted equilibrium reversal corollary. Under reversal of
   increments `RY(s)=Y(S)-Y(S-s)`, a stationary record does not generate an
   intrinsic thermodynamic arrow. No entropy-production convergence follows
   from an o(1) TV comparison alone.
4. State which bounded predictions are informative and which apparent
   fingerprints disappear under canonical normalization. No A7 uniqueness
   follows from class-level critical behavior.
5. Return one next open atom and stop. Do not append a drive, medium,
   conserved charge, geometry or memory cost to obtain the desired outcome.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
NEXT_OPEN_ATOM.md, REPLAY.md, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** one reproducible, evidence-scoped handoff with all
unpaid likelihood/physical claims unpromoted.
**Failure/refutation criterion:** mathematical integration is called storage,
filtering is called a physical observer, or a stationary classical channel
is presented as an emergent arrow/Born rule.
**What becomes possible:** an architect decision based on genuinely derived
observational content rather than a new assumption layer.
**Must NOT conclude:** physical source, units/locality/conservation/EV,
selector/role transfer, FIN-M0/B4/task335, QM/GR/SM/L_total or ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution and stop rules

Start RI-001, then RI-002 and RI-003. RI-004 is one separate, bounded proof
attempt; a precise remaining lemma is preferable to hundreds of new tasks.
Return RI-005 for architect selection. Do not run subsequent physics campaigns.

- Work only under the declared new execution directory. Preserve AGENTS,
  prior sources and reviews; return a proposed patch, not automatic promotion.
- No large count/tag state enumeration, broad N/g or threshold optimization,
  new noise/controller, physical graph, hardware, spending or messages.
- Tiny combinatorial/generator tests and one low-dimensional limiting-law
  computation are sufficient. Numerical precision escalation is not a proof.
- Per task export REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,
  NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. Separate execution/scientific
  status, supplied/derived/fitted objects, model/record/filter IDs, initial
  and time scope, arithmetic assumptions, errors, failures and next action.
- A theorem invocation requires the actual statement and hypothesis payment.
  A weak convergence theorem is not an observation-experiment theorem.
- All physical source, ontology, strict/legacy, selector, B4/FIN-M0/task335
  and empirical activation gates remain binding. No external data is invented.

## 7. Kill tests and what would count as progress

Kill only the claim that fails: a wrong record sigma-field, a boundary atom,
no informative finite statistic, an unpaid filter theorem, or a nonvanishing
full-record prediction advantage over the proposed coarse experiment.
Do not kill the correct microscopic/state/OBS theorem because a stronger
information claim fails. Do not evade a failure by changing the law.

Success would establish a controlled classical *predictive observation*
layer of the supplied statistical model. It would not derive nature's A7,
its kinetics or clock, physical spacetime, conserved energy, particles,
quantum theory or gravity. Those are separate gaps, not conclusions obtained
by renaming the state, record or conditional probability.
