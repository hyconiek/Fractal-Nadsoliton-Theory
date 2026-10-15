# FIN after LI: universality versus kernel-sensitive intrinsic observables

Date: 2026-10-15. Repository baseline: `380e9dc0` plus the LI intake.
Status: **RESEARCH_PLAN; KO-001–KO-005 have not been executed.**
Execution directory: `fin_kernel_observability_campaign_20261015/`.
Language: English for the execution agent.
Decision: [KO_DECISION.md](fin_parity_record_sufficiency_review_20261015/KO_DECISION.md).

## 1. The next scientific question

The local critical observation and one-target full-parity risk bridge is
sufficiently closed. Do not make further scalar precision or more tagged
variants the default programme. Ask instead:

> Which intrinsic observations distinguish A7 from other matched operators,
> and which successful predictions belong to a broader critical class?

This is a theoretical identifiability/universality test. It neither derives
a kernel from first principles nor licenses a physical experiment. Primary
FIN remains unchanged. Alternate operators below are explicit diagnostic
controls, not newly adopted microscopic laws or ontology.

## 2. Exact inputs and reused results

| Alias | Repository input |
|---|---|
| K01 | AGENTS.md, latest LI/RI/OBS/MT/UF guardrails |
| K02 | fin_parity_record_sufficiency_review_20261015/INTAKE_20261015.md, PROOF_COMPLETIONS.md, verification.json and CLAIM_REGISTER.json |
| K03 | fin_parity_record_sufficiency_review_20261015/KERNEL_FAMILY_SEED.md, KERNEL_FAMILY_SEED.json and kernel_family_certificate.py |
| K04 | fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md and INTAKE_20261015.md; fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md |
| K05 | fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md and UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md |
| K06 | fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md; fin_intrinsic_observation_campaign_20261015/OBS-C01_LABELLED_TRANSFER.md |
| K07 | fin_record_information_review_20261015/PROOF_COMPLETIONS.md and TAG_FEEDBACK_REMOVAL_LEMMA.md |
| K08 | fin_generative_operational_campaign_20261014/RG-E02/IDENTIFIABILITY_THEOREM.md; fin_generative_operational_review_20261014/INTAKE_20261014.md |
| K09 | fin_theory_synthesis_20261014/GEN_EVIDENCE_MAP.md and CORE_DERIVED_BRIDGES.md; fin_intrinsic_theory_review_20261014/INTAKE_20261014.md |
| K10 | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py; fin_rank7_followup/inputs/fin_handoff_audit/results.json |

Hash actual files. Existing static fast covariance, leading drift/noise,
stationary hidden OU and Fourier identifiability results are **reused**, not
new discoveries. The new scope is the critical local window, a changed
partial readout and matched-operator comparison. If an earlier result already
settles a substep in exactly that scope, bind it instead of rerunning it.

## 3. Frozen family and time conventions

Use q12 and the original labelled unit-attempt LOO heat bath. For primary
A7 retain exact g_c=12/Lambda6, iid-uniform initial law and the admitted
deterministic post-MT local window. N passes through multiples of12.

Comparison seed, already algebraically paid under inherited spectral boxes:

    A_epsilon=A7+epsilon(P3-P4), |epsilon|<=1/100.

P3,P4 are the rank-two real Fourier projectors. Trace, diagonal, rank,
centering, D12 covariance and Lambda6 are unchanged. The active eigenvalues
Lambda3 and Lambda4 shift by +epsilon and -epsilon. First compare only
epsilon=0 and+1/100. A continuous-family statement is optional and must state
pointwise versus uniform parameter scope.

The parity line and g_c are identical algebraically. **Transfer of the full
process/local-regime theorems is still a KO-001 proof obligation**, not a
consequence of the seed file's PASS flag. Do not borrow A7's full global
phase census, coexistence point, domain constants or B4 claims for a control.

Freeze a common deterministic observation origin after the admitted local
relaxation scales of both models, with future finite windows still before
their paid exits. This is an existence-level mathematical contract, not
a numeric burn-in or physical clock calibration.

The two analysis time scales must remain distinct:

- critical soft time: model t=sqrt(N)s;
- fast readout time: model t=O(1) after that same origin.

No physical length or dispersion relation is introduced by Fourier labels.

## 4. The new snapshot observation

Choose a predetermined block B_N of k_N=N/2 slots, independently of all
states/paths. Observe its **q12 labels at one time**, not only one tag's
parity. This is a new readout ID, not free information given to the LI observer.

For r=1,3,4 let `f_r(j)=sqrt(2) cos(2*pi*r*j/12)`, so its uniform mean is0
and variance1. Define

    F_r,N=sqrt(N) sum_j f_r(j) p_j,
    V_r,N=sqrt(N) (1/k_N) sum_(a in B_N) f_r(xi_a),
    U_r,N=V_r,N-F_r,N,
    gamma_r=1-Lambda_r/Lambda6 for r=3,4; gamma_1=1.

Lambda1 here is **zero in A7**, not the nonzero k1 eigenvalue of the original
Laplacian from which the retained construction was selected. Do not confuse
those spectra. All f_r are orthogonal to parity, so no soft mean is subtracted
by an estimated or fitted controller.

The candidate static joint law is

    F_r ~ Normal(0,1/gamma_r),
    U_r ~ Normal(0,c_alpha), c_alpha=(1-alpha)/alpha=1,
    V_r=F_r+U_r,

with Gaussian mode/private independence in the stated limit and independence
from the leading quartic soft coordinate. MT already gives the global
Gaussian covariance; the macroscopic sampling channel needs a separate
multivariate hypergeometric proof.

Do not reuse the small-block collision bound k^2/N here: k=N/2 is a different
regime. No finite-N independence is claimed.

## 5. Task cards

### KO-001 — Matched operators and exactly scoped parity degeneracy

**Title:** Separate universal leading parity behavior from A7-specific content.
**Scientific question:** do the two distinct matched operators share the
same leading critical parity law and fixed-resolution parity predictions?
**Why it matters:** a successful critical observation theorem cannot select
the kernel if the same theorem follows for a nontrivial comparator.
**Exact inputs:** K01–K10, especially the exact family certificate.
**Dependencies:** admitted LI chain; no new microscopic FIN premise.
**Method:**

1. Recheck rational projector identities and bind the exact spectral boxes.
   Keep primary A7 and each control under separate model IDs.
2. Audit precisely which UF/MT/OBS/LI hypotheses use the top parity mode,
   positive transverse gap, symmetry, smoothness, PSD and clock. Verify these
   for epsilon0 and+1/100; do not just copy a result's status string.
3. Pay any local-domain/trapping/initial-law extension needed by that audit.
   Existence of a sufficiently small domain may be enough; certified numeric
   N0 and a global saddle cover are not required. Record shared versus
   model-specific analysis domains and the common observation-origin scope.
4. Establish the common leading quartic/record law only where hypotheses
   are paid. For a fixed finite quantized parity record, show distribution
   convergence to the same finite table and its resulting ideal binary
   model-identification risk ->1/2.
5. **Do not** infer full-parity model-identification blindness from LI's
   one-target risk theorem. Distinguishing which operator generated data is
   a different target. Subleading/fine-timestamp information remains outside
   this finite-resolution no-go unless separately proved.

   A concrete optional stronger proof uses a fair **analysis model index**
   B selecting the two fixed operators, model-conditional autonomous
   backgrounds, and a common fair parity-path reference independent of both.
   Prove joint convergence `(B,Z_N,Y_N,L_N^B)->(B,Z,Y,L)` with B independent
   of the limit paths/likelihood, plus mean-one UI. Then apply the already
   paid LI finite-cylinder density argument with target B, not future sign.
   If every premise is paid, this would imply full-parity model-ID risk->1/2
   for the two frozen controls. The model prior is a testing convention,
   not a new FIN law. If this target extension is not proved, label it
   OPEN_CROSS_MODEL and retain only the finite-resolution conclusion.

6. Deduplicate old OU/CLT/Fourier identifiability results. State precisely
   which new observation scope remains before computing it.

**Deliverables:** MATCHED_MODEL_CONTRACT.json, STRUCTURAL_HYPOTHESIS_AUDIT.md,
LEADING_PARITY_SCOPE.md, DEDUPLICATION_LEDGER.md, INPUTS.sha256.
**Acceptance criterion:** a matched family with a paid leading/fixed-resolution
comparison, or the exact missing structural premise; no full-record no-go by analogy.
**Failure/refutation criterion:** the transverse gap, local window or symmetry
fails, or a supposedly shared leading law genuinely changes. Keep that result.
**What becomes possible:** a precise universality boundary and a meaningful
search for kernel-sensitive observations already present in the same model.
**Must NOT conclude:** necessity/nonnecessity of all FIN physics, natural A7
source, full statistical equivalence or unchanged global phase structure.
**Estimated computational class:** S/M, proof audit and tiny exact certificates.
**Priority:** P0.

### KO-002 — Derive the half-population fast snapshot channel

**Title:** Read existing fast fluctuations with their true sampling noise.
**Scientific question:** can a declared partial snapshot retain the fast
spectral information absent from the leading parity channel?
**Why it matters:** the new observable must be derived, not supplied with a
noise variance chosen to reveal A7.
**Exact inputs:** KO-001, K04/K06/K10, section 4's fixed readout.
**Dependencies:** local Gibbs/actual-parent transfer for the compared models.
**Method:**

1. Write the exact multivariate count-hypergeometric law of block labels
   conditional on the complete population count. Keep integer support,
   finite-population correction and all outcomes.
2. Apply MT's global joint quartic/fast Gaussian limit in the correct f_r
   normalization: E_u f_r^2=1 gives Var(F_r)=1/gamma_r, not the Euclidean
   occupation-coordinate precision inverse without its factor12.
3. Prove the conditional multivariate sampling CLT for k/N->alpha=1/2.
   Derive Var(U_r)->c_alpha and the cross covariances. Joint Gaussianity and
   conditional characteristic-function convergence, not zero covariance
   alone, must pay the claimed independent limiting errors.
4. Transfer weak laws to the actual iid-start observation by labelled TV.
   Do not move unbounded scaled moments through that error; bound decisions
   directly or prove moments natively in the stationary comparison.
5. Verify internal null-mode r1 gives Var(V1)=1/alpha=2. Compare r3,r4
   without fitting their noise, gain, block fraction or effective clock.
6. If claiming a minimum block size, prove that stronger resource claim.
   This task only requires a sufficient N/2 observation, not optimality.

**Deliverables:** EXACT_FAST_READOUT_KERNEL.md, JOINT_SNAPSHOT_LIMIT.md,
NORMALIZATION_AND_ERROR_LEDGER.json, SMALL_SUPPORT_CHECKS.json.
**Acceptance criterion:** a joint snapshot law with paid sampling/conditioning
scope and unchanged parent, or a precise failed CLT/moment step.
**Failure/refutation criterion:** sampling noise or cross dependence survives
with a different value, the readout secretly uses full counts, or labels are
chosen adaptively after seeing the state.
**What becomes possible:** a bounded ideal test of dimensionless kernel content.
**Must NOT conclude:** a physical half-universe instrument, spatial subsystem,
entanglement, fundamental conserved information or an independently sourced kernel.
**Estimated computational class:** M, combinatorics and local measure limits.
**Priority:** P0.

### KO-003 — Freeze one bounded ideal kernel contrast

**Title:** Make the spectral distinction an observable probability difference.
**Scientific question:** does one fixed event distinguish the two matched
operators in the derived readout limit?
**Why it matters:** covariance formulas alone are unbounded statistics;
a bounded event transfers through the actual-parent TV comparison cleanly.
**Exact inputs:** KO-001–002; epsilon0 versus+1/100; alpha=1/2.
**Dependencies:** accepted fast snapshot channel for both models.
**Method:**

1. Freeze `E_N={|V_3,N|>1}` before numerical output. This uses the new
   half-population q12 snapshot, not the old one-tag parity observer.
2. Derive `p_epsilon=2 Phi(-1/sqrt(1/gamma_3(epsilon)+1))`. Prove the null
   threshold boundary has zero limiting probability and P(E_N)->p_epsilon.
3. Since Lambda3 increases but Lambda6 is unchanged, prove strictly
   p_(+1/100)>p_0. Exact monotonicity suffices; do not make additional digits
   the campaign objective. Any decimal enclosure needs its actual arithmetic proof.
4. Give the limiting two-outcome law TV gap and the equal-prior ideal binary
   decision-risk formula `(1-|p_epsilon-p_0|)/2`. It is a model discrimination
   problem, not a forecast of each model's own future sign.
5. Separate this perfect-readout mathematical contract from RG-F's broad
   calibrated uncertainty classes. This does not reopen its binary-event
   overlap no-go, license cheap laboratory discrimination or supply calibration.
6. Report what ratios/coefficients are model-specific, symmetry-protected
   or class-level. Inferring a spectral ratio does not explain why nature
   selected it, rank7, g_c or the microscopic heat bath.

**Deliverables:** FROZEN_KERNEL_EVENT.json, BOUNDED_CONTRAST_THEOREM.md,
DIMENSIONLESS_FINGERPRINTS.md, NONCONCLUSIONS.md.
**Acceptance criterion:** a strict bounded ideal separation or a retained
failure, with no parameter/threshold/readout retuning under the old ID.
**Failure/refutation criterion:** the event laws coincide, the normalization
changes between models or success needs an unprovided instrument premise.
**What becomes possible:** an intrinsic kernel-sensitive prediction distinguished
from the universal critical mechanism, without claiming empirical validation.
**Must NOT conclude:** physical A7 uniqueness, natural constants, real detector
performance, finite-N sample complexity or a physical-source theorem.
**Estimated computational class:** S/M, analytical contrast; tiny numerical checks only.
**Priority:** P0.

### KO-004 — Conditional fast-time process gate, not another default campaign

**Title:** Derive a two-scale observed process only if a dynamic contrast is needed.
**Scientific question:** does the new readout admit the predicted sum of
global and private relaxation modes on O(1) model time?
**Why it matters:** a static Gaussian channel is not automatically its process,
and Fourier mode numbers are not physical wavevectors.
**Exact inputs:** KO-001–003, K04–K07 and deduplicated existing drift/CLT results.
**Dependencies:** terminal static kernel contrast and an explicit need for one
dynamic observable; otherwise return DEFERRED_NOT_REQUIRED and continue KO-005.
**Method:**

1. Work on one fixed fast-time horizon after the same origin. The critical
   Z coordinate is frozen in this limit, not dynamically replaced by zero.
2. Reuse the full fast linearization and MT Gaussian initial law. Pay their
   critical-preparation scope if an earlier stable-phase theorem excludes u.
3. Generalize the actual HB block-minus-population contrast identity to f_r.
   Derive the private drift/noise; do not reuse parity's pointwise constant
   jump bracket, which is false for general Fourier f_r.
4. Pay integrated global/private cross brackets by conditional block sampling
   and stationarity/coupling. Derive or refute independent Gaussian drivers
   and initial errors in the candidate joint system:

       dF_r=-gamma_r F_r dt+sqrt(2)dB_r,
       dU_r=-U_r dt+sqrt(2 c_alpha)dW_r,
       V_r=F_r+U_r.

5. Only if proved, export
   `Cov(V_r(t),V_r(0))=gamma_r^(-1)e^(-gamma_r|t|)+c_alpha e^(-|t|)`.
   V alone is not generally Markov. Use gamma1=1 as a model-clock control;
   rate ratios can remove a common clock gauge, not create physical seconds.
6. Freeze at most one two-time contrast. No new pole-fit sweep, white noise,
   spatial gradient or full critical-horizon fast-path limit is needed.

**Deliverables:** FAST_PATH_THEOREM_OR_GATE.md, BRACKET_AND_INITIAL_PAYMENT.md,
ONE_DYNAMIC_CONTRAST.json or DEFERRED.md.
**Acceptance criterion:** one paid joint fast process/two-time contrast or
a clearly deferred/open gate; no process inferred solely from covariance.
**Failure/refutation criterion:** nonvanishing coupling changes the process,
brackets are unpaid, or a known leading OU law is relabelled a new physical result.
**What becomes possible:** a process-level kernel-sensitive effective prediction.
**Must NOT conclude:** field propagation, physical mass gap, metric/space,
particle modes, conserved energy or a new fundamental clock.
**Estimated computational class:** M, proof gate; no large generator enumeration.
**Priority:** P1, conditional.

### KO-005 — Scope synthesis and architect stop

**Title:** What is universal, what identifies the model, and what remains physics?
**Scientific question:** does FIN now have kernel-sensitive intrinsic predictions
beyond its common critical observation law?
**Why it matters:** successful statistical reduction must not be confused with
selection of fundamental laws or empirical evidence.
**Exact inputs:** all KO results/failures and K01–K10.
**Dependencies:** terminal scoped dispositions; KO-004 may be explicitly deferred.
**Method:**

1. Export a DAG and a table of supplied/derived/fitted objects. Separate
   leading parity degeneracy, finite-resolution discrimination, fast snapshots,
   optional process laws and actual physical-source gaps.
2. Retain all failed contrasts and unresolved theorem hypotheses. Do not
   infer full parity model-ID blindness from LI's within-model sign-risk result.
3. Identify reused Gaussian/OU/static covariance and RG-E Fourier algebra.
   No credit for rediscovering old formulas without a new scope theorem.
4. State the remaining source/units/locality/conservation/empirical bridges.
   Their absence is not repaired by calling internal Fourier modes fields.
5. Return one evidence-based next question and stop. No automatic larger
   model, extra microscopic primitive, closed-source inventory, entropy/QM/GR
   programme or physical experiment follows.

**Deliverables:** HANDOFF.md, REPORT.md, CLAIM_REGISTER.json, NONCONCLUSIONS.md,
NEXT_OPEN_ATOM.md, REPLAY.md, MANIFEST.sha256 and proposed AGENTS patch.
**Acceptance criterion:** reproducible scoped comparison, including a correct
negative result or deferred process gate, with physical claims unpromoted.
**Failure/refutation criterion:** operator identification is called source
derivation, sampling tags become physical cells, or a control silently replaces A7.
**What becomes possible:** one architect decision about genuine model-specific
predictive content rather than more scalar precision.
**Must NOT conclude:** kernel/kinetic/gain/initial-source, physical units/
locality/conservation/EV, selector/role transfer, FIN-M0/B4/task335,
QM/GR/SM/L_total or ToE closure.
**Estimated computational class:** S.
**Priority:** P0.

## 6. Execution, stop conditions and scientific value

Start **KO-001–003**. KO-004 needs its stated activation; it does not block
the static bounded theorem. Return KO-005. Work only in the new directory,
preserve AGENTS and all historical packages, and return a proposed patch.

Per task export REPORT.md, RESULTS.json, PROOF_OR_OPEN_LEMMA.md,
NONCONCLUSIONS.md, REPLAY.md and INPUTS.sha256. Separate execution/scientific
status, model/readout IDs, source/reused/new objects, initial/time/parameter
scopes, arithmetic assumptions, errors, failures and next action.

No broad N/g/epsilon/threshold sweep, large count enumeration, new graph,
noise/drive/storage law, B4/FIN-M0/task335 activation, external experiment,
purchase or messaging. Read and pay any literature theorem actually used.

Kill the relevant interpretation if the matched-family hypotheses fail,
the sampling law changes the claimed contrast, or ideal observation cannot
distinguish the proposed models. Keep a correct partial theorem. Success
supports an intrinsic **model-specific statistical prediction**, not nature's
selection of the supplied operator or a fundamental theory of reality.
