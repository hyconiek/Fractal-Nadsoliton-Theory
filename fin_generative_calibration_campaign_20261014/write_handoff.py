#!/usr/bin/env python3
"""Generate human reports from frozen RG-F computational output; never changes scientific results."""
import json,hashlib,platform,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def read(p):return json.loads((ROOT/p).read_text())
def write(p,text):
 x=ROOT/p;x.parent.mkdir(parents=True,exist_ok=True);x.write_text(text.strip()+'\n',encoding='utf-8')
def jwrite(p,obj):write(p,json.dumps(obj,indent=2,sort_keys=True))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

proc=read('RG-F01/COMPOSITION_TESTS.json');sym=read('RG-F02/SYMMETRY_AND_ERASURE_CONTROLS.json')
stat=read('RG-F03/TOTAL_LAW_ERROR_BUDGET.json');decision=read('RG-F04/DECISION_RULE.json')
verify=read('INDEPENDENT_VERIFICATION.json');synt=read('RG-F04/SYNTHETIC_COMPOSITE_TESTS.json')

write('RG-F01/PROBABILITY_SCOPE.md',r'''# Probability scope — RG-F01

**Exact finite-model result.** Let `R[x,d]=1{d=(b-a) mod12}`, for `x=(a,b)` on 144 labelled states. For every specified zero-field circulant, equivariant generator `L`, the rows satisfy `L R = R B`; `B` has nonnegative off-diagonal rates and zero row sums. Therefore a single 12-state CTMC works for **every microscopic initial phase distribution**, after pushing it to the difference marginal. The heat-bath `B=2(1 p^T-I)` and `exp(tB)=exp(-2t)I+(1-exp(-2t))1 p^T`; the other two B matrices are distinct. The process remains exact only for this difference-fibre model.

**Exact instrument rule.** A fixed column-stochastic 13×12 matrix `K[o,d]` defines a nondisturbing classical emission (including `INVALID`) with branch update `v -> v diag(K[o,:])` after propagation. Branches stay unnormalized. Summing over all outcomes gives the propagated law; summing intermediate observations recovers the shorter record law. This constructs one positive, compositional likelihood for arbitrary finite sequences and mixed initial laws. No basis/parameter is refitted per record.

**Physical assumption.** A detector that disturbs the state needs an independently calibrated family `J_o[d,d']` (or a full 144-state successor kernel). `K` is not sufficient to represent disturbance or an invalid physical state that evolves differently. An observed noisy output string need not itself be a Markov chain. Phase-dependent emissions may destroy the 12-state instrument representation even though the parent kinetics remain lumpable.

**Arithmetic.** Serialized matrices and laws are double-precision numerical witnesses; positivity of the abstract CTMC is an exact mathematical consequence of the nonnegative rate definition. No full-operation interval certificate is claimed.
''')
write('RG-F02/QUOTIENT_VALIDITY_TEST.md',f'''# Quotient validity — RG-F02

The testable exact condition is equality of **whole fibre transition rows** for all phases `(a,b)` having the same `d=(b-a) mod12`. The ideal finite generators pass numerically to residual {sym['mathematical_zero_field_quotient_defect']:.3e}. This is a model verification, not a physical symmetry measurement.

The single declared negative control adds an absolute first-copy preference `h=+0.01` for label 0. It creates a maximum same-fibre transition-rate defect of **{sym['one_declared_0p01_absolute_first_copy_field']['max_block_rate_difference']:.12f}** per model time. The producer computes all 144 rows rather than assuming averaged lumpability. Do not generalize to other asymmetries.

A second counterexample emits a 1% invalid outcome only for phase `a=0`. It fails phase-independent 13×12 emission despite unchanged ideal dynamics. The **erasure channel** emits each of 12 valid values uniformly, never invalid; it passes the uniform g=0 valid-output histogram and destroys every g=3 interaction fingerprint. Thus a uniform static control alone is insufficient.

**Independent physical calibration required:** prepare phase-distinct absolute pairs at equal difference; retain raw absolute states, timestamps, controls and readouts; compare transitions into each of 12 difference fibres at a common calibrated rate; compare **joint outcome/successor** laws for phase shifts, including invalids and disturbance. If impossible, retain the 144-state physical parent or report an independently bounded defect. No physical observations have been provided.
''')
write('RG-F02/FALLBACK_MODEL_GATE.md',r'''# 12-state or 144-state fallback gate

1. If independent phase-resolved calibration confirms or rigorously bounds full fibre-generator equality and the complete joint instrument factors through difference, use the 12-state quotient, with any validated defect charged to the total law error budget.
2. If rates vary with absolute phase, phase conditioning/readout varies, or invalid records cause hidden-state feedback, retain the full 144-state parent (or a justified augmented physical state). Never enforce exact lumpability by stationary averaging.
3. If physical calibration lacks absolute phase access, label the quotient permission `OPEN_INPUT`, not PASS. Difference-only records cannot prove all 144 transition rows identical.
''')
write('RG-F03/CALIBRATION_TO_MODEL_MAP.md',r'''# Calibration-to-model map — RG-F03

The already reviewed **hypothetical** interaction comparison uses FIN and flat-P7, both at independently varying gain `g in [2.97,3.03]`, with at most `0.005` output-TV deviation **per hypothesis**. The predecessor intake supplied a complete-gain-cover lower bound `0.01119114` before output allowances; hence `0.00119114` after both `0.005` allowances. This F package inherits that conditional result; it does not contain the repository review's exact-fraction gain-cover witness and must not be cited as a new independent interval replication.

Additional unaccounted errors `e_FIN+e_flat` shrink the guaranteed gap to `max(0, 0.00119114-e_FIN-e_flat)`. All unbounded physical errors make the actual gate OPEN. The illustrative allocations in the JSON are not measured. Do not count the same detector uncertainty twice: allocations labelled additional must be disjoint from the already-budgeted `0.005` class.

For a heat-bath reset process, `mu(t)=exp(-2t)mu(0)+(1-exp(-2t))p`. A worst-case equilibrium error `eps_burn` requires `t >= (1/2)log(1/eps_burn)` **in model time**. Converting this to seconds needs an independent measured production clock and a check that the attempt rate is stable when coupling is enabled. Do not treat the earlier `t=.4/1.0` quench as stationary.

For two finite CTMCs on the SAME labelled state space with generators `L,L~` and identical readout semantics, Duhamel and TV contraction give `TV(mu exp(tL), mu~ exp(tL~)) <= TV(mu,mu~) + (t/2) max_x sum_y |L_xy-L~_xy|`. Concatenate this term across intervals; at each observation add the supremum TV over fixed outcome/**successor** kernels `J_o`. Stochastic maps contract TV. This is an unconditional algebraic bound **if** its numerical inputs and representation premises are supplied. A postselected probability requires separately bounded selection mass.

Uniformization `P=I+B/4` is stochastic for all nine nominal 12-state generators and offers positive finite propagation. The Poisson tail bound in the JSON is a high-precision diagnostic; matrix elements and errors of the actual compiled rates still require an independent interval enclosure. No all-operations outward floating-point or physical prediction enclosure is claimed.
''')
write('RG-F03/SEPARATION_GATE.md',r'''# Separation gate

**Mathematical conditional example:** inherited rho >= 0.00119114 under the stated gain/readout class; positive only if every additional, disjoint pairwise TV error sums to strictly less than rho. The sample confidence radius must also satisfy 2*epsilon < the remaining gap. A common fixed detector has a different model-class geometry than two independently variable nuisance detector laws; the latter is conservatively used here.

**Empirical status:** no measured gain, detector confusion, field/symmetry defect, invalid-state feedback, equilibration error, drift, clock stability, or dependence records. No physical margin is established. `NO_NEW_PHYSICAL_INPUT`.
''')
write('RG-F04/STATISTICAL_PROOF.md',r'''# Finite-data rule and proof — RG-F04

Let `k=13` include 12 difference symbols and `INVALID`. Let iid validation outcomes have actual distribution `q`; their empirical histogram is `qhat`. For each proper nonempty subset `S` of outcomes, Hoeffding's one-sided bound is `P(qhat(S)-q(S) >= eps) <= exp(-2 n eps²)`. As `TV(qhat,q)=max_S (qhat(S)-q(S))`, union bounding the `2^k-2` proper nonempty subsets yields:

`P(TV(qhat,q)>eps) <= (2^k-2) exp(-2 n eps²)`.

This holds for iid reset records, *not* arbitrary adjacent readouts of one trajectory. A measurable serial dependence factor cannot be replaced by a guessed scalar ESS. For validation allocation `alpha_val=0.025`, choose `eps=rho/3`, `rho=0.00119114` (hypothetical minimum model-set gap). A sufficient integer n is **40,278,676 genuinely independent full-alphabet readouts**. The count is intentionally conservative and is a **design condition, not achieved sampling or an apparatus feasibility claim**. Replacing rho with any verified smaller residual margin increases n quadratically.

Define two composite model sets `H_FIN` and `H_FLAT` of **all** observed distributions admitted by each hypothesis's own independently calibrated gains, readouts, invalid and drift law. Do not force nuisance settings to be common unless physically justified. Reject `H` only if a conservative **lower bound** on `inf_{p in H} TV(qhat,p)` exceeds eps. Return: FIN_ONLY, FLAT_ONLY, INCONCLUSIVE or NEITHER_MODEL_FITS. For a true `q in H`, on the simultaneous concentration event `TV(qhat,q)<=eps`, H cannot be rejected. Hence false rejection has probability at most `alpha_val` under the declared sampling assumptions. If `inf_{p in H_FIN,r in H_FLAT}TV(p,r)>=rho>2eps`, at most one can remain compatible with the *exact* test; a lower-bound implementation may conservatively return inconclusive.

**Executable conservative example:** `decision_rule.py` represents each model by the gain-indexed exact nominal 12-outcome law augmented with zero `INVALID`, then enlarges it by a radius-`0.005` TV output ball. A 241-point gain grid in `[2.97,3.03]` plus the analytic derivative bound `||dp_g/dg||_TV <= range(Arow)/8` (twice the bound used as safety slack in code) yields a **lower bound on distance to the entire continuous gain set**. This is a software model of *hypothetical* calibration; not a real detector characterization.

The one-event comparator is the nominal maximizer `S={0,2,5,6,7,10}`; no robust positive event gap is admitted under the full hypothetical gain/0.005 detector class from the conservative bound used. An iid Bernoulli Hoeffding test might be cheaper **if** a positive fixed-event minimax gap is independently established; do not promote endpoint samples to an interval theorem.

A 2.5% calibration error allocation is reserved but **no calibration coverage method or data exist**. A 5% overall experimental claim requires separate valid calibration confidence accounting, frozen independently of the validation data. The secondary heat-bath versus Metropolis diagnostic is descriptive here and must not be read as a separately powered hypothesis test.
''')
write('RG-F04/PREREGISTRATION_DRAFT.md',r'''# Prospective preregistration draft — NOT an empirical registration

**Target:** PHYS-014-type N2 pair, q12, two distinctly addressable devices. Only a proposed engineered realization; no evidence of working device exists.

**Primary:** FIN A7 versus fixed trace-matched flat P7 for the **complete stationary 13-outcome difference law**, gain `g3`. `INVALID` retained. Independently calibrate physical realized gains and source compilation across `g in [2.97,3.03]` or **replace** that hypothetical box with measured bounds before freezing a real test. Require certified bound of observed TV error per model ≤0.005 including detector effects allocated to that category. All additional errors reduce rho as stated in RG-F03.

**Preparation:** physically verify equilibration after resetting labels; distinguish g0 uniform calibration from actual g3 steady law. Heat-bath worst-case decay exp(-2t) in model time only; timing in seconds needs independent g0 relaxation and rate-stability checks. Alternative kin­etics require their own burn-in proof.

**Calibration:** disjoint, time-stamped phase-resolved pair preparations; fibre transition rows, true readout confusion/invalids, disturbance and control/compilation/stability, detector drift, independent clock, replicate/reset protocol. Freeze channel IDs, software hashes, data schema, model variant and error budgets before first validation outcome.

**Validation:** iid independently reset full 13-alphabet observations *or* a separately proved dependence-aware procedure. Conditional demonstration: if the calibration justifies rho ≥0.00119114 and validation alpha 0.025, 40,278,676 iid draws suffice for epsilon=rho/3 using the conservative full-alphabet union bound. This is not an achieved or automatically practical target. Alternative power-efficient event test demands its own proven minimax event gap.

**Algorithm:** `decision_rule.py` returns FIN_ONLY, FLAT_ONLY, INCONCLUSIVE or NEITHER_MODEL_FITS from conservative model-set rejection; never force model selection. Never tune gains/readouts on validation. Reserve alpha=0.025 for calibration but define its actual procedure before deployment. A smaller validated rho requires a revised n *before* data opening.

**Secondary diagnostic:** complete joint 13×13 two-time law after independently characterized g0 quench and physical time conversion. Heat bath versus Metropolis is a model comparison, not an empirical kinetic identification. Do not use the already-open predecessor .4/1 results as a blind holdout.

**Stop:** if the phase-resolved interface, calibration, independence or actual uncertainty cannot be paid, stop with NO_NEW_PHYSICAL_INPUT. No acquisition or external experiment is authorized by this draft.
''')

open_inputs='''# Unpaid physical inputs

1. Physical q12 two-device system and a proven mapping from realized device states to labelled a,b, with hardware/controller compilation hashes.
2. Individually addressable absolute-phase preparations and transition-rate/readout calibration sufficient to license 12-state lumpability; otherwise full 144-state instrument.
3. Independently measured gain, field/compilation and drift/error envelopes, with a source independent of the target predictions for each numerical allowance.
4. Complete 13-outcome detector confusion, invalid frequency, outcome-conditioned disturbance/successor state, latency and physical clock rate including rate-stability under g3.
5. Preparation/burn-in evidence at g3, and physical mixing/error estimate.
6. Independent replicate/trajectory IDs and defensible iid or rigorous dependence-aware sampling, calibration/validation partition and record custody.
7. Independently collected validation outcomes not used in compilation, calibration or model selection.
8. An independently justified A7 source atom if asking for a **natural** interaction, distinct from programming the supplied kernel. This does not block a separately approved engineered test.
'''
write('RG-F05/OPEN_INPUTS.md',open_inputs)
write('RG-F05/READINESS_GATE.md',r'''# RG-F05 gate and stopping decision

**Selected status: NO_NEW_PHYSICAL_INPUT.** A prospective, runnable finite-model and statistical design exists conditionally, but *no independently measured numerical calibration or external q12 records* were supplied. It would be incorrect to mark either empirical readiness or completed statistical testing. `DESIGN_READY_PENDING_CALIBRATION` is a useful description of the software, but the requested single primary disposition is NO_NEW_PHYSICAL_INPUT.

Minimum next input: one documented calibration package providing phase-resolved raw pair states/transition trajectories, full 13-outcome readout including invalid, outcome/successor disturbances, actual gain/field/controller stability, physical g0 clock and g3 rate-stability, burn-in, dataset custody, independent replicate definition, and error/concentration evidence. Re-evaluate the complete gap and finite-data rule before authorizing separate measurements. No G campaign automatically follows.

**Claim ladder:** M = supplied finite N2 generator, exact 12-state quotient, Poisson/TV/concentration theorems; PM = independently identified physical dynamical interface (unpaid); PB = frozen model-to-apparatus gain/clock/preparation/readout map (unpaid); PP = prospective full-law prediction conditional on PM/PB; EV = independent held-out real validation (none). Engineered programming of A7 is not source/ontology evidence.
''')
write('RG-F05/NONCONCLUSIONS.md',r'''# Nonconclusions for RG-F01–RG-F05

No physical hardware or dataset was supplied or measured; no physical readiness/validity or measured 0.005 detector error is established. No new independent physical or fundamental source for A7 is identified. The illustrative sample size is *not* data and not a power forecast for dependent trajectories. A numerical zero phase defect in a mathematical model does not certify physical symmetry. An instrument restricted to d and its exact quotient do not validate general absolute-label readouts or phase-breaking fields. Source mathematical inequalities and Poisson tail alone are not a full IEEE/real-device error certificate. There is no M3 universal memory claim, geometric/particle/QM/GR/ToE bridge, global B4 proof, FIN-M0, task335/task337 resolution or ontology update. RG-F06 is architect-only and not executed.
''')

files={
 'RG-F01':('One reusable positive N2 process and 13-output instrument','NUMERICAL_REPLAY_OF_EXACT_MODEL_IDENTITIES',
  'Intertwining L R=R B and complete 1/2/3-time laws supported in the finite, supplied zero-field circulant model',
  ['PROCESS_CONTRACT.json','INSTRUMENT_CONTRACT.json','COMPOSITION_TESTS.json','PROBABILITY_SCOPE.md']),
 'RG-F02':('Symmetry, erasure and phase-dependent instrument negative controls','MATHEMATICAL_OBSTRUCTION_AND_OPEN_PHYSICAL_CALIBRATION',
  'A single h=.01 field violates strong lumpability; phase-dependent emission and uniform erasure show distinct detector failures',
  ['SYMMETRY_AND_ERASURE_CONTROLS.json','CALIBRATION_CAPABILITIES.json','QUOTIENT_VALIDITY_TEST.md','FALLBACK_MODEL_GATE.md']),
 'RG-F03':('Total calibration-to-law uncertainty budget','CONDITIONAL_DESIGN_BOUND_AND_OPEN_INPUT',
  'Inherited rho=.00119114 is conditional on hypothetical gain+detector class; missing field/burn-in/clock/drift/physical constraints are OPEN',
  ['TOTAL_LAW_ERROR_BUDGET.json','CALIBRATION_TO_MODEL_MAP.md','OPEN_BOUND.json','PREDICTION_ENCLOSURES.json','SEPARATION_GATE.md']),
 'RG-F04':('Finite-alphabet rejection/test protocol','PROSPECTIVE_IID_TEST_THEOREM',
  'Under iid full-alphabet sampling, rho=.00119114, alpha_val=.025, eps=rho/3 a sufficient n=40278676; conditional only',
  ['DECISION_RULE.json','DATA_SCHEMA.json','STATISTICAL_PROOF.md','SYNTHETIC_TESTS.json','SYNTHETIC_COMPOSITE_TESTS.json','ONE_EVENT_COMPARISON.json','PREREGISTRATION_DRAFT.md']),
 'RG-F05':('No-input stop and architect handoff','NO_NEW_PHYSICAL_INPUT',
  'No independent real calibration/validation records; design can be submitted for separate review but empirical promotion stops',
  ['OPEN_INPUTS.md','READINESS_GATE.md','NONCONCLUSIONS.md']),
}
for t,(title,status,sentence,artifacts) in files.items():
  write(f'{t}/REPORT.md',f'''# {t}: {title}\n\nExecution: COMPLETED_LOCAL_SCOPE. Scientific disposition: **{status}**.\n\n{sentence}.\n\nLocal outputs: {', '.join(artifacts)}.\n\nThe calculation uses a supplied, finite N2 q12, g3/g0 model. Previous reviewed static gain-cover constants are inherited from the user-provided intake; missing physical records are neither synthesized nor inferred. Numeric tests are not empirical validation.\n''')
  jwrite(f'{t}/RESULTS.json',{'task':t,'execution_status':'COMPLETED_LOCAL_SCOPE',
     'scientific_verdict':status,'arithmetic_status':'exact mathematical identities/proofs; binary64 numerical tests; inherited conditional interval gain cover; no new all-ops outward physical bound',
     'scope':'N2 q12 12-state quotient, full 144-state checking, no hardware',
     'assumptions':['supplied circulant zero-field interaction, fixed gain3/g0','specified heat bath, Metropolis, Barker rate laws','ideal d-based noninvasive readout unless illustrative nuisance declared','iid independent records only for finite-data statistical theorem'],
     'supplied_objects':['rank-seven FIN A7 model','two trace-matched controls','hypothetical gain2.97–3.03','hypothetical .005 TV per-law output class','E06 design disposition'],
     'derived_objects':[sentence],'fitted_objects':[],'empirical_records_available':False,
     'error_budget':'See RG-F03. Remaining physical errors OPEN_INPUT; no achieved detector bound.',
     'discrepancies':'No numerical inconsistencies exceeding declared tolerances; physical assumptions not established.',
     'artifact_names':artifacts,'next_input':'Independent calibration records or architect decision RG-F06; no automatic experiment.'})
  write(f'{t}/NONCONCLUSIONS.md',f'''# {t} — nonconclusions\n\nThis output does not certify physical calibration, natural A7 source, physical generator time, independent empirical validation, or transfer to M3 memory, geometry, particles, QM/GR or ToE. Numerical source checks and conditional mathematical design cannot be relabelled as EV. Missing hardware / calibration data remains OPEN_INPUT.\n''')
  write(f'{t}/REPLAY.md',f'''# {t} — replay\n\nFrom campaign root: `OPENBLAS_NUM_THREADS=1 python run_calibration.py`; then `OPENBLAS_NUM_THREADS=1 python decision_rule.py` and `OPENBLAS_NUM_THREADS=1 python independent_verify.py`. Task-specific artifacts are under `{t}/`. No physical data is fetched or generated. `write_handoff.py` only writes reports; `finalize_package.py` hashes after these outputs.\n''')
  inputs=[ROOT/'source_model.py',ROOT/'source_context/FIN_N2_CALIBRATION_NEXT_CHAIN_USER.md',ROOT/'source_context/FIN_RG_E01_E05_FULL_HANDOFF_20261014.zip',ROOT/'positive_process.py',ROOT/'run_calibration.py',ROOT/'independent_verify.py']
  if t in ('RG-F04','RG-F05'):inputs.append(ROOT/'decision_rule.py')
  input_manifest='\n'.join(f'{sha(p)}  {p.relative_to(ROOT)}' for p in inputs)+'\n'
  write(f'{t}/INPUTS.sha256',input_manifest)

readiness={'status':'NO_NEW_PHYSICAL_INPUT','execution_scope':'F01-F05 complete local equations, simulations, negative controls and conditional statistical design',
 'quotient_mathematics':'admitted under supplied zero-field equivariant input','instrument_physical_validity':'OPEN_INPUT',
 'calibration_reuse':'hypothetical ±1% gain, .005 output TV each; not measured',
 'iid_sample_sufficient_under_hypothetical_budget':decision['finite_alphabet_iid_rule']['n_independent_records_sufficient_Hoeffding_union'],
 'physical_records':0,'independent_source_A7':False,'next_input':'one independently collected, documented N2 physical calibration/custody package described in OPEN_INPUTS.md',
 'architect_only_F06':'NOT_EXECUTED','automatic_experiment_authorized':False}
jwrite('RG-F05/READINESS.json',readiness)

claims=[
 {'id':'F-M01','type':'M','claim':'L R = R B exact zero-field circulant N2 and positive 12-state quotient','status':'EXACT_FINITE_MODEL_THEOREM; numeric checks'},
 {'id':'F-M02','type':'M','claim':'Heat-bath reset semigroup and dimensionless burnin e^-2t','status':'EXACT_FINITE_MODEL_THEOREM'},
 {'id':'F-M03','type':'M','claim':'Phase-dependent first-site field h=.01 breaks exact quotient','status':'NUMERICAL_CONSTRUCTIVE_WITNESS'},
 {'id':'F-M04','type':'M','claim':'Erasure channel has exact zero model information even when g0 histogram uniform','status':'EXACT_NO_GO_EXAMPLE'},
 {'id':'F-M05','type':'M','claim':'Gain-box robust static gap >=.00119114 under two .005 detector allowances','status':'INHERITED_CONDITIONAL_CERTIFICATE_FROM_USER_INTAKE; not repeated independently'},
 {'id':'F-M06','type':'M','claim':'Finite k=13 iid histogram TV concentration, n=40278676 sufficient at rho/3','status':'EXACT_PROBABILITY_THEOREM_WITH_HYPOTHETICAL_RHO'},
 {'id':'F-PM01','type':'PM','claim':'Physical q12 N2 circulant/zero-field controller and rate stability','status':'OPEN_INPUT'},
 {'id':'F-PB01','type':'PB','claim':'Calibrated phase/resolved 13-output instrument and physical clock/burnin','status':'OPEN_INPUT'},
 {'id':'F-PP01','type':'PP','claim':'Full static difference and secondary two-time predictions under calibrated model','status':'DESIGN_ONLY'},
 {'id':'F-EV01','type':'EV','claim':'External independent validation records','status':'NO_NEW_PHYSICAL_INPUT'},
 {'id':'F-SRC01','type':'PM','claim':'Independent physical or fundamental source of A7','status':'NO_NEW_SOURCE'}]
jwrite('CLAIM_REGISTER.json',claims)
write('NONCONCLUSIONS.md',(ROOT/'RG-F05/NONCONCLUSIONS.md').read_text())
write('RG-F06_ARCHITECT_ONLY.md',r'''# RG-F06 reserved — not executed

Architect may choose to: (1) audit independent real calibration records; (2) revise instrument / fall back to 144-state physical parent for a paid defect; (3) submit a complete engineered test for separate human approval; (4) inspect one genuinely independent A7 source atom; or (5) preserve NO_NEW_PHYSICAL_INPUT and stop. Under this handoff **no new physical input was supplied**, so empirical promotion stops. No automatic RG-G is authorized.
''')
write('HANDOFF.md',f'''# FIN full handoff — N2 calibration validity and falsification, RG-F01–RG-F05\n\nDate 2026-10-14. **Status: NO_NEW_PHYSICAL_INPUT / CONDITIONAL_TEST_DESIGN**. RG-F06 is architect-only and was NOT executed.\n\n## Results\n\n- F01: reused and checked exact positive difference quotient, nine combinations of three matrices × three kinetics. Independent 144-state replay checks pass; 13 outcomes incl INVALID compose for arbitrarily long finite sequences under a fixed difference-based nondisturbing emission.\n- F02: exact model symmetry gate, single h=.01 phase-breaking witness {sym['one_declared_0p01_absolute_first_copy_field']['max_block_rate_difference']:.12f} in transition-rate defect; invalid-phase dependence and a uniform-valid erasure detector are explicit counterexamples.\n- F03: inherited **conditional/hypothetical** gain and detector margin **0.00119114**; disjoint added errors must cost less than this before statistical power. No real calibration bounds and no all-operation dynamic interval cert.\n- F04: complete finite-data iid confidence theorem for 13 outcomes. A conservative example with eps=rho/3 and validation alpha=.025 requires **40,278,676 independent records** if that rho is physically validated; not measured or automatically feasible. The continuous-gain composite rejection example distinguishes synthetic ideal FIN/flat and correctly declines erasure data.\n- F05: `NO_NEW_PHYSICAL_INPUT`. The mathematical/statistical workflow is a design, NOT physical validation.\n\n## Evidence\n\nIndependent verifier: **{verify['checks_passed']}/{verify['checks_total']} PASS**; max independent 144-vs-12 joint-law coordinate discrepancy {verify['max_144_vs_12_joint_coordinate_error']:.3e}. Binary64 reproduction is not a source/physics certificate.\n\n## Source provenance\n\nThe original prior RG-E package, as mounted from this conversation, is included under `source_context` unmodified; `source_model.py` is its copied local mathematical reimplementation rather than versioned original `safe_seed.py` bytes. The new uploaded plan is copied byte-for-byte. The RG-E intake supplied by user in conversational text was used as a reviewed claim summary, **not** an independently downloaded file or reverified repo SHA. Repository baseline `08b6d235` is user-reported, not verified against live Git in this run. No predecessor review GAIN_BOX_CERTIFICATE or numerical calibration records are present: its reported sharp bound is inherited explicitly as a conditional premise.\n\n## Replay order\n\n`python run_calibration.py`, `python decision_rule.py`, `python independent_verify.py`, `python write_handoff.py`, `python finalize_package.py`; `python validate_records.py path/to/raw.csv` is syntax-only and does not certify data or iid. See REPLAY.md.\n\n## Disposition\n\nReturn RG-F01–F05 artifacts, full probability arrays, proofs, raw-record schema, independent checks, source manifests and no-input certificate for the architect's RG-F06 review. **Do not create an experiment, spend money, move to a new campaign or promote A7 as a physical source.**\n''')
write('RESEARCH_SUMMARY_PL.md',f'''# Podsumowanie badań FIN — RG-F01–RG-F05\n\n**Werdykt: NO_NEW_PHYSICAL_INPUT — projekt badania gotowy do oceny, bez rzeczywistych pomiarów.**\n\n1. Potwierdzono numerycznie na generatorze 144-stanowym dodatnią 12-stanową redukcję N2, również dla trzech kinetyk i trzech modeli. {verify['checks_passed']}/{verify['checks_total']} testów niezależnego weryfikatora przeszło pomyślnie.\n2. Pole o wielkości 0,01 zależne od bezwzględnej etykiety powoduje defekt redukcji około {sym['one_declared_0p01_absolute_first_copy_field']['max_block_rate_difference']:.8f}. Jednorodny sygnał detektora może całkowicie wymazać ślad oddziaływania, nawet gdy kontrola g=0 wygląda poprawnie.\n3. Z poprzedniego intake przejęto warunkowy margines TV co najmniej 0,00119114. Nie potwierdza on jakości rzeczywistych urządzeń. Dodatkowe błędy kalibracji mogą usunąć rozróżnialność.\n4. Wyprowadzono statystyczną regułę decyzji dla 13 możliwych wyników oraz czterech rezultatów: tylko FIN, tylko P7, nierozstrzygnięte, nie pasuje żaden model. Konserwatywna liczba 40 278 676 jest wymaganą **liczbą niezależnych odczytów w przykładzie warunkowym**, nie wynikiem eksperymentu.\n5. Brak rzeczywistych kalibracji, pomiarów i niezależnego źródła A7. W RG-F05 wydano certyfikat zatrzymania promocji empirycznej. RG-F06 pozostaje wyłącznie dla architekta.\n''')
write('REPLAY.md',r'''# Reproducible replay (local, NO device access)

From root with Python3, numpy, scipy and mpmath installed:

```
OPENBLAS_NUM_THREADS=1 python run_calibration.py
OPENBLAS_NUM_THREADS=1 python decision_rule.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
python write_handoff.py
python finalize_package.py
```

`validate_records.py some_file.csv` audits only syntax/required fields and never treats records as independent unless externally established. `independent_verify.py` uses the 144-state labelled generator and its own branch propagation; it does NOT import positive_process nor run_calibration. The derived reporting and manifest stages are deterministic for a fixed local environment. To compare saved calculations across environments, use tolerances stated in independent verifier rather than byte-identical locks. The source context is immutable.

No formal IEEE/mpmath all-operations proof, remote checkout, physical measurements or experimental preregistration is supplied.
''')

jwrite('ENVIRONMENT.json',{'python':sys.version.split()[0],'platform':platform.platform(),'numpy_scipy_mpmath_needed':True,
 'status':'local numerical only','date':'2026-10-14'})
print('Reports generated; scientific output unchanged.')
