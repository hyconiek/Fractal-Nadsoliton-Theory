#!/usr/bin/env python3
"""Produce transparent static reports and manifest; run after numerical scripts.
Only authorial research reports; not a model proof producer.
"""
import json,hashlib,math,platform,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
J=lambda p: json.loads((BASE/p).read_text())
write=lambda p,s: (BASE/p).write_text(s.strip()+'\n',encoding='utf8')
S=J('RG-E02/SHAPE_PREDICTIONS.json');D=J('RG-E03/SEPARATION_BOUNDS.json')
K=J('RG-E04/MULTITIME_PREDICTION_DESIGN.json')
V=J('INDEPENDENT_VERIFICATION.json')
assert V['all_passed']
p1='FIN_A7_vs_FULL_POTTS_TRACE';p2='FIN_A7_vs_FLAT_P7_TRACE';p3='FULL_POTTS_TRACE_vs_FLAT_P7_TRACE'
fmt=lambda x: f'{x:.12g}'
B=D['N2_pairwise_nominal_tv']; M=D['M3_ordered_pairwise_nominal_tv']

write('RG-E01/SOURCE_INPUT_LEDGER.md', '''# Source and operational input ledger — RG-E01

Typed source classification (model M, physical premise PM, bridge PB, prediction PP, empirical validation EV):

| Object | Actual role | Availability |
|---|---|---|
| 12-label alphabet; M3 three labelled copies; g=5.145228719489142 | M supplied, not discovered | Specified in intake/plan |
| Exact-decimal strict-kernel parameters; projected Fourier support {3,4,5,6} | M supplied, selection open | Versioned safe_seed formula, reproduced locally |
| M3 heat-bath at one attempt per copy; first-two-label observation | M supplied | Simulation law and mathematical observation only |
| Previous RG-D finite laws/moments | M computational; conditional arithmetic | Prior campaign embedded unchanged as ZIP |
| PHYS-014 two q12 nodes, g=3, difference-only readout | PM/PB prospective design | Design-only repository document |
| 24 physical channels; actual hardware selection, weights, temperature and state logs | PM/PB | Missing |
| g=0 physical relaxation measurements, physical time mapping | PB | Missing |
| Readout confusion, invalid/no-hot/multi-hot mass, measurement disturbance/latency | PB | Missing |
| Independent measured mediator K,c or independent selection rule for A7 | Independent source PM | Missing — NO_NEW_SOURCE |
| Physical validation records independent of programming/calibration | EV | Missing |

An engineered realization with *programmed* A7 is a legitimate engineering target but not an independent A7 source. A Gaussian K,c factorization selected after seeing A7 cannot count as source independence. No device, contact, external experiment or purchase is inferred.

Direct public file sources were viewed, but their bytes were not materialized into this working directory: `AGENTS.md`, `PHYS-002/safe_seed.py`, `PHYS-004/{COUNTERMODELS.json,CALIBRATION_CONTRACT.json}`, `PHYS-005/SOURCE_GRAPH.md`, `PHYS-014/REVISED_PREREGISTRATION.md`, and `fin_generative_memory_review_20261014/RG_D06_DECISION.md`. Current git HEAD and their byte hashes were **not** verified. The attached latest instruction and predecessor archive are exact local input bytes and are hashed.
''')
write('RG-E01/CAPABILITY_MATRIX.md', '''# M3 versus PHYS-014 N2 capability matrix

| Capability | M3 memory parent | PHYS-014 N2 design | Current physical evidence |
|---|---|---|---|
| Copies | Three labelled q12 variables | Two q12 nodes | Neither physical run supplied |
| Gain | 5.145228719489142 | 3, plus g=0 | Compilation/design only |
| Allowed static readout | ordered (sigma1,sigma2) | difference (j-i) mod 12 | N2 proposed, not demonstrated |
| Hidden-variable delta prep | requires accessible sigma3 | no sigma3 | Not provided |
| Dynamics | leave-one-out heat bath | hardware law not yet identified | no physical rates |
| Time | attempt-rate dimensionless | physical time requires g0 calibration | missing physical g0 record |
| Invalid detection | idealized absent | explicitly retain | no measured mass |
| Two-time law | verified synthetic prior M3 response | new distinct N2 prospective design | no measured multi-time record |

**No model transfer:** N2 positive test cannot validate M3 hidden-memory coefficients. Retain both contracts under different model identifiers. The only physically proposed readout is N2 difference. Keep E02's M3 conditioning algebra as a hypothetical future requirement, not a supported N2 operation.
''')
write('RG-E01/OPEN_INPUTS.md', '''# RG-E01 OPEN_INPUTS

1. A physical pair of q12 nodes with 24 accepted channels and frozen channel IDs/weights; the existing PHYS-014 proposal is not a record of realized hardware.
2. Independent g=0 clock measurements and evidence that normalized attempt rates do not drift after coupling is applied.
3. Independently measured detector confusion/no-hot/multi-hot, invalid mass, disturbance, timestamp resolution and drift; record all outcomes.
4. Separate calibration and validation trajectories with dependencies/ESS and timestamps.
5. For M3 conditional log-odds specifically: third independent q12 node, conditioning/preparation access and a validated transition/conditional estimator. Not supported by N2.
6. For independent A7 source: one source object/selection rule fixed independently of target A7, rather than a fitted K,c factorization.

No empirical/fundamental claim may bypass these source and bridge gates.
''')
write('RG-E02/IDENTIFIABILITY_THEOREM.md', r'''# Conditional contrast and identifiable quotient — RG-E02

## Theorem 1 (M3 mixed log-odds, model identity)
Let a *supplied* symmetric constant-diagonal q12 matrix A, M=3 labelled copies, and Gibbs exponent (g/6)n^T A n define a finite parent. Given the other two labels b,c, the heat-bath conditional distribution of the first label j is

R_bc(j)=exp[(g/3)(A_jb+A_jc)+h_j]/Z_bc.

For every j,l,b,b',c (same c),

D(j,l;b,b';c)=log[R_bc(j)/R_bc(l)]-log[R_b'c(j)/R_b'c(l)]
             =(g/3)(A_jb-A_lb-A_jb'+A_lb').

**Proof.** Take logarithms of each odds ratio; normalization cancels. On subtracting, the same c contribution and any context-independent supplied one-body field h_j-h_l cancel. The remaining alternating four-entry contrast is as shown. The result uses no attempt clock, but the *availability of conditional probabilities* requires a physical preparation/estimation protocol. End.

Consequently, a complete array of mixed contrasts determines the double-centered matrix (g/3) H A H, with H=I-11^T/q, and not the unconstrained matrix A. The double-centering is a gauge choice, not a source for g or A.

## Exact equivalences / no-go
- A -> A + c 11^T changes all Gibbs logweights by a constant; the normalized law is unchanged, even if the field is fixed to zero.
- For mixed contrasts, A -> A + 1u^T + u1^T leaves all D unchanged. In the full Gibbs model this transform adds g sum_a u_{sigma_a}, and a field h -> h-gu gives exactly the same law. Without a field uncertainty, this second transform need not leave unconditional probabilities unchanged.
- (g,A) -> (g/s,sA) for s>0 leaves the model Gibbs laws and R exactly unchanged; gain and interaction scale are not separately recoverable from normalized conditional probabilities.
- With fixed circulant A, centered log R_00 has Fourier nonzero modes proportional to Fourier modes of A. Ratios among nonzero modes cancel g and the common scale.
- A fully uncalibrated readout can have identical output distributions for every source state, destroying all observational identifiability. That is a *counterexample*, not merely insufficient signal-to-noise.

## Distinct static N2 rule
For two copies N=2 with the same Gibbs normalization, the pair difference d has P(d) proportional to exp[(g/2)A_d] when A is circulant. This rule is a separate N2 model at the PHYS-014 g=3, **not** the stationary law or dynamics of the M3 memory campaign. It has its own invisible constant shift/scale degeneracy. PHYS-014's actual proposed observable is difference-only, not M3 R_bc(j).

## Domain
All algebra is exact given the supplied finite model. Numerical Fourier ratios use the supplied exact-decimal formula in binary64 and are not derived nature-selected constants. If one-body fields vary with conditioning, the cancellation premise fails and requires a new declared nuisance model.
''')
write('RG-E02/GAUGE_REGISTER.json', json.dumps({
 'exact_invisible_constant_matrix_shift':'A_new=A+c*ones_12x12, normalized stationary law unchanged',
 'exact_invisible_gain_interaction_rescaling':'(g,A)->(g/s,s*A) at s>0; fixed labels and field',
 'mixed_contrast_row_column_gauge':'A_new=A+1*u^T+u*1^T; all mixed D unchanged',
 'full_gibbs_field_compensation':'if h_new=h-g*u, the above transform leaves whole finite Gibbs law unchanged',
 'identifiable_matrix_combination':'(g/3)*(I-11^T/12) A (I-11^T/12), under full conditional assay',
 'circulant_shape_observable':'nonzero Fourier ratios of centered log R00 or calibrated N2 log P(d)',
 'unknown_clock':'stationary probabilities cannot determine attempt frequency',
 'unknown_readout':'unrestricted stochastic readout admits complete erasure',
 'status':'mathematical model equivalences; physical identification OPEN_INPUT'},indent=2))
write('RG-E02/OBSERVATION_OBSTRUCTIONS.md', '''# RG-E02 observation restrictions

- The M3 contrast `D` uses *true full conditionals* at externally fixed b,b' and c. With only the ordered Y=(sigma1,sigma2) readout but inaccessible sigma3, the conditional `R_bc` is not directly measurable. An inference of R from the projected stationary 144-vector alone would impose extra model assumptions, not supply this measurement.
- The N2 PHYS-014 difference-only readout **does** admit a 12-label stationary shape test at g=3. It cannot implement the M3 hidden-delta/memory experiment.
- Unknown one-body fields and unknown per-category detector confusion are nuisance parameters, not implicit zeros. In a circulant ideal zero-field N2 model they can be frozen as assumptions; physical enforcement is open.
- Constant A shifts, gain-scale rescaling and field compensated row/column gauges are explicit equivalences, not numerical failures.
- Without any bounded detector channel, a constant-output channel maps all hypotheses to the same observed law. No uniform positive model-separation lower bound survives arbitrary readouts.
''')
write('RG-E03/FALSIFICATION_TABLE.md', f'''# RG-E03 fixed-ideal discrimination table

These numbers refer to the **ideal N2 difference law at g=3**, under trace-matched fixed matrices, not hardware data. The raw TV interval endpoints are independently computed under the mpmath.iv assumption.

| Fixed hypotheses | Ideal TV | Per-hypothesis max common TV error if neither has other error | Lower bound under *hypothetical* gain ±1% + instrument TV≤0.005 |
|---|---:|---:|---:|
| FIN / full Potts | {B[p1]:.12f} | {B[p1]/2:.9f} | {D['N2_nominal_plus_hypothetical_robust_separation'][p1]['lower_bound_if_calibrated']:.9f} |
| FIN / flat P7 | {B[p2]:.12f} | {B[p2]/2:.9f} | {D['N2_nominal_plus_hypothetical_robust_separation'][p2]['lower_bound_if_calibrated']:.9f} |
| Full Potts / flat P7 | {B[p3]:.12f} | {B[p3]/2:.9f} | {D['N2_nominal_plus_hypothetical_robust_separation'][p3]['lower_bound_if_calibrated']:.9f} |

**Readout nuisance rule:** if each allowed observed-law set is within TV radius r_i of ideal p_i, any pair has TV ≥ max(0,TV(p_i,p_j)-r_i-r_j). This follows from the triangle inequality and is a *guaranteed* minimax lower bound under that declared ball constraint, not an exact minimax optimum. Unknown arbitrary detector maps have lower bound zero. The ±1% and 0.005 figures are hypothetical design tolerances: NO measured calibration is present.

**Statistical error:** for 12 iid effective observations, a crude union-Hoeffding guarantee for TV empirical error≤eps at confidence≥1-alpha is n_eff≥(12²/(8 eps²)) log(24/alpha). At eps=0.002 and alpha=0.05 this gives {math.ceil(12**2/(8*.002**2)*math.log(24/.05)):,} effective iid observations, only a conservative sufficient design number, not required/achieved sample size. Correlated trajectories demand separate category-specific mixing/ESS bounds or a predeclared valid block procedure; raw record counts are not iid samples.

An additional trace-preserving mode perturbation PERT_34_10P changes FIN's N2 law by TV≈{D['FIN_vs_trace_perturbation_N2_tv']:.12f}. Arbitrarily smaller perturbations exist, so fixed-alternative separation never proves uniqueness in the unrestricted rank-seven class.
''')
write('RG-E03/NEGATIVE_RESULTS.md', '''# Scoped negative results

1. **NO_SEPARATION (unrestricted readout):** exact erasure-channel counterexample makes all observed laws equal, even with distinct latent distributions.
2. **FIN versus flat P7:** nominal N2 TV≈0.011787, but the conservative ±1% gain + TV≤0.005 readout design allowances do not certify a positive separation. This is lack of *guaranteed* separation under the chosen broad hypothetical ball bound, not a proof that the true nuisance-constrained optimum is zero.
3. **No unique A7 from three models:** the allowed rank-seven perturbation family has arbitrarily near alternatives with TV→0; trace and Fourier-support restrictions leave free ratios.
4. **No natural A7 source:** both engineering compilation and Gaussian target factorization supply the target interaction, not predict it independently.
5. **No M3 conditional assay from PHYS-014:** difference-only N2 readout has fewer capabilities and a different gain, model and experimental aim.
''')
write('RG-E04/ERROR_BUDGET.md', '''# RG-E04 dynamic design error and kinetics scope

Exact model identities for all proposed q12 pair kinetics:
- A two-copy Gibbs equilibrium `pi(a,b)∝exp[(g/2)A_ab]` (constant diagonal).
- A single-copy heat-bath clock: at rate 1 resample a label from `softmax_j[(g/2)A_{j,other}]`.
- Competing proposals to each of the 11 different labels, with Metropolis acceptance `min(1,exp(delta))` and attempt rate `11/12` per copy; Barker acceptance `1/(1+exp(-delta))` and attempt rate `22/12` per copy.
- At g=0 the two adjusted competing generators are exactly the same as heat-bath: per offdiagonal single-site transition rate `1/12`; the nonconstant single-label Fourier correlation is `exp(-t)` in dimensionless time. This is a **defined model calibration**, not evidence that any device retains the same physical clock when interaction is switched on.

The prospective **N2** prediction fixes uniform preparation at g=0, switch to g=3 at model t=0, a noninvasive difference readout at t1=0.4, and another at t2=1.0, with outcomes (0..11, INVALID). The full 13x13 joint is computed from the 144-state parent; no M3 ROM is used. All nominal invalid mass is zero because the *ideal* instrument is defined that way.

**Full physical uncertainty budget is OPEN**, not estimated from simulations:
1. actual finite channel compilation/gain and interaction drift;
2. detector confusion, latency, back-action and invalid/no-hot/multi-hot mass;
3. g=0 calibration time conversion and evidence of rate constancy under g=3;
4. preparation fidelity, nonstationarity, record dependence/ESS and adaptation effects;
5. numerical expm_multiply floating roundoff (independently checked with dense expm; no outward matrix-function certificate);
6. alternative reversible kinetics beyond the three *declared* laws.

No physical PASS or uniform error interval is assigned. The previous M3 heat-bath hidden spectral bound I≤H≤3I is **not** transferred to N2 or alternative kinetic laws. Full multinomial/time-series likelihood analysis requires data with timestamps, invalid records and independently frozen nuisance models.
''')
write('RG-E05/DATA_REQUIREMENTS.md', '''# RG-E05 physical bridge data intake — prospective, unexecuted

**Calibration dataset**, recorded before unblinding validation: 24 or more physical channel IDs, node split, temperature, baseline and coupled rates, all 144 compiled weights, g=0 relaxation trajectories, measured clock/temperature estimates, detector confusion including invalid/no-hot/multi-hot states, readout latency/back-action, physical drift, timestamp precision and replicate/trajectory IDs. Check PHYS-014's predeclared channel quality box and g=0 control.

**Validation dataset**, independent of calibration: consecutive raw timestamped pair-state observations (not only differences), the mapped difference and invalid flag for every attempt, controller/gain command logs, readout quality, metadata for channel dropouts, and mutually identified trajectories/blocks for ESS. The prospective N2 two-time joint requires recording both observations in a way that permits noninvasive filter modelling; disturbance must be independently bounded or a new instrument map explicitly declared before analysis.

**Planning rules:** freeze alternative matrices (trace-match), g and nuisance sets from calibration only; reserve validation records untouched; no post-test retuning; handle multiple candidate comparisons, event selection and serial correlations. Accept physical engineering bridge only under a predeclared TV or likelihood gate with paid total error. Successful programmed A7 does not make A7 independent of the apparatus.

**Source atom:** for nature-selected A7 a distinct measured/predictive source K,c or selection rule must be fixed without consulting the desired A7 fingerprint. None provided here. The physical validation branch requires a human-authorized separate experiment; local computations do not instantiate it.
''')

reports={
'RG-E01':('Operational and source gate','SCOPED_INTERFACE_AND_OPEN_INPUT','M3 has no physical interface; N2 PHYS-014 is a design-only distinct pair protocol. No independent A7 source.',
  ['RG-E01/OPERATIONAL_CONTRACT.json','RG-E01/SOURCE_INPUT_LEDGER.md','RG-E01/CAPABILITY_MATRIX.md','RG-E01/OPEN_INPUTS.md']),
'RG-E02':('Identifiable shape and gauge','EXACT_CONDITIONAL_THEOREM_AND_NOMINAL_SHAPE','Full M3 conditional contrasts recover centered gA, not gain or gauge; physical hidden conditioning is unavailable. N2 difference law is separately predicted.',
  ['RG-E02/IDENTIFIABILITY_THEOREM.md','RG-E02/GAUGE_REGISTER.json','RG-E02/SHAPE_PREDICTIONS.json','RG-E02/OBSERVATION_OBSTRUCTIONS.md']),
'RG-E03':('Robust separation and stop tests','IDEAL_FIXED_ALTERNATIVE_SEPARATION; ACTUAL_CALIBRATION_OPEN','Ideal fixed N2 FIN vs flat P7 TV≈0.011787; unbounded readout gives NO_SEPARATION. Hypothetical detector/gain sets do not guarantee FIN–flat P7 separation.',
  ['RG-E03/SEPARATION_BOUNDS.json','RG-E03/STATIC_INTERVAL_CERTIFICATE.json','RG-E03/DISCRIMINATION_CONTRACT.json','RG-E03/FALSIFICATION_TABLE.md','RG-E03/NEGATIVE_RESULTS.md']),
'RG-E04':('Kinetic discrimination design','VALID_NUMERICAL_FINITE_DESIGN; PHYSICAL_CLOCK_OPEN','Three g0-normalized reversible N2 kinetics share equilibrium yet differ on full 13x13 joint law at g=3; a physical rate/readout interface is missing.',
  ['RG-E04/KINETIC_CONTRACTS.json','RG-E04/MULTITIME_PREDICTION_DESIGN.json','RG-E04/INSTRUMENT_MAPS.json','RG-E04/ERROR_BUDGET.md']),
'RG-E05':('Evidential claim gate','DESIGN_ONLY_NO_NEW_SOURCE_NO_NEW_PHYSICAL_INPUT','N2 preregistration draft and stop certificate exported; no empirical experiment, source selection or M3 transfer.',
  ['RG-E05/PREREGISTRATION.json','RG-E05/STOP_CERTIFICATE.json','RG-E05/PHYSICAL_CLAIM_PATH.json','RG-E05/DATA_REQUIREMENTS.md'])}
for task,(title,verdict,summary,paths) in reports.items():
  write(f'{task}/REPORT.md',f'''# {task} — {title}

**Execution:** COMPLETED_WITH_SCOPED_VERDICT. **Scientific verdict:** `{verdict}`.

{summary}

## Evidence basis

Local exact model algebra plus recomputed finite numerical predictions; N2 ideal stationary separation additionally receives a mpmath.iv conditional interval certificate. Distinguish mathematical `M`, prospective physical premise `PM`, unadmitted physical bridge `PB`, numerical prediction `PP`, and absent experimental evidence `EV`. The underlying q12 kernel, selection support, gain, kinetic laws, and ideal observation maps are **supplied** inputs.

## Files

'''+ '\n'.join(f'- `{p}`' for p in paths)+'''

## Limitations

No independently sourced A7, channel calibration, natural apparatus kinetic law or experimental validation. Nothing here establishes particle/space/field/QM/GR/ToE physics. Full finite model is accepted as *mathematical reference*, not a measured physical system.
''')
  outdict={'task_id':task,'execution_status':'COMPLETED_WITH_SCOPED_VERDICT','scientific_verdict':verdict,
      'model_ids':['N2_q12_g3_difference_design','M3_q12_g5_conditional_mathematical'] if task in ('RG-E01','RG-E02') else ['N2_q12_g3_difference_design'],
      'supplied_objects':['rank-seven FIN Fourier matrix from exact-decimal strict kernel','gain g, observation map, comparison matrices, kinetic options'],
      'derived_objects':[summary], 'fitted_objects':[],
      'arithmetic_status': 'finite numerical + conditional mpmath.iv N2 statics; exact mathematical identities' if task in ('RG-E02','RG-E03') else 'finite numerical/deductive design only',
      'scope': summary,'assumptions':['mathematical Gibbs law and declared generators','real physical error and readout not established'],
      'evidence_availability':{'model_source':'local reconstructed versioned formulas plus embedded prior archive', 'physical_validation_data':False,'independent_A7_source':False},
      'error_budget':{'input':'exact-decimal FIN static interval in RG-E03 only, conditional on mpmath.iv', 'dynamics':'independent float matrix-exponential replay; no rigorous outward full error', 'physical':'OPEN_INPUT'},
      'discrepancies':['none observed beyond floating cross-replay; independently calibrated physical data unavailable'],
      'dependencies':['predecessor RG-D scoped intake and current attached operational plan'],
      'next_action':'Architect RG-E06 decision; do not perform automatically','artifacts':paths}
  (BASE/f'{task}/RESULTS.json').write_text(json.dumps(outdict,indent=2,sort_keys=True)+'\n')
  write(f'{task}/NONCONCLUSIONS.md','''# Nonconclusions

No natural selection of A7, no empirical source or validation, no fundamental clock, geometric dimension, fields/particles, quantum theory, gravity, Standard Model, ToE, or a general theorem of rank-seven uniqueness. A programmed target is not a physical derivation. N2 does not establish M3 memory; conditional input/calibration assumptions are not verified device capabilities.
''')
  write(f'{task}/REPLAY.md','''# Replay

From the campaign root with Python 3, NumPy, SciPy, mpmath available:

```bash
OPENBLAS_NUM_THREADS=1 python run_operational.py
OPENBLAS_NUM_THREADS=1 python interval_static.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
```

`run_operational.py` regenerates prospective model probabilities and diagnostics; `interval_static.py` separately reconstructs the supplied strict-kernel static probabilities in outward intervals (conditional on mpmath.iv); `independent_verify.py` reconstructs FIN with an independent X X^T basis and dense matrix exponentials. The `make_reports.py` script reconstructs report and manifest metadata, so do **not** run it before verifying the distributed source MANIFEST. All reference data are mathematical/synthetic.
''')
  inp=[BASE/'source_context/FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014_USER.md',BASE/'source_context/FIN_RG_D01_D05_FULL_HANDOFF_20261014.zip',BASE/'source_model.py',BASE/'run_operational.py',BASE/'interval_static.py',BASE/'independent_verify.py']
  write(f'{task}/INPUTS.sha256','\n'.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(BASE)}' for f in inp))

write('HANDOFF.md',f'''# FIN RG-E01–RG-E05 — full bounded operational identifiability handoff

**Date:** 14 October 2026. **Architect branch:** `SOURCE_AND_OPERATIONAL_IDENTIFIABILITY_GATE`. **Executed:** RG-E01–RG-E05 only. **RG-E06:** not executed (architect-only).

## Executive verdict

`DESIGN_ONLY / NO_NEW_SOURCE / NO_NEW_PHYSICAL_INPUT`, with mathematical model theorems and independently reproduced finite numerical design contrasts. The full parent is a mathematical reference; no real channel measurement or empirical validation was performed. All `model/physical-premise/bridge/prediction/evidence` roles are recorded.

## 1. Actual input contract

Current user-attached operational plan and the complete predecessor RG-D handoff ZIP are preserved verbatim under source_context, with SHA-256 in the manifest. Current repository files were reviewed via public raw pages (`PHYS-002`, `PHYS-004`, `PHYS-005`, `PHYS-014`, `RG_D06_DECISION`, `AGENTS.md`), but a full byte-verifiable git checkout/HEAD was not available; source code is a declared independent reimplementation, not a reproduced repository executable. Prior user-pasted intake is treated as user-supplied scientific context, not an independently audited original byte file.

## 2. Main new scientific results

| Fixed ideal hypotheses | N2 difference TV at g=3 | M3 ordered TV at old g=5.145... |
|---|---:|---:|
| FIN vs full Potts | {B[p1]:.12f} | {M[p1]:.12f} |
| FIN vs flat P7 | {B[p2]:.12f} | {M[p2]:.12f} |
| Full Potts vs flat P7 | {B[p3]:.12f} | {M[p3]:.12f} |

N2 *static* TV bounds were recomputed from exact-decimal strict-kernel inputs with an independent 65-decimal-digit interval implementation. The 1e-15-scale **public outward endpoints are conditional on the mpmath.iv backend** and a supported fixed finite model, not full IEEE proof. The M3 numbers are numerical comparisons only; the two parents must not be conflated.

Full M3 conditional log-odds contrasts exactly identify only a centered scaled interaction quotient when hidden b/c conditioning exists; this observation is NOT available in the PHYS-014 N2 proposal. Clock and gain-scale freedoms remain. FIN nonzero Fourier ratios relative to mode k3: {', '.join(f'{v:.9f}' for v in S['shape']['FIN_A7']['active_ratios_to_first'])}. Flat P7 ratios are all one and full-Potts has extra positive low modes, under the specified fixed models.

For completely uncalibrated detector readout, a constant-output stochastic channel proves minimax separation zero. Under *hypothetical* ±1% gain and ≤0.005 detector-TV design allowances, the conservative FIN/flat-P7 guarantee is zero. With no gain or detector uncertainty, the ideal FIN/flat P7 difference is {B[p2]:.12f}; half of it ({B[p2]/2:.9f}) is the equal per-model TV-error separation budget. Not a claim that the hardware achieves it.

N2 kinetics with **same supplied equilibrium law** and model-g0-normalized rates yield distinct *two-time* 13×13 joint laws at model t=0.4 and t=1.0, starting from g0 uniform preparation and switching to g=3. At FIN, TV(heat-bath, Metropolis)={K['N2_two_time']['pairwise_kinetic_joint_TV']['FIN_A7']['heat_bath_vs_metropolis']:.12f}, TV(heat-bath, Barker)={K['N2_two_time']['pairwise_kinetic_joint_TV']['FIN_A7']['heat_bath_vs_barker']:.12f}. The ideal invalid outcome is retained explicitly with zero mathematical probability; its **physical** rate is unknown.

## 3. Verification / deliverable scope

- Independent N2 generator built by a different dense implementation and FIN X X^T rather than the producer's circulant shortcut: **{V['checks_passed']}/{V['checks_total']} controls pass**. The reference propagation uses dense `scipy.linalg.expm`; producer uses sparse `expm_multiply`.
- Separate N2 interval producer (mpmath.iv conditional) outputs probability and pairwise TV endpoint enclosures; comparison claims are exact-model, not experiment-level.
- Each task RG-E01–RG-E05 contains `REPORT.md`, `RESULTS.json`, `NONCONCLUSIONS.md`, `REPLAY.md`, `INPUTS.sha256` and its required special artifacts. Script `run_operational.py`, source `source_model.py`, static independent interval `interval_static.py`, `independent_verify.py` and complete output arrays are included. Prior RG-D is embedded unmodified.
- The local original source byte/input hashes are recorded; public source pages viewed via web are **not** represented as locally verified original SHA.

## 4. Final boundary and one executable next step

`RG-E05/PREREGISTRATION.json` is a prospective **N2 engineered** design only, not an approved experiment. It requires independent PHYS-014 actual 24-channel calibration/weights, physical g0 clock, detector/error/invalid logs and disjoint validation before empirical status can be discussed. The separate independent A7 source atom is still missing. Architect RG-E06 must decide whether to retain `PHYSICAL_TEST_DESIGN_ONLY` / `NO_NEW_PHYSICAL_INPUT` or authorize a separately scoped operational measurement; executor has made no architectural decision.

Do not open task335, redo B4/FIN-M0, invent new dimensions/particles, or infer QM/GR/ToE.

See `CLAIM_REGISTER.json`, `NONCONCLUSIONS.md`, `REPLAY.md`, `SOURCE_ACCESS.md`, `RG-E03/FALSIFICATION_TABLE.md`, `RG-E04/ERROR_BUDGET.md`, and `RG-E05/STOP_CERTIFICATE.json`.
''')
write('SOURCE_ACCESS.md','''# Source access / provenance

The instruction in `source_context/FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014_USER.md` is the **actual uploaded input** for this task. The predecessor ZIP is the actually available previous artifact. Both preserve original user-supplied bytes, with hashes in SOURCE_INPUTS.sha256.

Public raw URLs reviewed on 2026-10-14:
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/AGENTS.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-004/COUNTERMODELS.json
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-004/CALIBRATION_CONTRACT.json
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-005/SOURCE_GRAPH.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_03/PHYS-014/REVISED_PREREGISTRATION.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_generative_memory_review_20261014/RG_D06_DECISION.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_generative_memory_review_20261014/HIDDEN_REFRESH_GAP_AND_MOMENTS.md

These pages were read through a web tool, not downloaded into a git checkout. No remote git SHA, repo byte-level manifest or current AGENTS.md byte hash is claimed. We do not modify any upstream code or source package. `source_model.py` is a new independently authored finite formula implementation, and its divergence from source (if any) is a model/input risk to be checked by any later repository intake.
''')
write('NONCONCLUSIONS.md','''# Global nonconclusions

- No physical apparatus was constructed or tested; PHYS-014 was reviewed as a DESIGN_ONLY document, not an experimental result.
- No independently sourced A7 kernel, rank selection, or fundamental space/time/matter/geometry mechanism was derived.
- No M3 memory-law physical validation follows from N2 difference-observation predictions.
- Neither an exact zero minimax with *bounded calibrated* nuisance nor a guaranteed N2 FIN-vs-flat separation under unknown realistic apparatus errors is claimed. Zero worst-case discrimination is proved only when readout is unrestricted; broad hypothetical ball bounds can be inconclusive.
- M3 mixed conditionals assume physically unavailable hidden-variable conditioning; no physical R_bc sample was collected.
- The 13x13 kinetic predictions are finite synthetic probabilities from declared models. Similar equilibrium laws do not validate a physical kinetic clock or each alternate law.
- mpmath.iv static endpoints are conditional on its implementation, not independently proof-assistant verified. Dense/sparse SciPy generator replay is floating numerical, not all-operation outward interval propagation.
- Serial correlation, invalid events, physical drift, state preparation and detector confusion are unknown, not silently set to measured zero.
- Nothing establishes universal selection of A7, spatial dimension, new forces, particles, SM, QM, GR, or a ToE.
''')
write('RG-E06_ARCHITECT_ONLY.md','''# RG-E06 (NOT EXECUTED)

The executor cannot choose the next architecture. Evidence for architect:

1. A finite mathematical prospective N2 PHYS-014 difference and two-time law is fully specified, but apparatus/clock/readout calibration is missing.
2. The M3 mixed-log-odds theorem requires extra hidden conditioning unavailable in PHYS-014; do not transfer it.
3. FIN versus fixed flat P7 ideal N2 difference is nominally TV≈0.011787, but unknown detector and gain nuisance can defeat discrimination; unrestricted-detector no-go is exact.
4. No independent physical A7 source atom is provided; engineered programming is not derivation.

**Eligible architect dispositions under supplied plan:** PHYSICAL_TEST_DESIGN_ONLY or NO_NEW_PHYSICAL_INPUT (depending on decision policy). No choice is made here; architect must freeze one next bounded action and stop criterion.
''')
write('REPLAY.md','''# Minimal replay and verification

1. Verify distributed artifact integrity **before** running any writer: `sha256sum -c MANIFEST.sha256` from this directory (Linux/macOS).
2. Python 3 with NumPy/SciPy/mpmath:

```bash
OPENBLAS_NUM_THREADS=1 python run_operational.py
OPENBLAS_NUM_THREADS=1 python interval_static.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
```

3. Inspect JSON under RG-E02, RG-E03, RG-E04; verify 58/58 checks and conditional N2 static enclosures.
4. Regenerate static narrative/manifest only if desired: `python make_reports.py`. Replays may change floating serialization across OS/numpy/SciPy; compare model laws with numeric tolerances, not byte identity.
5. For independent actual experiment, do not execute without separate authorization and calibration. This pack has no physical readout implementation.
''')
claims=[
 {'id':'E-M01','type':'M','status':'EXACT_THEOREM','claim':'M3 mixed log-odds equal g/3 alternating four-entry interaction contrast, independent of clock, c and fixed fields','evidence':'RG-E02/IDENTIFIABILITY_THEOREM.md'},
 {'id':'E-M02','type':'M','status':'EXACT_THEOREM','claim':'gain-scale, matrix constant-shift and field-compensated row-column gauge equivalences','evidence':'RG-E02/GAUGE_REGISTER.json'},
 {'id':'E-M03','type':'M','status':'CONDITIONAL_INTERVAL','claim':'ideal N2 static pairwise TV enclosures for three fixed supplied matrices','evidence':'RG-E03/STATIC_INTERVAL_CERTIFICATE.json'},
 {'id':'E-M04','type':'M','status':'EXACT_THEOREM','claim':'Unrestricted detector erasure yields no positive minimax separation','evidence':'RG-E03/NEGATIVE_RESULTS.md'},
 {'id':'E-M05','type':'M','status':'NUMERICAL_EVIDENCE','claim':'three g0-normalized N2 reversible kinetics have different ideal full two-time joint laws','evidence':'RG-E04/MULTITIME_PREDICTION_DESIGN.json'},
 {'id':'E-PM01','type':'PM','status':'OPEN_INPUT','claim':'PHYS-014 calibrated 24-channel two-q12 platform actually exists and operates with the proposed law','evidence':'RG-E01/OPEN_INPUTS.md'},
 {'id':'E-PB01','type':'PB','status':'OPEN_INPUT','claim':'independent detector mapping and g0 physical clock make model laws operational','evidence':'RG-E04/ERROR_BUDGET.md'},
 {'id':'E-PP01','type':'PP','status':'DESIGN_ONLY','claim':'prospective 12-label static and 13x13 two-time N2 ideal predictions','evidence':'RG-E05/PREREGISTRATION.json'},
 {'id':'E-EV01','type':'EV','status':'NO_RECORDS','claim':'independent physical validation comparisons absent','evidence':'RG-E05/STOP_CERTIFICATE.json'},
 {'id':'E-SOURCE01','type':'PM','status':'NO_NEW_SOURCE','claim':'No independent mediator K,c or kernel selection rule fixed without A7 target','evidence':'RG-E01/SOURCE_INPUT_LEDGER.md'}]
(BASE/'CLAIM_REGISTER.json').write_text(json.dumps(claims,indent=2)+'\n')
(BASE/'ENVIRONMENT.json').write_text(json.dumps({'python':sys.version.split()[0],'platform':platform.platform(),
 'numpy':__import__('numpy').__version__,'scipy':__import__('scipy').__version__,'mpmath':__import__('mpmath').__version__,
 'verifier_passed':V['checks_passed'],'verifier_total':V['checks_total'],
 'note':'Floating libraries affect byte-level replay; source hashes are only local actual inputs'},indent=2)+'\n')
orig_files=[BASE/'source_context/FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014_USER.md',BASE/'source_context/FIN_RG_D01_D05_FULL_HANDOFF_20261014.zip']
write('SOURCE_INPUTS.sha256','\n'.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(BASE)}' for f in orig_files))
# Avoid hashing the manifest itself or directories. This is a post-execution file snapshot.
allfiles=sorted(p for p in BASE.rglob('*') if p.is_file() and p.name not in ('MANIFEST.sha256',) and '__pycache__' not in p.parts)
write('MANIFEST.sha256','\n'.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(BASE)}' for f in allfiles))
print('REPORTS_WRITTEN',len(allfiles),'FILES_MANIFESTED','INDEPENDENT',V['checks_passed'],'/',V['checks_total'])
