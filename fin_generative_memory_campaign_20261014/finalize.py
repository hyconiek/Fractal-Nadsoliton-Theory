"""Emit all research reports/ledgers without changing frozen RG-D04 files.
The report phase is separate from code for generating predictions and opening tests.
"""
import hashlib,json,platform,sys,zipfile,os,time
from pathlib import Path
import numpy as np,scipy,mpmath
from model import *
ROOT=HERE
read=lambda f:json.loads((ROOT/f).read_text())
def md(p,s):(ROOT/p).write_text(s.strip()+'\n')
def jwrite(p,data):write(p,data)

def main():
    checks=model_diagnostics(build(A7()));obs=read('RG-D05/VALIDATION.json');mom=read('RG-D03/MEMORY_REDUCTION.json');inputs=read('RG-D02/INPUT_ENCLOSURES.json');pre=read('RG-D04/PREREGISTRATION.json')
    frozen=read('RG-D04/FROZEN_PREDICTOR.json');old=read('RG-D03/TRAINING_PREDICTIONS.json');verify=read('INDEPENDENT_VERIFICATION.json')
    model={'model_id':'FIN_A7_3copy_leave_one_out_q12_exact_decimal_target','q':12,'M':3,'microstates':1728,'visible_states':144,'zero_control_field':True,'attempt_rate_per_label':1,'gain_exact_decimal':'5.145228719489142','W_kernel':{'omega':'0.18575','phi':'0.1625','eta':'1.8','circle_constant':'mathematical pi'},'coupling':'A7 Fourier modes 3,4,5,6, constant diagonal, supplied not independently sourced','stationary_pi':'normalize exp((g/6) n^T A7 n), n counts of three labelled copies','observed_piY':'sum over hidden sigma3 of SAME three-copy pi','L':'sum of three independent-rate-one leave-one-out categorical heat-bath refresh operators, self-refresh included in uniformization P','uniformization':'L=3(P-I), P stochastic including self-refresh','observable':'Y=(sigma1,sigma2) in order','coarse_equilibrium_generator':'U^T S U mapped to probabilities with sqrt piY; not Gibbs-sampler replacement','U':'U[x,y]=sqrt(pi[x]/piY[y]) when Y(x)=y','S':'diag(sqrt(pi)) L diag(1/sqrt(pi))','memory_projection':'Q=I-UU^T, A=U^T S U, C=QSU, H=-QSQ restricted to Q','g0':'all independent-refresh control; memory coupling zero','diagnostics':checks,'source_roles':{'prior_C_campaign':'embedded unchanged archived ZIP','new_instruction':'source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md','intake':'source_context/INTAKE_20261014.md','repository_AGENTS_and_review':'not included as byte-verified checkout; accepted through attached intake; available local previous archive as separate provenance'},'source_repository_commit':'NO_CURRENT_VERIFIED_GIT_CHECKOUT','notes':'A7 computed from published versioned formula, interval compared to exact-decimal target.'}
    jwrite('RG-D01/MODEL_CONTRACT.json',model)
    pr={'preparation_id_to_exact_parent_p0':{'conditional_00':'pi[(0,0,h)]/piY(0,0) for h=0..11, no fitted hidden value','delta_0':'p0=1 at (0,0,0)','delta_1':'p0=1 at (0,0,1)','delta_2':'p0=1 at (0,0,2), withheld until D05 for source references'},'training':['conditional_00','delta_0','delta_1'],'heldout_preparation':'delta_2','readout_single':'144 probabilities of ordered Y=(sigma1,sigma2)','readout_joint':'first Y00 vs NOT_Y00 (unclassified complement retained) at .35 and all 144 Y at .8; 288 joint outcomes','initial_coordinates':'z0=p0/sqrt(pi), a0=U^T z0, b0=Q z0','no_refit':True}
    jwrite('RG-D01/PREPARATIONS.json',pr)
    md('RG-D01/NOTATION_MAP.md',f'''# Corrected kinetic and memory notation

Old C03 fields `M0_coupled = {0.009952015225620}` and `M1_derivative = {-0.018039981780484}` correspond to the scalar equilibrium quadratic channel for the UNNORMALIZED indicator of visible Y00: **K(0)** and **K'(0)**, respectively, not integrated M0/M1.

Here we use the normalized resolved basis vector `e00` (orthonormal in the reversible coordinates), so `K(0)={mom['FIN_A7']['K0_00']:.15g}` and `K'(0)={mom['FIN_A7']['Kprime0_00']:.15g}`. Multiplying by `piY00={build(A7())['piY'][0]:.15g}` gives the previous unnormalized channel values (up to numerical rounding).

The actual integrals of the scalar memory channel computed by a **numerical CG solve** are `M0=∫ K(t)dt ≈ {mom['FIN_A7']['M0_integral_00_numeric']:.15g}` and `M1=∫ t K(t)dt ≈ {mom['FIN_A7']['M1_first_integral_moment_00_numeric']:.15g}`. These are not interval-certified global integrals; no independent rigorous Q-sector spectral-gap/tail certificate is supplied.

For all claims S is the parent **M=3** reversible generator. The M=4→M=3 static marginal is a different experiment and is NEVER substituted for piY.
''')
    md('RG-D01/REPORT.md',f'''# RG-D01 — one-parent contract

**Execution:** COMPLETED_WITH_SOURCE_PROVENANCE_CONSTRAINT. **Scientific verdict:** SAME_THREE_COPY_PARENT_BOUND.

Frozen q12 M3, gain and input kernel, 1728 microstates, 144 resolved. Projection and source identity are numerically checked. Max detailed-balance symmetry residual of S: {checks['S_symmetry_inf']:.3g}; `U^T S² U − A² = C^T C` max residual {checks['semigroup_defect_identity_inf']:.3g}. All preparations are analytic functions of pi and hidden label. Original AGENTS.md/full remote git checkout unavailable as byte-verified inputs; no SHA attribution to latest commit.
''')
    md('RG-D01/NONCONCLUSIONS.md','# RG-D01 nonconclusions\n\nNo physical A7 source, real clock, new intervention law, verified checkout SHA, spatial geometry or fundamental closed-system ontology. The parent law is supplied mathematics, not laboratory kinetics.\n')
    md('RG-D01/REPLAY.md','# RG-D01 replay\n\n`OPENBLAS_NUM_THREADS=1 python -c "from model import *; print(model_diagnostics(build(A7())))"`\n\nCompare `RG-D01/MODEL_CONTRACT.json` and `RG-D02/INPUT_ENCLOSURES.json`. No source files are modified.\n')

    # RG-D02 numerical / reference accounting
    reference_err={'execution_status':'CONDITIONALLY_BOUNDED_REFERENCE','target_reference_TV':1e-6,'frozen_reference_allowance_TV':1e-8,'basis':'L=3(P-I), P stochastic; any reference finite-time probability is convex Poisson sum of stochastic iterates; tensor/matrix implementations crosschecked independently','max_predeclared_t':1.5,'z_max':4.5,'cutoff':65,'rigorous_exact_formula_poisson_tail_upper_at_z4p5':'2.83153110260728e-52','source_interval_status':inputs['input_status'],'conditional_rate_max_coordinate_error':1e-12,'conditional_initial_probability_max_coordinate_error':1e-12,'P_row_L1_error_bound_from_rate_enclosures':1.2e-11,'source_input_TV_error_bound_at_t_le_1p5':3.3e-11,'float_matvec_allowance_TV':1e-9,'float_poisson_sum_and_libm_allowance_TV':1e-9,'reference_total_TV_allowance':1e-8,'slack_in_reference_allowance':True,'floating_assumptions':['binary64 round-to-nearest and standard small relative error for +,-,*,/ and exponential computations','categorical rates computed by stable softmax are within the declared 1e-12 per-coordinate source interval, as independently checked against mpmath.iv','numpy positive 12-term reductions and 65-term sums accumulate within conservative 1e-9 TV budgets','the mpmath.iv interval engine correctly encloses exact-decimal kernel transendentals'],'unpaid_gate':'no mechanically checked IEEE/libm proof or formal proof assistant; acceptance is conditional on stated arithmetic/runtime assumptions, not an unconditional bit-level theorem','comparator':'sparse P.T@row propagation vs independent tensor-contract product differ <=1e-12 on all heldout cases; checks embedded in RG-D05/VALIDATION.json'}
    jwrite('RG-D02/ERROR_CERTIFICATE.json',reference_err)
    md('RG-D02/REPORT.md',f'''# RG-D02 — finite-window source propagation

**Execution:** EXECUTED. **Verdict:** CONDITIONAL_REFERENCE_ALLOWANCE_1e-8_TV.

Both 1728-state stochastic P and the independent tensor tensor-contraction refresh reproduce equilibrium; self-refresh included in P. The all-outcome 144-state training distributions were evaluated at t=.15,.5,1,1.5 under conditional equilibrium and hidden delta0/delta1. The frozen validation dates/preparation were not opened in the training stage. A 65-term Poisson cutoff at z≤4.5 admits exact mathematical tail ≤2.832e−52, using the geometric-ratio remainder. The source exact-decimal kernel, 144 rate vectors, all 1728 stationary probabilities, and conditional initial distribution were enclosed by mpmath.iv; all source values differ by <1e−12 per scalar from double proposals under that engine.

Nominal all-operation reference allowance: 1e−8 TV (requirement ≤1e−6). The allowance depends on the explicit IEEE/libm/reduction hypotheses in ERROR_CERTIFICATE.json; this is not a formal verified floating-point computation. Double reference tensor-versus-sparse agreement is a crosscheck, not an all-operation proof.
''')
    md('RG-D02/NONCONCLUSIONS.md','# RG-D02 nonconclusions\n\nNo laboratory source validation or formally verified implementation of exp/softmax/IEEE. Poisson tail by itself would not bound input/rounding error. Interval certification is conditional on mpmath.iv.\n')
    md('RG-D02/REPLAY.md','# RG-D02 replay\n\n`OPENBLAS_NUM_THREADS=1 python interval_inputs.py` regenerates conditional outward rate/input enclosures. `OPENBLAS_NUM_THREADS=1 python train.py` regenerates training numerical distributions **only**; this command does not open reserved validations. See `model.uniform_tensor` for independent contraction and `model.uniform_p` for the sparse-P comparison.\n')

    # D03
    md('RG-D03/BASIS_CERTIFICATE.md',r'''# Memory construction and analytic bound

The single-parent microscopic chain is reversible. Each per-copy conditional heat-bath operator is an L²(pi)-orthogonal projection `P_i`, so `L=Σ_i (P_i−I)` and `−3I ⪯ S ⪯ 0`. For resolved `U`, with `Q=I−UUᵀ`, set `A=UᵀSU`, `C=QSU`, `H=−QSQ|Ran Q`. It follows exactly that `UᵀS²U−A²=CᵀC`, and the hidden block satisfies `b'=Ca−Hb`. Elimination gives both `K(t)=Cᵀexp(−tH)C` and **preparation forcing** `Cᵀexp(−tH)b(0)`, not just the stationary kernel.

The frozen predictor uses one deterministic *preparation- and observation-conditioned* Krylov/Galerkin construction: start from `z0=p0/sqrt(pi)` and fixed event projector `D_E` onto states with `Y=(0,0)`; generate `Q S^k z0` and `Q S^ell D_E S^k z0`, sorted by total degree, with twice-reorthogonalized Gram–Schmidt and fixed depth=6, cap=32. No poles, amplitudes, transition rates or probabilities are fitted against trajectories. It yields ≤27 hidden directions on selected FIN preparations, appended to the 144 resolved coordinates. The coefficients `T=VᵀSV`, `w0=Vᵀz0`, `readout=sqrt(piY)*w[:144]` are entirely source-derived. V depends algorithmically on the *declared preparation*, so T is **not** a single universal 171-state generator across all preparations. For each prepared response, T is symmetric negative semidefinite; there is no global CTMC-positivity theorem.

For reference S and compressed T with spectra within [−3,0], the **shifted Taylor comparison** follows from factoring `exp(tS)=exp(−1.5t) exp[t(S+1.5I)]`. For any z,w=Vᵀz and truncation degree d,

`||exp(tS)z−V exp(tT)w||₂ ≤ exp(−1.5t) Σ(k=0..d) t^k/k! ||(S+1.5I)^k z−V(T+1.5I)^k w||₂ + (||z||₂+||w||₂) Pr[Poisson(1.5t)>d]`.

Probability TV after readout is ≤ half the displayed microscopic L² difference by weighted Cauchy–Schwarz. For binary first readout E/complement and full later Y, upper-bound the two branches by one first-stage error, a total filter-projection leakage `||(I−VVᵀ)D_E V w1||₂`, and half the two post-filter propagation errors. All finite sums and tails are saved in `RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json`. These are conditional numerical source-model envelopes with separately declared arithmetic margin, not proof-assistant verified intervals.

At each application of observation D_E the compressed state is explicitly reprojected `w_E=VᵀD_E Vw`; `w_not=w−w_E`. This accounts for measured-event filtering and the complement; it is not multiplication of two independently fitted marginal distributions. Positivity is checked for all 36 exported laws, but positive semidefiniteness of T alone is NOT a general Markov-positivity theorem.
''')
    md('RG-D03/PROCESS_SCOPE.md',f'''# Process scope and preparation

Train: previously open t=.15,.5,1,1.5 and first event Y00 at .15, Y at .5. Tested global 14/17/23/29-dimensional preparation-complete bases: largest training TV at 29 dims for FIN = {old['FIN_A7']['training_by_order']['6']['max_train_TV']:.8g}, but no strong source-based bound for unseen delta2. This negative outcome is retained.

The frozen response candidate is a different **preparation-conditioned channel-aware Krylov rule**, chosen before heldout outputs; depth6 produces 27 hidden basis coordinates for FIN. It uses the same microscopic generator at every preparation; the hidden preparation b0 is NOT fitted and the first readout filter is treated as an operator. The entire 144-outcome t-law and 288-outcome two-time law have been independently checked for probability normalization and positivity in the declared tests. No universal one-model-for-all-preparations Markov process is claimed.
''')
    md('RG-D03/REPORT.md',f'''# RG-D03 — derived memory

**Execution:** EXECUTED. **Verdict:** ACCURATE_PREPARATION_CONDITIONED_FINITE_RESPONSE; GLOBAL_COMPACT_CLOSURE_NOT_ESTABLISHED.

Two production lanes are retained: global preparation-complete 29-mode basis (good training but nontransferable tight bound to delta2), and frozen source-derived channel-aware Krylov rule of 27 modes. Both retain initial forcing, visible-to-hidden coupling, finite-window observation filtering, and g=0 negative control. FIN channel K0={mom['FIN_A7']['K0_00']:.12g}; K'(0)={mom['FIN_A7']['Kprime0_00']:.12g}; integrated M0={mom['FIN_A7']['M0_integral_00_numeric']:.12g}; integrated M1={mom['FIN_A7']['M1_first_integral_moment_00_numeric']:.12g}. M0/M1 here from numerical CG, residual {mom['FIN_A7']['CG_residual_l2']:.3g}, without global interval spectral-tail certificate. Memory also occurs in the two controls.
''')
    md('RG-D03/NONCONCLUSIONS.md','# RG-D03 nonconclusions\n\nDo not equate the preparation-conditioned Krylov predictor with a universal fixed-rate CTMC, unique hidden state reconstruction, fields/particles, physical temporal inertia, geometry, or FIN-specific mechanism. The integrated full memory moments have no rigorous global tail certificate. Positivity checked only on the exported finite test laws.\n')
    md('RG-D03/REPLAY.md','# RG-D03 replay\n\n`OPENBLAS_NUM_THREADS=1 python train.py` for global training; `OPENBLAS_NUM_THREADS=1 python memory_moments.py` for numerical integrated memory; study `bounds.py`, `model.channel_conditioned_basis` and `RG-D03/BASIS_CERTIFICATE.md` for the transfer-proof scope.\n')

    # D04 (avoid changes to frozen objects)
    md('RG-D04/TOTAL_ERROR_BUDGET.md',f'''# Predeclared total finite-law allowance (nonempirical)

Pre-open lock is `RG-D04/LOCK.json`, which binds four files: the frozen predictor, the 144/288-outcome predictions, a-priori per-law Taylor/filter error envelopes, and the preregistration. No withheld parent reference was accessed by `freeze.py`; only predictions, source-operator residuals and training documents were read.

The full one-time and joint-law model errors are bounded by the pre-open *analytic* formulas in `RG-D04/MODEL_ERROR_ENVELOPES_PREOPEN.json`. Worst model term for FIN = {pre['frozen_bound_pretest']['FIN_A7']['max_model_TV_upper_formula']:.12g} TV. The distinct exact-decimal-input + independent reference + numeric-source allowances are recorded in `RG-D02/ERROR_CERTIFICATE.json`. Conservative frozen aggregate margin = **model term + 1e−8 reference/source + 1e−7 additional numerical margin**, which is at most {pre['frozen_bound_pretest']['FIN_A7']['max_total_nominal_TV_upper']:.12g} for FIN. All are less than the fixed 0.01 primary target.

**Arithmetic status:** bound is conditional on interval math correctness and standard binary64 elementary operations. The 1e−7 margin is deliberately explicit rather than fabricated as a rigorous interval error proof for SciPy matrix `expm`. For nonnegative predicted finite outcomes, arithmetic/roundoff is also audited by independent serial rereads; no clipping/renormalization occurs.

**Predicted scope:** all 144 outcomes at t=.35 and .8 for each of four preparations and models, and 288 ordered binary/full joint outcomes at t1=.35, t2=.8, retaining the first-event complement. Model coefficients are rederived by a *frozen algorithm* per preparation, so this is source-based computational response prediction, not a preparation-independent autonomous closure.
''')
    md('RG-D04/REPORT.md',f'''# RG-D04 — freeze and countermodels

**Execution:** PREVALIDATION_FROZEN; **scientific verdict:** ALGORITHM_FREEZE_WITH_CONDITIONAL_NUMERICAL_ENVELOPES. Frozen three model classes (FIN A7, trace matched full Potts, flat rank7) share q/g/kinetics/readout, have no new-g sweep and no fitted kinetic poles. Freeze/lock SHA256 verified by the D05 opener and stdlib verifier. FIN maximal pre-open a priori ROM TV envelope {pre['frozen_bound_pretest']['FIN_A7']['max_model_TV_upper_formula']:.9g}. The reference allowance 1e−8 and further +1e−7 numerical margin are separate.
''')
    md('RG-D04/NONCONCLUSIONS.md','# RG-D04 nonconclusions\n\nNo certified universal kernel, new source atom, identity of physical A7, global Markov positivity or fully mechanically proven FP budget. Controlled synthetic test is not measured evidence.\n')
    md('RG-D04/REPLAY.md','# RG-D04 replay\n\n`OPENBLAS_NUM_THREADS=1 python freeze.py` recomputes frozen prediction files from training inputs and overwrites the original lock, therefore **do not rerun in an audited existing package**. Instead copy the whole research folder to a scratch directory, run the production sequence, compare the newly frozen SHA-256 to the archived immutable originals, then run `verify_frozen.py`.\n')

    # D05
    fin=obs['model_tests']['FIN_A7']
    delta2=fin['delta_2'];vrom=obs['FIN_max_TV_rom'];mark=obs['FIN_max_TV_markov']
    md('RG-D05/REPORT.md',f'''# RG-D05 — reserved finite-window computational validation

**Execution:** TEST_OPENED_AFTER_LOCK_VERIFIED. **Scientific verdict:** 36/36 FINITE-MODEL_COMPUTATIONAL_TESTS_PASS_CONDITIONAL_ARITHMETIC.

Reserved t=.35/.8 and `hidden delta2` reference outputs were constructed **only in verify_frozen.py** after checking all four D04 SHA-256 hashes. Three models × four preparations × (two complete 144-outcome one-time laws and one complete 288-outcome joint law) = 36 independently compared finite laws. FIN largest numerical TV = **{vrom:.12g}**, against largest equivalent stationary-average Markov TV **{mark:.12g}**; overall largest compressed TV across all models {obs['aggregate_max_TV_rom']:.12g}. FIN hidden delta2: TV single at .8 = {delta2['one_time']['0.8']['TV_ROM_vs_parent']:.12g}, and joint TV .35→.8 = {delta2['two_time']['TV_ROM_vs_parent']:.12g}.

All 36 frozen ROM outputs are nonnegative and normalized within floating tolerance; the 288-law first readout includes its complement. The numerical discrepancies are smaller than the **pre-open** conditional finite-spectral error envelopes plus declared rounding allowance; no post-test tuning. Both sparse-state and tensor refresh propagations agree to ≤1e−12 in the scored one-time comparisons. For all scored laws ROM numerically outperforms the autonomous averaged-Markov baseline. This is a synthetic finite-generator computation, not a physical prediction.

**Critical scope:** predictor's 27-mode hidden basis for FIN is derived afresh by the same fixed Krylov rule from the **declared initial preparation and fixed readout**. It is not a single preparation-independent 27-mode memory law. Source/gain/clock are supplied, not independently derived.
''')
    md('RG-D05/NONCONCLUSIONS.md','# RG-D05 nonconclusions\n\nNo empirical validation, new gain holdout, universality theorem, single fixed reduced kinetic operator for every preparation, geometric locality, physical measurement of memory or independent source A7. Tests reuse the same supplied generator and gain; passing a synthetic test cannot establish nature’s laws.\n')
    md('RG-D05/REPLAY.md','# RG-D05 replay\n\n`OPENBLAS_NUM_THREADS=1 python verify_frozen.py` recomputes heldout source references only after verifying the original preregistration hash. `python independent_verify.py` separately verifies the original immutable predictions and scored full laws **without importing the generator**. This command overwrites only post-open validation outputs, not any pre-open file.\n')

    # all task required RESULTS.json metadata + input hashes
    statuses=[('RG-D01','PARENT_MODEL_FROZEN','local-source reconstruction',list(model),['C06 checkpoint and source intake'],'D02/D03'),('RG-D02','CONDITIONAL_REFERENCE_TV_1E-8','mpmath.iv input plus IEEE/libm assumptions',list(reference_err),['D01 model/interval sources'],'D04'),('RG-D03','PREPARATION_CONDITIONED_REDUCTION','numerical finite-L plus conditional Taylor/CG',list(mom),['D01 model and training'],'D04'),('RG-D04','HASHED_PREREGISTRATION_PASS','source-derived Galerkin + conditional bounds',list(pre['frozen_bound_pretest']),['D02/D03 training'],'D05 frozen validation'),('RG-D05','36_OF_36_COMPUTATIONAL_PASS','numerical full laws with predeclared conditional bound',list(obs['model_tests']),['locked D04','D02 independent reference'],'RG-D06 architect only')]
    for key,verdict,arith,derived,deps,nxt in statuses:
        refs=['model.py','bounds.py','RG-D01/MODEL_CONTRACT.json','RG-D02/INPUT_ENCLOSURES.json','RG-D04/LOCK.json','source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md','source_context/INTAKE_20261014.md']
        hashobj={p:sha(ROOT/p) for p in refs}
        result={'execution_status':'EXECUTED' if key!='RG-D04' else 'PREVALIDATION_LOCKED_THEN_VALIDATED','scientific_verdict':verdict,'model_id':'FIN_A7_q12_M3_g5.145228719489142_and_trace_matched_controls','input_hashes':hashobj,'assumptions':['supplied Fourier A7, exact decimal source tuple, M3 leave-one-out heat bath','source/FP arithmetic conditional on stated assumptions','fixed preparation and observation protocol','same preexisting gain for all'], 'derived_objects':derived,'supplied_objects':['A7 formula and gain','three-copy exchangeable mean-field model','observation/readout contract'],'fitted_objects':[],'arithmetic_status':arith,'error_budget':{'reference_TV':1e-8,'extra_arithmetic_TV':1e-7,'model_envelope':'see pretest frozen MODEL_ERROR_ENVELOPES_PREOPEN.json'},'discrepancies':[],'dependencies':deps,'next_action':nxt}
        jwrite(f'{key}/RESULTS.json',result)
        hashed_inputs={p:sha(ROOT/p) for p in refs}
        (ROOT/key/'INPUTS.sha256').write_text(''.join(f'{h}  ../{p}\n' for p,h in sorted(hashed_inputs.items())))

    # root general records
    md('NONCONCLUSIONS.md', '''# Scope and nonconclusions

- No independently sourced physical A7, new-g holdout, laboratory data, material particle/field theory or spacetime derivation.
- No FIN-M0 executable source; no reopened task335, task337, global B4 capacity proof or PHYS-014 apparatus work.
- The earlier RG-C01–C05 archive is embedded without modification; the user-provided intake is independent contextual review and corrects a prior overclaim of C05 dynamic prediction.
- The derived near-exact finite-window projection is not proof of a preparation-independent finite memory closure; it is a **fixed-algorithm, preparation-and-readout-conditioned response reduction**. Adapting deterministic Krylov vectors per new initial preparation remains a substantial modelling/computation requirement.
- mpmath.iv is not formally proved, and the numeric floating-point envelopes rely on declared binary64/exp/scipy conditions, not a theorem-prover certified program.
- Full exact infinite-time memory M0/M1 computed numerically through CG without independent global gap and tail bound.
- No automatic RG-D06 architectural decision. No universality from matching countermodels; if anything, three positive controls weaken A7-specific claims.
''')
    research_overview={'RG-D01':'same three-copy parent and preparation map','RG-D02':'conditional input enclosures and independent uniformization reference','RG-D03':'static-to-kinetic projection, explicit memory forcing, adapted Krylov ROM, memory integrals','RG-D04':'frozen all-probability predictions, bound protocols and SHA lock','RG-D05':'36 heldout finite-window computations, complete 144/288 outcomes, verifiers','RG-D06':'ARCHITECT_RESERVED','source_status':'no verified full new-commit SHA, embedded previous C campaign unchanged','model_conditions':'synthetic, same-g, fixed source'}
    jwrite('CLAIM_REGISTER.json',research_overview)
    md('SOURCE_PROVENANCE.md',f'''# Exact source and provenance ledger

- New research plan: `source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md`, SHA-256 `{sha(ROOT/'source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md')}`.
- User-supplied intake: `source_context/INTAKE_20261014.md`, SHA-256 `{sha(ROOT/'source_context/INTAKE_20261014.md')}`.
- Previous independently produced research archive, retained byte-for-byte: `source_context/FIN_RG_C01_C05_FULL_HANDOFF_20261014.zip`, SHA-256 `{sha(ROOT/'source_context/FIN_RG_C01_C05_FULL_HANDOFF_20261014.zip')}`.
- This environment has no byte-verified working checkout of current `hyconiek/Fractal-Nadsoliton-Theory` nor actual remote HEAD SHA; source intake reported `9e507988` as repository reference, **not established as current HEAD here**.
- Prior archive contains reconstructed safe-seed Fourier construction and model code; the new model was rebuilt independently using the formulas and compared with the previous A7 to about 2e-15 in binary64.
- Exact-decimal source envelope for A7/rates/pi/conditional preparation established conditionally by new interval engine `interval_inputs.py` and saved in `RG-D02/INPUT_ENCLOSURES.json`.
- Only the supplied instruction and intake were considered binding for new work, without replacing any historical repository files or claiming to modify AGENTS.md.
- Frozen predictor/evidence: `RG-D04/LOCK.json` with four SHA256 entries, generated *before* `verify_frozen.py` opened the reserved test.
''')
    md('RG-D06_ARCHITECT_ONLY.md','''# RG-D06 — ARCHITECT RESERVED (NOT EXECUTED)

The computational result supports a strong **preparation-conditioned finite generator response reduction** with explicit observation filtering. A globally fixed low-order memory closure remains unproved. Independently supplied A7 and its physical source remain missing; nothing singles FIN out against controls.

**Questions for the architect:** Is on-demand state/preparation-conditioned Krylov an acceptable compact representation, or is one fixed memory kernel needed? Can a preparation-independent bound be established for a nontrivial preparation set under a ≤32-direction cap? Is the next step a source/measurement contract rather than additional same-parent synthetic tests?

**No automatically selected branch; no proposed fundamental physics promotion.**
''')
    md('HANDOFF.md',f'''# FIN memory campaign — FULL HANDOFF (2026-10-14)

**Execution:** RG-D01 through RG-D05 complete in the declared finite computational scope. RG-D06 **NOT EXECUTED / ARCHITECT ONLY**. The current result is a synthetic finite-generator mathematical/computational test, not physics experimentally validated.

## Scope and correction

Single q=12 **M=3** labelled-copy parent, 1728 microscopic and 144 visible states. Existing g=5.145228719489142, zero control, leave-one-out heat-bath, attempt rate one per copy. Unlike the earlier C05 static consistency result, this is the deferred dynamic finite-window prediction. Source prior-archive/intake are included unchanged and current Git SHA not verified.

## Major results

1. **D01:** Parent stationary marginal, kinetic generator, reversible embedding and initial states frozen. The identity `UᵀS²U − A² = CᵀC` holds to {checks['semigroup_defect_identity_inf']:.3g} numerical max difference.
2. **D02:** Independently tensor-propagated finite laws using stochastic Poisson uniformization; exact-decimal A7/conditional rates/source pi enclosed with mpmath.iv. Poisson truncation tail for m=65,z≤4.5 <2.832e−52. Conditional reference allowance 1e−8 TV (IEEE/libm assumptions disclosed; no formal FP verification).
3. **D03:** Nonlumpable projected process augmented by derived hidden response. FIN scalar normalized-memory `K(0)={mom['FIN_A7']['K0_00']:.12g}`, `K'(0)={mom['FIN_A7']['Kprime0_00']:.12g}`, numerical integrated `M0={mom['FIN_A7']['M0_integral_00_numeric']:.12g}` and `M1={mom['FIN_A7']['M1_first_integral_moment_00_numeric']:.12g}`; all four are different quantities. Global 29-basis transfer had an OPEN general-prep error; selected **fixed algorithm, preparation/observation-conditioned 27-hidden-mode Krylov rule** with no fitted poles.
4. **D04:** Preregistration with complete 144/288-probability predictions, conditional analytic Taylor/filter envelopes and frozen SHA-256 before withheld reference outputs. FIN max prior model TV term `{pre['frozen_bound_pretest']['FIN_A7']['max_model_TV_upper_formula']:.12g}` and total reference/numeric allowance separately.
5. **D05:** 36/36 frozen full-law tests for FIN, Potts, flat rank7 pass ≤0.01 and their predeclared conditional numerical envelopes. FIN max TV `{obs['FIN_max_TV_rom']:.12g}`, FIN max averaged Markov TV `{obs['FIN_max_TV_markov']:.12g}`; across controls max ROM TV `{obs['aggregate_max_TV_rom']:.12g}`. Delta2, `.35/.8`, and exact binary/full time filtering were only opened after lock verification. No retuning.

## Scientific qualification

**Mathematically derived memory dependence and a highly accurate preparation-conditioned Galerkin response** have been demonstrated for a supplied finite model. This is not a single universal 144+27 state CTMC, and coefficients are rederived by the same frozen algorithm from each declared preparation and observation projector. Passing controls demonstrate no unique FIN origin of the mechanism. The model/rate inputs are ungrounded in independent laboratory physics. Source intervals are conditional on mpmath.iv; generic floating arithmetic is under declared allowances, not a mechanically checked proof. Integrated full memory moments are numerically solved, not interval-certified.

## Reproducibility / handoff

Start at `README.md`, read `REPLAY.md`, `CLAIM_REGISTER.json`, `NONCONCLUSIONS.md`, `SOURCE_PROVENANCE.md`. `RG-D04/LOCK.json` and `MANIFEST.sha256` bind all work. Inspect `RG-D05/RAW_MODEL_OUTPUTS.json` and pretest predictions for **complete probability vectors**, not only a reported PASS flag. Run stdlib `python independent_verify.py` then `sha256sum -c MANIFEST.sha256`. For full independent regeneration run on **a copy** of the original research package: `interval_inputs.py`, `train.py`, `memory_moments.py`, `freeze.py`, `verify_frozen.py`, `independent_verify.py`; do not re-freeze in audited original.

## Stop point

Return package to architect for RG-D06 decision. Do not activate hardware spending, B4, FIN-M0, task335, physical quantum/gravity/ToE, or new source propositions merely on the basis of these synthetic computations.
''')
    md('README.md','''# FIN RG-D01–RG-D05 research package (2026-10-14)

Self-contained as executable Python and complete JSON training/test distributions, immutable prior RG-C archive and both newly supplied Markdown attachments. Python 3 with `numpy`, `scipy`, `mpmath` required for scientific regeneration; pure-stdlib `independent_verify.py` checks frozen records and all scored laws without importing model code.

Read in order: `HANDOFF.md`, `RG-D01/MODEL_CONTRACT.json`, `RG-D02/ERROR_CERTIFICATE.json`, `RG-D03/BASIS_CERTIFICATE.md`, `RG-D04/PREREGISTRATION.json` & `LOCK.json`, `RG-D05/VALIDATION.json`, and `NONCONCLUSIONS.md`.

The previous campaign's review corrects earlier C05 overpromotion: first static event was not a dynamic holdout; this new campaign explicitly does a finite-time dynamic heldout test. Not a Git commit or hardware experiment.
''')
    md('REPLAY.md','''# Reproduction steps

Extract ZIP, `cd fin_generative_memory_campaign_20261014`.

1. Verify saved files and freeze without altering them: `python independent_verify.py` and `sha256sum -c MANIFEST.sha256`.
2. For a full regeneration, make a SECOND working-directory copy so you do not overwrite frozen preregistration. In that copy run (in this strict order):

```
export OPENBLAS_NUM_THREADS=1
python interval_inputs.py
python train.py
python memory_moments.py
python freeze.py           # freeze must PRECEDE any validation reference
python verify_frozen.py    # independently opens reserved t=.35/.8, delta_2
python independent_verify.py
```

The data-producing scripts are deterministic for fixed numpy/scipy/BLAS. Hashes of float results are environment-sensitive; scientific equivalence means numerical/envelope agreement, NOT automatically byte-for-byte frozen SHA across different versions. The original saved LOCK+MANIFEST bind the delivered artifact.

For exact recheck of Poisson-reference correctness, `model.uniform_tensor` is distinct from sparse `model.uniform_p`. A minimum floating/reference accounting disclosure is at `RG-D02/ERROR_CERTIFICATE.json`. Independent verifier reads raw frozen predictions and reference vectors and recomputes scores without model imports.
''')
    jwrite('ENVIRONMENT.json',{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mpmath.__version__,'OPENBLAS_NUM_THREADS_for_commands':1,'date':'2026-10-14','operating_system':platform.platform()})
    # Root README contains model discrepancies as admission
    jwrite('PROVENANCE.json',{'source_plan_sha256':sha(ROOT/'source_context/FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md'),'source_intake_sha256':sha(ROOT/'source_context/INTAKE_20261014.md'),'previous_archive_sha256':sha(ROOT/'source_context/FIN_RG_C01_C05_FULL_HANDOFF_20261014.zip'),'frozen_sha256':read('RG-D04/LOCK.json')['pretest_files_sha256'],'model_not_remote_sha_verified':True})
    md('RESEARCH_SUMMARY_PL.md',f'''# Wyniki badania pamięci FIN – 14 października 2026

Badanie dotyczy **modelu skończonego** (3 kopie, 12 etykiet, 1728 stanów), a nie empirycznego potwierdzenia fizyki FIN.

- D01: zgodne prawo stacjonarne i generator tego samego rodzica.
- D02: referencyjna propagacja Poissona, dane wejściowe w przedziałach `mpmath.iv`, modelowany budżet referencji 1e−8 TV (warunkowy numerycznie).
- D03: jawna pamięć `K(t)` i wymuszenie od nieekwilibrialnego przygotowania `b(0)`. Obliczono oddzielnie K(0), K'(0) i numeryczne momenty całkowe M0/M1. Uniwersalne zamknięcie ze stałą bazą NIE zostało udowodnione.
- D04: przed otwarciem odłożonych wyników zapisano algorytm, predykcje 144/288 elementów, oszacowania i ich SHA-256.
- D05: **{obs['n_tests_meeting_bound_and_0.01']}/{obs['n_test_laws']}** testów dla FIN oraz dwóch kontrmodeli spełniło TV ≤0.01 i zadeklarowane oszacowania. FIN max TV modelu z pamięcią **{obs['FIN_max_TV_rom']:.10g}**, a analogicznego modelu Markowa **{obs['FIN_max_TV_markov']:.10g}**.

Ważne: dobra redukcja jest **warunkowana przygotowaniem i kanałem odczytu**; dla nowej preparacji baza Kryłowa jest generowana ponownie według niezmiennej reguły bez dopasowywania parametrów. Nie jest to odkrycie jednego niezależnego od stanu procesu Markowa ani dowód unikalności A7. Zgodnie z instrukcją RG-D06 zostaje zarezerwowany do decyzji architektonicznej.
''')
    # Verify locked objects have not moved
    lock=read('RG-D04/LOCK.json')
    assert all(sha(ROOT/p)==h for p,h in lock['pretest_files_sha256'].items()),'PREOPEN LOCK WAS MODIFIED!'
    # Manifest generated only AFTER docs, excluding manifest itself, archive and caches.
    files=sorted(p for p in ROOT.rglob('*') if p.is_file() and not ('__pycache__' in p.parts or p.name=='MANIFEST.sha256'))
    (ROOT/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(ROOT).as_posix()}\n' for p in files))
    print('COMPLETE',len(files),'manifest records; D05 numerical max',obs['aggregate_max_TV_rom'],'verified',verify['overall_pass'])
if __name__=='__main__':main()
