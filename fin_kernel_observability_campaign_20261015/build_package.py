from __future__ import annotations
import json,hashlib,datetime
from pathlib import Path
R=Path(__file__).resolve().parent

def put(n,s):
 p=R/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s.strip()+'\n',encoding='utf-8')
def data(n,x):put(n,json.dumps(x,ensure_ascii=False,indent=2))

a=json.loads((R/'INDEPENDENT_VERIFICATION.json').read_text())
claims={
 'KO-001':{'status':'CONDITIONAL_ANALYTICAL_MATCHED_FAMILY_AND_FIXED_RECORD_NO_GO','model':'A7 vs A7+0.01(P3-P4)','results':['exact rational projector identities and rank','same g_c and leading quartic parity potential','comparator local MT transfer by unchanged gap/symmetry/local-gradient inequalities conditional on admitted MT proof','fixed finite parity-record model-ID Bayes error tends 1/2'],'remaining':['original K02/K03 independent review bytes unavailable','OPEN_CROSS_MODEL full parity model-index experiment','no transferred global phase census']},
 'KO-002':{'status':'CONDITIONAL_JOINT_HALF_BLOCK_SNAPSHOT_LIMIT','model':'q12 LOO HB local conditional Gibbs/actual iid after MT window','results':['exact multivariate hypergeometric block-label law','conditional finite-population CLT with correction (1-alpha)/alpha','joint independent limiting Gaussian fast and sampling noise; independence from quartic soft mode','Var V1=2 null-mode check'],'remaining':['requires admitted comparator MT local window and inherited Gaussian law','no finite-N rate or unbounded moment transfer through TV']},
 'KO-003':{'status':'CONDITIONAL_STRICT_BOUNDED_IDEAL_MODEL_CONTRAST','model':'epsilon 0 vs +0.01; pre-registered event abs(V3)>1','results':['strict p(+.01)>p(0)','asymptotic bounded-event TV gap','equal-prior single-event Bayes model-ID risk'],'remaining':['ideal point models; no calibration or finite-N rate','tiny advantage may disappear under physical nuisance']},
 'KO-004':{'status':'DEFERRED_NOT_REQUIRED','model':'O(1) fast time process proposed only','results':['no dynamic theorem asserted; static bounded contrast suffices for KO question'],'remaining':['full fast-path bracket and tightness not paid, not activated']},
 'KO-005':{'status':'HANDOFF_SCOPED_CONDITIONAL','model':'FIN primary remains A7; comparison is diagnostic','results':['new vs reused register; exactly scoped predictions and nonconclusions','one physical identifiability input contract returned; no automatic new campaign'],'remaining':['physical source, initial preparation and readout calibration absent']}}
for task,v in claims.items():
 base=task
 data(f'{base}/RESULTS.json',{'task':task,'execution_status':'COMPLETED' if task!='KO-004' else 'DEFERRED','scientific_verdict':v['status'],'model_id':v['model'],'new_results':v['results'],'unpaid':v['remaining'],'error_budget':{'finite_N':'no explicit N0 or rate','physical':'no independently sourced instrument','arithmetical':'see INDEPENDENT_VERIFICATION.json','proof':'analytic conditional on predecessor reviews'},'fitted_objects':[],'next_action':'architect intake only'})
 put(f'{base}/REPORT.md', f'# {task} report\n\n**Execution:** '+('not activated; conditional gate deliberately deferred.' if task=='KO-004' else 'completed within the frozen mathematical model.')+'\n\n**Scientific verdict:** `'+v['status']+'`.\n\n## Findings\n\n'+ '\n'.join('- '+t for t in v['results'])+'\n\n## Open/negative scope\n\n'+'\n'.join('- '+t for t in v['remaining'])+'\n\nSee the task-specific proof and parent `HANDOFF.md`; no external experiment or new physical law was run.')
 put(f'{base}/NONCONCLUSIONS.md','# '+task+' — nonconclusions\n\n'+'\n'.join('- '+t for t in v['remaining'])+'\n- No independent physical source for A7, no physical clock/units, locality, conserved charge, measurement apparatus, quantum/gravity/SM/ToE result, strict selector or empirical activation.\n- Model or observation comparison is not a replacement microscopic FIN law.')
 put(f'{base}/REPLAY.md',f'# {task} — replay\n\n`python independent_verify.py` from the package root reproduces the rational Fourier identities, source-parameter interval eigen estimates, Gaussian bounded-event probability intervals and the finite-sample hypergeometric checks. It checks numerical identities, **not** analytic process convergence or the local MT lemma. The analytical arguments are in the corresponding KO-specific Markdown files. No random seed or Monte Carlo run was used.')
 put(f'{base}/PROOF_OR_OPEN_LEMMA.md',{'KO-001':'# KO-001 proof contract\n\nExact rational matching: `STRUCTURAL_HYPOTHESIS_AUDIT.md`, §1. Conditional local theorem transfer: same file §§2–5. The previous MT path proof is an explicitly inherited conditional premise. OPEN_CROSS_MODEL remains open for full microscopic parity model identification.', 'KO-002':'# KO-002 proof contract\n\n`EXACT_FAST_READOUT_KERNEL.md` contains the precise conditional multivariate hypergeometric law, finite covariance, characteristic-function proof sketch of the sampling CLT, stable Gaussian limit and labels-to-count coupling. Joint entire-population fast Gaussian limit uses inherited MT conditional Gibbs. No full process convergence asserted.', 'KO-003':'# KO-003 proof contract\n\n`BOUNDED_CONTRAST_THEOREM.md` derives Gaussian continuity-set probability, strict derivative sign and exact Bernoulli total-variation/Bayes-risk identities. It is conditional on KO-002 weak convergence.', 'KO-004':'# KO-004 open-gate contract\n\nThis test is not required after successful static KO-003. The omitted fast-time path theorem remains a **proposal**. See `DEFERRED.md`. No diffusion process theorem, bracket cancellation or dynamic covariance is admitted.', 'KO-005':'# KO-005 disposition contract\n\nHandoff only, depending on the KO-001–003 scoped proofs and the explicit KO-004 deferral. See parent HANDOFF/CLAIM_REGISTER/NEXT_OPEN_ATOM, with no new physical source or experiment.'}[task])

data('KO-001/MATCHED_MODEL_CONTRACT.json',{'two_model_ids':['A7-primary','A7-plus-0.01-P3-minus-P4-control'],'primary_law':'labelled q12 LOO heat bath; unit attempt rate; denominator N','N_sequence':'N positive multiples of 12','gain':'12/Lambda6','epsilon':['0','0.01'],'kernel_difference':'epsilon(P3-P4)','unchanged':['trace','each diagonal entry','rank7','PSD','centered','D12 symmetry','top simple parity eigenvalue Lambda6','parity line free energy','critical soft drift and leading bracket'],'changed':['Lambda3 + epsilon','Lambda4 - epsilon'],'eigenvalues_nominal':{'Lambda3':1.9614068619764455,'Lambda4':2.19956884933321,'Lambda5':2.298606272079097,'Lambda6':2.3421820411463004},'premises':['post-MT local window for each model','quartic/fast Gibbs CLT as accepted','model-specific local analysis domain'],'independent_cert':'root INDEPENDENT_VERIFICATION.json'})
put('KO-001/LEADING_PARITY_SCOPE.md','# KO-001 — leading parity blindness scope\n\nBoth specified kernels have identical `g_c`, parity-line drift `tanh(m)−m`, its `−m³/3` coefficient, uniform qv 2 and local quartic starting law. On **fixed finite parity-record alphabets**, all limiting cell probabilities coincide if the incoming UF/MT/OBS claims and comparator local hypothesis audit are accepted. Equal-prior binary model discrimination risk then converges to 1/2. The LI theorem concerns predicting future sign under one model; by itself it does **not** prove full-parity information-theoretic blindness to which model was selected. That cross-model stronger property remains open.')
put('KO-001/DEDUPLICATION_LEDGER.md','# KO-001 — deduplication\n\nReused from prior packages: q12 A7 construction and spectra; previous rank7 Fourier static identifiability algebra; MT local Gibbs quartic/fast Gaussian and canonical path; UF critical diffusion, cancellation; OBS integrated one-tag channel; LI within-model bounded sign-forecast risk.\n\nNew in KO: frozen ε=0,+.01 control with independent rational projector verification, comparator-local hypotheses audit, *cross-model* fixed-record asymptotic no-go in that narrowly defined observation, and the separately derived N/2 block observation and event contrast. These are not evidence that nature selects the A7 kernel.')
put('KO-001/INPUTS.sha256','')
put('KO-002/JOINT_SNAPSHOT_LIMIT.md','# KO-002 — joint snapshot theorem and dependency boundaries\n\nLet `X_N=(Z_N,F_1,N,F_3,N,F_4,N)` be the full population statistic at deterministic post-MT time with `X_N⇒(Z,F_1,F_3,F_4)` under admitted local-Gibbs MT. Conditional on the complete count n, the finite-population block error `U_N` satisfies stable convergence to `G∼N_3(0,I)`, because the finite label vectors are bounded, `n/N→u`, and the finite-sampling fraction tends to 1/2. Therefore `(X_N,U_N)⇒(Z,F_1,F_3,F_4,G)` with G independent of X and the F independent among distinct Fourier harmonics and from Z. The full population f_r normalization produces `Var F_r=1/gamma_r`. Hence `V=F+U` has independent Gaussian coordinates with `Var V_r=1+1/gamma_r` and Z independent. The result is an asymptotic **single snapshot**, not a fast-time path.\n\nThe actual iid-start labelled model receives the weak limit by the incoming MT-labelled TV coupling for each model, provided KO-001 comparator local transfer is independently accepted. Nonbounded sample moments do not transfer by TV alone. See `EXACT_FAST_READOUT_KERNEL.md` for hypergeometric support and CLT argument.')
data('KO-002/NORMALIZATION_AND_ERROR_LEDGER.json',{'alpha':'1/2','features':['f1= sqrt(2) cos(2pi j/12)','f3=sqrt(2) cos(6pi j/12)','f4=sqrt(2) cos(8pi j/12)'],'uniform_feature_covariance':'Identity_3','finite_conditional_covariance':'N(N-k)/(k(N-1)) Sigma_n','limit_private_covariance':'Identity_3','fast_full_precision':'12-g_c Lambda_r=12 gamma_r','full_mode_covariance':'Var(F_r)=12/(12 gamma_r)=1/gamma_r','null_r1':'Lambda1(A7)=0, gamma1=1, Var(V1)=2','p_to_u':'from MT local quartic/fast limit','sampling_independence':'asymptotic stable CLT, not finite-N independence','errors':['actual model vs Gibbs local-window TV o(1), no explicit rate','finite conditional hypergeom covariance correction O(1/N) if p fixed','population Sigma_n-I tends to0 in probability, no transferred unbounded moments']})
data('KO-002/SMALL_SUPPORT_CHECKS.json',a['finite_hypergeom_N12_k6'])
put('KO-003/DIMENSIONLESS_FINGERPRINTS.md','# KO-003 — dimensionless information\n\n- `gamma3=1−Lambda3/Lambda6` and `gamma4=1−Lambda4/Lambda6` are model-sensitive dimensionless ratios.\n- The private half-block asymptotic variance equals 1 independent of kernel because the sample fraction is fixed at 1/2.\n- Soft critical parity diffusion and the λ=0 quartic law are identical across the two matched controls (conditional on their local theorem).\n- The event probability `Pr(|V3|>1)` depends on gamma3 and separates the two ideal models.\n- A ratio or identification contrast is not a derivation of A7, its kernel normalization, physical length, field dynamics or physical units.\n- The r1 null mode with Var V1=2 tests arithmetic/readout normalization, not a new fast eigenfrequency.')
data('KO-003/FROZEN_KERNEL_EVENT.json',{'event_id':'KO-E3-ABS-1-ALPHA-HALF-V1','definition':'|V_3,N|>1','freeze_time':'before calculation; supplied KO plan','model_A':'A7 epsilon=0','model_B':'A7 epsilon=+1/100(P3-P4)','N':'12 multiples','gain':'g_c for both','readout':'fixed half-population q12 label snapshot at common deterministic post-MT origin','limit_laws':{'A':a['probability_A7_interval'],'B':a['probability_comparator_interval']},'ideal_binary_event_tv_gap':a['event_probability_gap_interval'],'equal_prior_bayes_risk':a['ideal_equal_prior_risk_interval'],'scope':'ideal asymptotic mathematical sample; no finite-N accuracy/physical calibration'})
put('KO-004/FAST_PATH_THEOREM_OR_GATE.md','# KO-004 gate result\n\n`DEFERRED_NOT_REQUIRED` under the plan activation criterion: strict KO-003 static model contrast already obtained. No OU process independence/bracket proof was performed. Do not promote the proposed dynamic covariance or any physical wave label.')
put('KO-004/BRACKET_AND_INITIAL_PAYMENT.md','# KO-004 deferred payment ledger\n\nNot proved: actual HB block-minus-population generator drift at O(1) time, full predictable bracket for r∈{1,3,4}, mixed fast/block cross brackets, tightness, conditional initial distribution and joint independence of Brownian drivers. The static KO-002 result only fixes one-time covariance and sampling dependence.\n\n**Decision:** DEFERRED_NOT_REQUIRED; no process theorem follows.')
put('KO-005/DEPENDENCY_DAG.md','# KO-005 — dependency graph\n\nSupplied A7 and exact rank7 source → rational perturbation P3−P4 → independent spectra + gap → comparator local MT coercivity/corridor/trapping (conditional on inherited review) → conditional Gibbs soft×fast Gaussian → labelled exchangeability + exact fixed half-block hypergeometric law → conditional sampling CLT → joint Gaussian snapshot → fixed continuity event |V3|>1 → strict ideal model-ID contrast.\n\nSeparate chain: same parity eigenvalue/symmetry/clock → UF critical soft path + OBS fixed-tag channel → common finite parity-record table → equal-prior model-ID risk→1/2 at fixed finite resolution. LI within-model forecast sufficiency does not automatically imply any stronger model-ID no-go.\n\nKO-004 fast-time dynamic process gate was not triggered. No arrow leads from these mathematical tests to physical A7 selection or empirical calibration.')
put('KO-005/REUSED_VS_NEW.md','# KO-005 — classification\n\n| Object | Status |\n|---|---|\n| A7/X/ω,φ,η, q12 & genuine LOO HB | supplied predecessor |\n| parity pitchfork, quartic law, stable transverse gap | reused UF/MT |\n| local conditional Gibbs and iid post-MT law | reused/conditional transferred to comparator |\n| OBS integrated tag record; LI within-model sign risk | reused (not re-proved) |\n| rational P3−P4 exact operator certificate | newly verified in KO |\n| comparison-operator local hypothesis audit | new conditional analytical extension |\n| fixed finite parity-record model-ID no-go | new scope, conditional |\n| half-block conditional hypergeometric and stable CLT | new observation law |\n| strict bounded Gaussian event model contrast | new result, conditional |\n| fast-time dynamic observation process | deferred; not proved |\n| natural/independent A7 source, detector and calibration | absent |')

put('HANDOFF.md', '''# FIN KO-001–KO-005 — kernel-sensitive intrinsic observables

**Date:** 2026-10-15. **Execution:** `fin_kernel_observability_campaign_20261015/`. **Starting decision:** the user-provided LI intake accepts the scoped one-target LI forecast-risk theorem and selects `KERNEL_SENSITIVE_TWO_SCALE_OBSERVABLES`. **Terminal disposition:** `ACCEPT_CONDITIONAL_STATIC_KERNEL_SENSITIVE_IDEAL_CONTRAST; FAST_TIME_PROCESS_DEFERRED`. This is a *submitted model-mathematics handoff*, not a new independent repository intake, formal-machine proof, physical experiment or empirical validation.

## New principal result

For two diagnostic operators `A_0=A7` and `A_{.01}=A7+.01(P3−P4)` with original q12 labelled unit-clock LOO heat bath, fixed `g_c=12/Lambda6`, identical top parity eigenvalue and same local-window assumptions, the critical soft/parity leading laws coincide. Fixed **finite quantized parity** records cannot asymptotically distinguish the two ideal model IDs (equal-prior risk tends to 1/2), but the declared **half-population q12 snapshot** can.

For a predetermined deterministic set of `k=N/2` label slots and `f_3(j)=√2cos(2π3j/12)`, define `V_3,N=√N k^{-1}∑_{a∈B}f_3(ξ_a)` after the post-MT local relaxation. The exact sampling kernel is multivariate hypergeometric. Conditional finite-population CLT + the inherited full fast Gaussian gives

`V_3,N ⇒ N(0,1+1/γ3(ε))`,  `γ3(ε)=1−(Lambda3+ε)/Lambda6`.

The frozen bounded event `E_N={|V_3,N|>1}` has limiting probabilities `p0≈0.708441279664269` and `p+.01≈0.711615925301417`; **strictly** `p+.01>p0`. Its asymptotic Bernoulli TV gap is `Δ≈0.00317464563714729`; ideal equal-prior single-event operator-ID risk is `(1−Δ)/2≈0.498412677181426`. These are interval-computed ideal mathematical numbers, not measurements, real experimental precision, a practical N threshold or calibrated discrimination.

## Dispositions

| Card | Status | Result |
|---|---|---|
| KO-001 | conditional analytic proof/audit + exact projector identities | two trace/diagonal/rank/parity-matched operators; local MT assumptions audited, shared quartic/record law; finite-quantization model-ID no-go |
| KO-002 | conditional theorem with self-contained hypergeometric/sampling proof | independent Gaussian global fast F and sampling U in snapshot limit; `Var(V1)=2` internal null test |
| KO-003 | conditional bounded event theorem | strict model-sensitive event separation, exact formula/risk |
| KO-004 | `DEFERRED_NOT_REQUIRED` | no O(1)-fast-time path theorem needed to establish static contrast; none claimed |
| KO-005 | handoff | dependency DAG, nonconclusions and one next source/calibration input contract |

## What is actually independent here

`independent_verify.py` provides 38 passing checks, including exact `Fraction` arithmetic for the rational P3/P4 projectors, positive-gap interval spectral checks, direct q12 N=12,k=6 exact hypergeometric covariance enumeration and floating auxiliary consistency, and interval-erf probability checks. These are checks of finite algebra/arithmetic, not 38 independent proofs of the asymptotic propositions. The full analytical statements depend on predecessor MT/UF/OBS mathematical results as identified in each task file.

## Sources and custody

The source ZIP packages for UF, MT, OBS, RI and LI and the uploaded KO plan were present in this session and are SHA-256 bound in `SOURCE_INPUTS.sha256` and `PROVENANCE.json`; none is copied into the new output package. The actual versioned source `theory_core.py` was inspected from MT/UF archives; eigenvalues and the kernel construction are independently reconstructed in the new verifier. The user's LI review/intake is a **pasted claim** accepted as an input contract, not a locally verified review file; the original LI and KO review proof-completion/KERNEL_FAMILY_SEED files are not in the mounted ZIPs. The exact seed projectors are instead reconstructed algebraically. Do not misrepresent these inaccessible review bytes as checked.

## Boundaries and next stop

This is local mathematical observability under supplied A7, kinetic law, gain, preparation and ideal half-block readout. No conclusion follows about physical realization, fundamental source of A7, locality/geometry, units/seconds, conserved energy, particles, source/selector transfer, MI/entropy equivalence, quantum/gravity/SM/ToE closure or empirical discrimination. No additional kernel/threshold/epsilon sweeps, fast-time dynamic proof, hardware or messaging were run. `NEXT_OPEN_ATOM.md` returns **one** missing independent input for architect choice; no follow-on execution is authorized.
''')
put('REPORT.md','''# KO campaign scientific report

**Question:** which FIN critical observation laws are shared by matched operators, and can an intrinsic partial observation detect the fast spectral difference?

**Answer:** the soft parity critical law and every previously fixed finite quantized parity record converge to common limits (conditional on admitted local theory), yet the newly frozen N/2 q12 block snapshot has a strictly kernel-dependent Gaussian law. A pre-fixed binary threshold event has a small but nonzero asymptotic model-ID advantage. The operator spectrum is rechecked using exact rational projector identities and interval evaluations. The q12 half-block conditional law is exactly multivariate hypergeometric, and the Gaussian limit is justified by a characteristic-function finite-population argument, not by importing a small-sample collision estimate.

**Conditional dependency:** the new matched comparator uses the same spectral gap and symmetry that underlie inherited MT local Gibbs/trapping and UF/OBS path theorems. This package gives their comparator-specific inequality audit, but is not independent peer review of the entire inherited theorem chain. The LI pasted intake states admission for one-target risk; it is not a new proof of model-index full-parity blindness.

**New mathematical discrimination:** `p0≈0.70844128`, `p+.01≈0.71161593`, `Δ≈0.00317465`, `R*≈0.49841268`. These are ideal asymptotic numbers and cannot be interpreted as an actual A7 physical test.

**Defer KO-004:** static test answers the frozen question. The two-OU process/covariance candidate is not a proved FIN fast-path law in this campaign.

**Verification:** standalone verifier with 38 passing checks; all outputs and hashes are reproducible by REPLAY. No large-N simulation/empirical measurement.

**Next:** retain a single open physical/operational source-and-observation input contract; no automatic new campaign.
''')
put('NONCONCLUSIONS.md','''# KO-001–KO-005: scope and nonconclusions

- No natural mechanism or independent physical evidence selecting A7 over A7+.01(P3−P4).
- Ideal fixed-gain mathematical contrast is not practical identification under nuisance classes, no finite-N N0/rate or sample-size claim.
- Fixed finite parity-record model-ID no-go is **not** full microscopic parity-timestamp experiment equivalence. LI's sign-forecast risk is a different target.
- No broad universality over arbitrary operators, ranks, kinetics, initial laws, times or spatial scales.
- The half-population set is a predetermined subset of exchangeable labels, not a physical spatial region, geometry, entangled subsystem or apparatus.
- Conditional Gibbs is a locally prepared comparison, not the global equilibrium of the primary finite-N law; escaped paths are not postselected away.
- Static snapshot convergence does not prove the O(1) fast-time OU process, driver independence at process level or dynamic covariance formula.
- No physical clock, locality, units, mass, charge, conserved energy, empirical/EV claim, selector/role transfer, FIN-M0/B4/task335, SM/L_total, quantum mechanics, general relativity or ToE closure.
- The numerical interval verification is not a formal theorem prover or peer-review of predecessor MT/UF/OBS/LI proofs.
''')
put('NEXT_OPEN_ATOM.md','''# Next single unresolved atom — KO-SOURCE-OBS-001

**Question for architect, not an automatically activated research campaign:** Can one identify an *independently sourced and specified* real preparation/observation channel whose readout is demonstrably sensitive to the dimensionless fast ratio `gamma3` and whose nuisance-parameter uncertainty is below the theoretical ideal contrast between the two frozen model hypotheses?

**Required new input, currently absent:** (i) independently selected physical candidate/source for an A-like coupling and its supported uncertainty set, before fitting any FIN fingerprint; (ii) a physically defined sampling/readout map corresponding to the mathematical half-block q12 observation or a newly specified competing readout; (iii) physical time/units and independently calibrated noise/preparation uncertainty; (iv) prospective distinct-model outcome with controls and a predeclared acceptance threshold. An inability to supply this input is a *source/identifiability gate*, not a request to increase the number of ideal digits.

**Why next:** KO-003 already proves a nonzero difference under ideal fixed mathematical kernels. Its small Δ≈0.00317465 does not survive unquantified physical nuisance by logical implication; the existing RG-F broad-class no-go is not superseded.

**Do not execute without new independent evidence and human choice.** No new microscopic rule, geometry, dynamical controller, source claim or experiment is authorized.
''')
put('REPLAY.md','''# KO standalone replay

Requirements: Python 3, mpmath and numpy. All inputs are deterministic. No network, fitting, random draws or long finite-N trajectories.

1. From the extracted package root run `python independent_verify.py`.
2. It writes `INDEPENDENT_VERIFICATION.json` and prints a summary of eigenvalues, event probabilities and check count.
3. To validate final package file hashes run `python verify_package.py` (the builder/verify files are excluded from the output manifest only where explicitly stated in the manifest format).
4. Inspect KO-001/STRUCTURAL_HYPOTHESIS_AUDIT.md for comparator transfer premises, KO-002/EXACT_FAST_READOUT_KERNEL.md for conditional combinatorics/CLT, KO-003/BOUNDED_CONTRAST_THEOREM.md for the continuity event, and KO-004/DEFERRED.md for the nonactivated process proof.

Assertion types: rational projector and finite hypergeometric arithmetic (exact); interval transcendental kernel/probability bounds (mpmath.iv at 65 dps); small binary64 internal consistency (numerical tolerances). This does *not* machine-prove the asymptotic local Gibbs, martingale convergence, MT polynomial mixing or the previous LI risk theorem.

`SOURCE_INPUTS.sha256` binds supplied predecessor archives and the attached plan as inputs **outside** this handoff. The output ZIP contains only new KO work.
''')
put('PROPOSED_AGENTS_PATCH.md','''# Proposed AGENTS.md patch — not applied

2026-10-15 KO intake candidate: Two explicitly matched diagnostic rank-7 q12 operators A7 and A7+.01(P3−P4) retain the same parity eigenvalue, trace, diagonal, and conditional critical quartic/one-tag integrated observation laws. At fixed *finite* parity observation resolution, model-ID risk tends to 1/2 under the inherited local theory. A new predetermined N/2-block q12 snapshot has an exact hypergeometric kernel and, conditional on MT local Gibbs and labelled transfer, Gaussian fast/private limit with Var(V3)=1+1/gamma3. Event |V3|>1 is strictly operator-sensitive (ideal limiting probability difference ≈0.003174645637). KO-004 O(1) fast-time process remains deferred/unproved; cross-model full-parity histories remain OPEN. No physical A7 source, units, calibration, locality, conservation or ToE promotion. Apply only after independent mathematical/source-scope intake.
''')
put('RESEARCH_SUMMARY_PL.md','''# KO-001–KO-005 — podsumowanie po polsku

**Główne odkrycie matematyczne:** znaleziono dwie różne macierze oddziaływań o tej samej dominującej niestabilności FIN. Zwykły uśredniony sygnał jednego elementu ma w obu modelach ten sam graniczny rozkład. Natomiast jednorazowy odczyt stanów **połowy populacji** może rozpoznać różnicę w szybszych fluktuacjach.

Wybrano zdarzenie `|V3|>1` *przed* liczeniem. Jego przewidywane prawdopodobieństwo wynosi około 70,8441% dla A7 oraz 71,1616% dla operatora kontrolnego. W idealnym pojedynczym teście rozpoznawania, który z tych dwóch modeli wygenerował zdarzenie, błąd zmniejsza się tylko z 50% do około 49,8413% — różnica jest matematycznie dodatnia, ale bardzo mała.

Weryfikator sprawdził 38 skończonych własności, w tym dokładną geometrię projektorów i poprawkę hipergeometryczną. Dowód rozkładu granicznego korzysta z warunków wcześniejszych twierdzeń MT/UF/OBS oraz argumentu CLT dla próbkowania bez zwracania. Nie ma nowego pomiaru przyrody ani dowodu, skąd bierze się sama macierz A7.

**KO-004 celowo pominięto:** statyczny odczyt już odróżnia operator; nowy proces dynamiczny nie był wymagany. Następna sensowna kwestia wymaga **niezależnego** modelu fizycznego lub zweryfikowanego kanału pomiarowego — nie kolejnych cyfr w tej samej symulacji.
''')

data('CLAIM_REGISTER.json',{'campaign':'KO-001–KO-005','baseline_user_supplied':'380e9dc0 + user-stated LI intake','terminal_scientific_verdict':'CONDITIONAL_STATIC_KERNEL_SENSITIVE_IDEAL_CONTRAST; KO-004 DEFERRED','tasks':claims,'numerical_verification':a['assertions_passed'],'exact_vs_approximate':{'rational_projectors':'exact using fractions','finite_N12 hypergeometric r3/r4':'exact rational product enumeration','spectrum and Gaussian event':'mpmath.iv intervals; conditional on exact source construction','finite other numeric matrices':'binary64 consistency only','large_N_process_limit':'analytic proof submission dependent on accepted MT/UF/OBS theorem inputs'},'not_claimed':['full parity model-ID no-go','physical kernel selection','general observation sufficiency','fast dynamic theorem','finite_N performance threshold','real detector']})

# Source archive provenance. Mounted review files absent; cite user-supplied text as scoped premise only.
inputs={
'KO_PLAN_UPLOADED':'/mnt/data/Wklejony Markdown(20261015-105555).md',
'UF_ARCHIVE':'/mnt/data/files/FIN_UF_001_UF_005_FULL_HANDOFF_20261014.zip',
'MT_ARCHIVE':'/mnt/data/files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip',
'OBS_ARCHIVE':'/mnt/data/files/FIN_OBS_C01_C05_FULL_HANDOFF_20261015.zip',
'RI_ARCHIVE':'/mnt/data/files/FIN_RI_001_RI_005_ONLY_HANDOFF_20261015.zip',
'LI_ARCHIVE':'/mnt/data/FIN_LI_001_LI_005_ONLY_HANDOFF_20261015.zip'}
from hashlib import sha256
hashes={k:sha256(Path(v).read_bytes()).hexdigest() for k,v in inputs.items()}
put('SOURCE_INPUTS.sha256','\n'.join(f'{h}  {k}' for k,h in hashes.items()))
data('PROVENANCE.json',{'source_sha256':hashes,'source_roles':{'KO_PLAN_UPLOADED':'new task plan only, not results','UF_ARCHIVE':'source authorial critical diffusion package','MT_ARCHIVE':'source authorial local Gibbs and canonical path','OBS_ARCHIVE':'source authorial joint observation','RI_ARCHIVE':'source authorial record information','LI_ARCHIVE':'source authorial one-target forecast likelihood'},'not_accessed_as_bytes':['fin_parity_record_sufficiency_review_20261015/INTAKE_20261015.md','fin_parity_record_sufficiency_review_20261015/PROOF_COMPLETIONS.md','fin_parity_record_sufficiency_review_20261015/KERNEL_FAMILY_SEED.md','fin_parity_record_sufficiency_review_20261015/KERNEL_FAMILY_SEED.json','fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md','fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md'],'external_input_assumed':'user pasted LI intake verdict ACCEPT_SCOPED_FULL_PARITY_FORECAST_RISK_LIMIT; version 380e9dc0. This pasted statement is not a byte-verified repository artifact.','repository_changes':[],'predecessors_in_output_zip':False})

for t in claims:
 put(f'{t}/INPUTS.sha256','\n'.join(f'{h}  {k}' for k,h in hashes.items()))

print('Wrote core KO research content')
