"""Write observational contracts and stage documents; does NOT author mathematical proof.
Metadata is a deterministic rendering of claims present in the separately reviewed proof files.
"""
from __future__ import annotations
import json, hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def put(p,obj):
    (ROOT/p).parent.mkdir(parents=True,exist_ok=True)
    (ROOT/p).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
put('LABELLED_AND_READOUT_CONTRACT.json',{
 'model':'q12 A7, labelled unit-clock LOO HB', 'version':'OBS-C01-v1',
 'gain':'g_N=g_c+lambda/sqrt(N), exact implicit g_c=12/Lambda6',
 'N_scope':'N divisible by 12 as in inherited MT',
 'generator':'sum_a sum_j softmax_j[(g_N/N)A(n-e_{xi_a})]*(f(xi^a->j)-f(xi))',
 'slots':'fixed, predetermined without trajectory selection','parity':'sigma_j=(-1)^j',
 'projection':'count generator from the actual labelled generator',
 'initial':'iid uniform, then deterministic t_N after MT polynomial local mixing',
 'observation_origin':'t_N>=t_loc(N); t_N+sqrt(N) S<= exp(c N) for fixed c<b and fixed S',
 'exchangeability':'exact at each deterministic time; conditional count orbit uniform',
 'orbit_lift_TV':'isometry from count TV to full labelled orbit-mixture TV',
 'path_coupling':'same labelled semigroup after t_N; mismatch<=MT TV; reflect only outside-count edges',
 'physical_source':'NO_NEW_SOURCE','measured_data':'NONE','fitted_objects':[]})
put('HYPERGEOMETRIC_KERNEL.json',{
 'model':'snapshot k fixed tags, population plus count M', 'formula':'comb(M,h)*comb(N-M,k-h)/comb(N,k)',
 'support':'max(0,k-(N-M))<=h<=min(k,M)',
 'block_parity_mean':'m_B=2h/k-1','mean':'m_N=2M/N-1',
 'variance':'(1-m_N^2)*(N-k)/(k*(N-1))',
 'two_disjoint_blocks':'first hypergeom × second hypergeom conditional on h1',
 'cross_covariance':'Cov(m_B1,m_B2|m_N)=-(1-m_N^2)/(N-1)',
 'tie_cases':'h=k/2 or M=N/2 recorded, not silently relabelled',
 'self_report':'exact combinatorics, no independence at finite N'})
put('SNAPSHOT_CHANNEL_CONTRACT.json',{
 'initial_law':'post-MT local equilibrium (not global Gibbs)',
 'collective':'Z_N=N^(1/4)m_N => Z with density proportional exp(chi*lambda*z^2/2-z^4/12)',
 'three_regimes':{'k_small':'k=o(sqrt N), block data asympt independent of Z in TV for bounded tests',
 'k_critical':'k/sqrt N->alpha>0: N^(1/4)m_B=>Z+alpha^(-1/2)*G',
 'k_large':'sqrt N=o(k), k=o(N): N^(1/4)(m_B-m_N)->0 in L2'},
 'conditional_noise':'G~N(0,1) independent of Z in the limit only',
 'two_blocks':'limit Zi=Z+alpha_i^(-1/2)Gi, G1,G2 independent conditional Z',
 'sign_risk_lambda0':'E Phi(-sqrt alpha *abs(Z)); sign ties handled explicitly',
 'entropy_M_I':'OPEN', 'numerical_results_file':'NUMERICAL_RESULTS.json'})
put('TIGHTNESS_TEST_LEDGER.json',{
 'regime':'k/sqrt N->alpha>0, local stationary t_N, critical time s=t/sqrt N',
 'private_contrast':'U_N=sqrt k*(m_B-m_N)',
 'generator':'L U_N=-U_N+O(sqrt k/N)',
 'martingale_bracket':'Gamma(U_N)->2',
 'limit_modeltime':'stationary OU dU=-U dt+sqrt 2 dW, U(0)~N(0,1)',
 'critical_increment_spacing':'h_N=1/sqrt N ->0',
 'critical_increment_variance':'2*(1-exp(-1))/alpha',
 'conclusion':'RAW_PATH_NOT_TIGHT_IN_D_J1',
 'proof_basis':'stationary OU two-time limit + continuity at a.e. deterministic time of any cadlag J1 limit',
 'exceptions':'not a claim for fixed k after integration; static times valid; no physical device'})
put('JOINT_BRACKET_CONTRACT.json',{
 'full_generator':'frozen labelled LOO HB',
 'v':'m_B-m_N', 'drift':'L v=-v+e_N, abs(e_N)<=C/N uniformly',
 'c_a':'1_{a in B}/k-1/N', 'h_a':'2*(1-S_a*q_{xi_a})',
 'q_i':'sum_j sigma_j softmax_j((g_N/N) A(n-e_i))',
 'Gamma_m':'sum_a h_a/N^2',
 'Gamma_v':'sum_a c_a^2 h_a',
 'Gamma_m_v':'sum_a c_a h_a/N',
 'scaled_limits':{'Z':'2s','Y':'2s/k','cross':'0'},
 'jumps':{'Z':'O(N^-3/4)','Y':'O(N^-1/4)'},
 'prior_critical_theorem':'UF actual-generator under MT random local equilibrium initial law',
 'independence':'Brownian drivers independent from diagonal deterministic joint bracket; Z and Y not independent',
 'finite_sample_interval_proof':'NONE','physical_instrument':'NONE'})
put('INTEGRATED_RECORD_CONTRACT.json',{
 'ID':'OBS_FIXED_K_TAGGED_TIME_INTEGRAL_V1','primary':'HB only',
 'k':'fixed finite integer >=1',
 'Y_N':'N^-1/4 integral_{t_N}^{t_N+sqrt N*s} m_B(t) dt',
 'state':'Z_N=N^(1/4)m_N(t_N+sqrt N*s)',
 'limit':'dZ=(chi*lambda*Z-Z^3/3)ds+sqrt(2)dB; dY=Z ds+sqrt(2/k)dW',
 'initial':'Z(0)~rho_lambda stationary 1d; Y(0)=0',
 'cross_noise':'B,W independent innovations',
 'finite_N_error_terms':{'endpoint':'<=4*N^-1/4','LOO_force_integral':'<=C*S*N^-3/4'},
 'evidence_level':'model theorem submission, conditional on inherited MT/UF',
 'kinetic_controls':'MET/BAR open; do not apply HB private corrector'})
put('CONTROL_SCOPE.json',{'HB':'submitted joint path theorem, assumptions stated', 'MET':'OPEN generator Poisson equation and tagged corrector',
 'BAR':'OPEN generator Poisson equation and tagged corrector',
 'label_swap':'outside-primary comparator leaves count process unchanged but changes tag history'})
claims=[
 ('OBS-C01-A','EXACT_IDENTITY','labelled S_N covariance and uniform conditional count-fibre orbit'),
 ('OBS-C01-B','CONDITIONAL_ANALYTIC','MT count TV window lifts to labelled snapshots and paths under same generator'),
 ('OBS-C02-A','EXACT_COMBINATORIAL','hypergeometric tagged snapshot and two-block kernel'),
 ('OBS-C02-B','ANALYTIC_SUBMITTED','k~alpha sqrt N Gaussian channel; small/large k regimes'),
 ('OBS-C02-C','ANALYTIC_SUBMITTED','Bayes sign threshold and bounded risk, symmetric quartic prior'),
 ('OBS-C03-A','ANALYTIC_SUBMITTED','raw k~sqrt N J1 non-tightness via OU modeltime increment'),
 ('OBS-C04-A','EXACT_IDENTITY','L(m_B-m_N)=-contrast+O(1/N), uniform in configuration'),
 ('OBS-C04-B','ANALYTIC_SUBMITTED','fixed-k integrated joint classical SDE with bracket-derived 2/k noise'),
 ('OBS-C05-A','OPEN_PHYSICAL','no actual source/device/units/detector/custody evidence'),
 ('OBS-C05-B','OPEN_THEOREM','MET/BAR fixed-tag integrated martingale corrector not proved')]
put('CLAIM_REGISTER.json',{'campaign':'OBS-C01-C05', 'claims':[{'id':a,'status':b,'claim':c,'scope':'mathematical q12 labelled model, post-MT local sector'} for a,b,c in claims],
 'proof_verification':'independent bounded tests; not mechanized/stochastic-proof execution','fitted_objects':[]})
stages={
 'OBS-C01':('Exact labelled/orbit lift','LABELLED_AND_READOUT_CONTRACT.json','OBS-C01_LABELLED_TRANSFER.md','COMPLETED_EXACT_AND_CONDITIONAL_TRANSFER'),
 'OBS-C02':('Finite hypergeometric and Gaussian snapshot','HYPERGEOMETRIC_KERNEL.json','OBS-C02_SNAPSHOT_CHANNEL_THEOREM.md','COMPLETED_EXACT_AND_ANALYTIC_SUBMITTED'),
 'OBS-C03':('Raw critical path non-tightness','TIGHTNESS_TEST_LEDGER.json','OBS-C03_RAW_PATH_NO_GO.md','COMPLETED_SCOPED_NO_GO_SUBMITTED'),
 'OBS-C04':('Fixed-tag integrated record','JOINT_BRACKET_CONTRACT.json','OBS-C04_INTEGRATED_RECORD_THEOREM.md','COMPLETED_ANALYTIC_SUBMITTED_HB_ONLY'),
 'OBS-C05':('Handoff / architect stop','CLAIM_REGISTER.json','HANDOFF.md','COMPLETE_HANDOFF_STOP')}
inputs=['source_context/OBS_PLAN_USER.md','source_context/MT_INTAKE_USER.md','source_context/theory_core.py','source_context/FIN_MT_FULL_HANDOFF.zip','source_context/MT-003_CONDITIONAL_GIBBS.md','source_context/MT-004_REFLECTED_MIXING_WINDOW.md','source_context/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md']
(ROOT/'SOURCE_INPUTS.sha256').write_text(''.join(f'{sha(p)}  {p}\n' for p in inputs))
for i,(name,title,contract,proof,verdict) in enumerate(((s,*v) for s,v in stages.items()),start=1):
    d=ROOT/name; d.mkdir(exist_ok=True)
    put(f'{name}/RESULTS.json',{'task':name,'title':title,'execution_status':'COMPLETED', 'scientific_verdict':verdict,
      'supplied_objects':['q12 A7 exact-decimal kernel input','labelled LOO HB generator','accepted MT/UF asymptotic theorems','iid uniform first law'],
      'derived_objects':[proof,contract], 'fitted_objects':[],
      'initial_scope':'iid uniform reference; deterministic observation t_N after polynomial MT mixing; post-window law near rho_lambda',
      'time_scope':'t_N+sqrt(N)*S<=exp(cN), fixed S, c<b, N divisible by 12',
      'arithmetic_status':'finite combinatorics exact; generator/risk numeric diagnostics; analytic convergence proofs submitted',
      'discrepancies':['no entire MT-review proof text bundled; inherited as supplied accepted intake'] if i==1 else [],
      'next_action':'INDEPENDENT_ARCHITECT_REVIEW' if i==5 else 'NEXT_TASK',
      'contract_file':contract,'primary_proof_file':proof})
    (d/'REPORT.md').write_text(f'# {name} — {title}\n\nStatus: **{verdict}**.\n\nPrimary evidence: [`{proof}`](../{proof}). Machine-readable contract: [`{contract}`](../{contract}). Main source input is the supplied exact labelled LOO heat bath and the MT/UF admitted local-window theorem. No physical data or fitted objects.\n')
    (d/'PROOF_OR_OPEN_LEMMA.md').write_text(f'# {name} proof register\n\nAnalytic details: [`{proof}`](../{proof}).\n\nProof acceptance is a **submission to independent intake**, not a machine check. Archived input assumptions are explicit.\n\n' + ('For C04, MET/BAR Poisson correctors remain OPEN, but the primary HB lane is paid analytically.\n' if name=='OBS-C04' else 'No unsourced microscopic repair is authorized.\n'))
    (d/'NONCONCLUSIONS.md').write_text(f'# {name} nonconclusions\n\nSee [global NONCONCLUSIONS.md](../NONCONCLUSIONS.md). No observer/device, fundamental physics, measured error bar, information conservation, quantum mechanics or locality follows from this {name} model result.\n')
    (d/'REPLAY.md').write_text(f'# {name} replay\n\nFrom campaign root run `python run_observation.py && python independent_verify.py`. Inspect `{proof}` and `{contract}`. These tests verify formulas/finite arithmetic only.\n')
    (d/'INPUTS.sha256').write_text(''.join(f'{sha(p)}  ../{p}\n' for p in inputs)+f'{sha(proof)}  ../{proof}\n')
print('stage metadata written',len(stages))
