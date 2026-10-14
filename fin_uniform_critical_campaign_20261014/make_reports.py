#!/usr/bin/env python3
"""Deterministically creates authorial task documentation, not a proof verifier."""
import json,hashlib,shutil
from pathlib import Path
p=Path(__file__).resolve().parent
plan=Path('/mnt/data/Wklejony Markdown(8).md')
intake=Path('/mnt/data/Wklejony Markdown (2)(3).md')
# copy exact uploaded inputs if available on original execution, retain those provided in package when replaying
for src,name in [(plan,'PLAN_UF_USER.md'),(intake,'KI_INTAKE_USER.md')]:
 if src.exists() and not (p/name).exists():shutil.copy2(src,p/name)
D=json.loads((p/'DIAGNOSTICS.json').read_text())
V=json.loads((p/'INDEPENDENT_VERIFICATION.json').read_text())
modules=[
('UF-001','Exact symmetric-start spectral/initial contract','UF-001_UNIFORM_CONTRACT_AND_INITIAL_PROOF.md','CONTRACT_PASS_PREMISE_CONDITIONAL','Accepted unique top k6 eigenvalue and exact source input; iid-uniform preparation and fourth-moment proof',{'critical_gain_nominal':D['critical_gain_proposal'],'fast_gap_nominal':D['fast_gap_proposal']}),
('UF-002','Microscopic parity characteristics and selection rule','UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md','DERIVED_CHARACTERISTICS_AND_LIPSCHITZ_GATE_PASS','Symmetry and logarithmic-mean mobility pay the MET kink; no exact parity lumpability',{'uniform_bracket_coefficients':{k:v['uniform_m_bracket_coeff'] for k,v in D['rules'].items()},'critical_cubic_coefficient_exact':'-1/3'}),
('UF-003','Unbiased-start quartic critical process','UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md','SCOPED_THEOREM_SUBMITTED_FOR_INDEPENDENT_INTAKE','Local energy, tube and martingale proof with accepted spectrum as an inherited premise',{'scaled_fast_occupation_order':'O(N^-1)','fast_tube_exit_order':'O(N^-1/10)','limit_generator':'(chi*lambda*z-z^3/3)*d/dz+d2/dz2'}),
('UF-004','Local quartic and full-Gibbs limit orders','UF-004_QUARTIC_VS_FULL_GIBBS.md','LOCAL_GLOBAL_MISMATCH_PROVED_WITH_ACCEPTED_COEXISTENCE','Limiting diffusion stationary density is not global finite-N Gibbs stationary distribution',{'F_delta0_nominal':D['at_gc_full_gibbs_competitor']['F_delta0'],'lambda0_quartic_variance_nominal':V['quartic_variance_lambda0']}),
('UF-005','Final theorem/evidence handoff','HANDOFF.md','HANDOFF_READY_FOR_ARCHITECT_AND_INTAKE','No automatic next branch; separate source/input role from analytical consequences',{'independent_arithmetic_checks_passed':V['passed'],'independent_arithmetic_checks_total':V['total']})]
proof_global={
  'execution_status':'COMPLETED',
  'proof_status':'ANALYTIC_PROOF_SUBMITTED_NOT_EXTERNAL_PEER_REVIEW',
  'supplied_objects':['q12','A7 exact-decimal construction','LOO heat-bath generator','diagnostic MET/BAR declared rates','iid-uniform initial ensemble','critical gain sequence gN'],
  'fitted_objects':[],
  'arithmetic_status':'exact analytic identities with inherited interval spectral premise; finite scalar diagnostics binary64 only',
  'physical_evidence':'NONE',
  'expected_errors':'no certified finite-N constants or TV rate; asymptotic O estimates with analytic existence constants',
  'limits':['finite scaled S','initial iid uniform','bounded lambda','fixed algebraic source matrix']}
for code,title,proof,status,scoping,values in modules:
 q=p/code;q.mkdir(exist_ok=True)
 rec=dict(proof_global)
 rec.update(task=code,task_title=title,scientific_verdict=status,scope=scoping,derived_objects=values,dependencies=['accepted unique-top-mode spectrum','existing actual LOO generator'],next_action='Independent intake, not automatic experimental or new campaign')
 (q/'RESULTS.json').write_text(json.dumps(rec,indent=2,ensure_ascii=False)+'\n')
 (q/'REPORT.md').write_text(f'# {code} — {title}\n\n**Execution:** COMPLETED. **Scientific verdict:** `{status}`.\n\n{scoping}.\n\n**Main scientific proof:** [`{proof}`](../{proof}). Source conditions, error orders and nonconclusions remain in the shared proof texts and claim register. This report is an authored status summary, not independent proof validation.\n')
 (q/'PROOF_OR_OPEN_LEMMA.md').write_text(f'# {code} proof pointer\n\nThe complete analytical derivation appears in [{proof}](../{proof}). Check the inherited spectral/source input independently before promoting the theorem to an accepted reviewed result.\n\n**Status:** {status}. Numerical tests corroborate algebra but do not prove a process limit.\n')
 (q/'NONCONCLUSIONS.md').write_text(f'# {code} limits\n\n{scoping}. No external physical measurement, kernel source, fundamental units or empirical closure. This task does not establish global equilibrium from the local diffusion, finite-N uniform error rates, or unique primordial sign selection.\n')
 (q/'REPLAY.md').write_text('Run `python uf_analysis.py` and `python independent_verify.py` from the parent directory (requires numpy/scipy). Read the parent proof document. Compare input hashes before interpreting numerical outputs. The scripts are model diagnostics and consistency checks, not an automated theorem prover.\n')
 (q/'INPUTS.sha256').write_text(''.join(hashlib.sha256((p/n).read_bytes()).hexdigest()+f'  ../{n}\n' for n in sorted(['source_context/theory_core.py','PLAN_UF_USER.md','KI_INTAKE_USER.md','DIAGNOSTICS.json','INDEPENDENT_VERIFICATION.json',proof])))
claims=[
 {'id':'UF-01','class':'M','status':'EXACT_IDENTITY_WITH_ACCEPTED_SPECTRUM','claim':'At g_c top parity k6 is the single null tangent mode and the other ten occupation directions are stable.'},
 {'id':'UF-02','class':'M','status':'ANALYTIC_PROOF_SUBMITTED','claim':'iid-uniform preparation gives EZ0^2=O(N^-1/2) in the rescaled coordinate and fast fourth moment O(N^-2)'},
 {'id':'UF-03','class':'M','status':'ANALYTIC_PROOF_SUBMITTED','claim':'All three frozen kinetic rules have cubic coefficient -1/3 and parity bracket 2 at uniform criticality; MET uses Lipschitz mobility, not C3 rates.'},
 {'id':'UF-04','class':'M','status':'ANALYTIC_PROOF_SUBMITTED','claim':'Fast energy LV<=-kV+C/N on an interior domain, scaled occupation O(N^-1), fast-tube exit O(N^-1/10)'},
 {'id':'UF-05','class':'M','status':'ANALYTIC_PROOF_SUBMITTED_FOR_INDEPENDENT_INTAKE','claim':'At iid-uniform start Z_N weakly converges at finite horizons to dZ=(chi lambda Z-Z^3/3)ds+sqrt(2)dB; same for the two declared controls.'},
 {'id':'UF-06','class':'M','status':'ANALYTIC_CONSEQUENCE','claim':'Limiting diffusion stationary density is proportional to exp(chi lambda z^2/2-z^4/12), not global finite-N Gibbs equilibrium.'},
 {'id':'UF-07','class':'M','status':'CONDITIONAL_ACCEPTED_COEXISTENCE','claim':'Full Gibbs occupation at critical g is not concentrated near uniform as N tends to infinity; input uses earlier exact coexistence witness.'},
 {'id':'UF-08','class':'PP/EV','status':'NOT_ESTABLISHED','claim':'No physical critical gain clock preparation, kernel source, external records or experimental validation.'}
]
(p/'CLAIM_REGISTER.json').write_text(json.dumps(claims,indent=2,ensure_ascii=False)+'\n')
(p/'PROPOSED_AGENTS_PATCH.md').write_text('## Proposed scoped addition (not applied)\n\nUF-001–UF-005: under accepted q12 spectral boxes, exact LOO HB and iid uniform starts at g_N=g_c+lambda N^-1/2, an analytical proof has been submitted for a quartic critical path limit with generator (chi lambda z-z³/3)∂z+∂zz, including local fast-mode cancellation and tube-exit estimate O(N^-1/10). The frozen E04 g0-clock MET/BAR controls admit the same leading limit via a local Lipschitz logarithmic-mean argument. This remains **intake review pending** (not mechanized or peer-reviewed). The limiting diffusion invariant quartic law is not the full finite-N Gibbs equilibrium; prior exact coexistence and a nominal one-label competitor demonstrate the global mismatch. No physical A7, gain, initialization or time source is established.\n')
(p/'REPLAY.md').write_text('''# Reproduction steps\n\n1. Check `MANIFEST.sha256` before running any writer; it hashes the original frozen reports and all bundled inputs.\n2. Create an environment with Python 3, numpy and scipy (see ENVIRONMENT.json).\n3. Execute `python uf_analysis.py`, then `python independent_verify.py`. Compare numerical values to DIAGNOSTICS.json and INDEPENDENT_VERIFICATION.json within environment-appropriate floating tolerances; exact binary64 hash replay is not promised.\n4. The independent verifier reconstructs A from circulant spectral projectors without importing the producer's source model; both use the same mathematical definitions.\n5. Review UF-001 to UF-004 proof documents for analytical obligations. Script checks do not prove stopped process convergence.\n6. Do not automatically rerun `make_reports.py`, which re-authors verdicts and hashes. The previous KI ZIP and original user plan/intake are preserved unchanged.\n\nNo live remote repository access, external experiment, or Git mutation occurred.\n''')
(p/'README_PL.md').write_text('''# FIN UF-001–UF-005 — skrót\n\nNowa teoria krytyczna dotyczy **bezstronnego startu iid-uniform**, przy q=12 i g_N=g_c+lambda N^-1/2. Amplituda parzystości skaluje się N^-1/4, a czas N^1/2. Przedstawiono analityczny dowód lokalnego i skończonohoryzontowego limitu `dZ=(chi lambda Z-Z³/3)ds+sqrt(2)dB` wraz z oszacowaniem wyjścia z tuby szybkich modów O(N^-1/10). Trzy zadane kinetyki mają identyczne współczynniki wiodące w tym szczególnym punkcie; Metropolis wymaga argumentu lokalnie Lipschitzowskiego.\n\nGranica ta NIE jest globalną równowagą modelu. Niezależne źródło A7, zegar fizyczny i pomiary pozostają nieznane. Dowód jest autorskim argumentem do przeglądu, nie weryfikacją formalną/recenzją.\n''')
print('CREATED',len(modules),'tasks')
