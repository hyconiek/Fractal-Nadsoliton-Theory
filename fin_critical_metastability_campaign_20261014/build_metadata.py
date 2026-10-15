#!/usr/bin/env python3
"""Authorial reports and manifests, not a scientific-prover script."""
from pathlib import Path
import json,hashlib,os
ROOT=Path(__file__).resolve().parent
stages=[
('MT-001','LOCAL_CONNECTED_SUBLEVEL_PROOF_SUBMITTED',['MT-001_COERCIVITY_CONNECTIVITY.md'],['source_context/MT_PLAN_USER.md','source_context/UF_INTAKE_USER.md','source_context/theory_core.py'],
 'Under inherited gamma>0.0186 and exact gc, rational epsilon=.01, delta=.0001, h=1e-10 give a cylinder-separated sublevel. Discrete gradient descent and a central logarithmic corridor connect all local types for large N in 12Z. IID-uniform start is covered.'),
('MT-002','ACTUAL_IID_TRAPPING_THEOREM_SUBMITTED',['MT-002_IID_TRAPPING.md'],['MT-001_COERCIVITY_CONNECTIVITY.md','source_context/UF_INTAKE_USER.md'],
 'Reversible invariant weight W>=iid-uniform permits semigroup domination, yielding P(exit<=T)<=mu0(Dc)+2NT(N+1)^12 exp(-Nb), c<b.'),
('MT-003','JOINT_CONDITIONAL_GIBBS_LIMIT_SUBMITTED',['MT-003_CONDITIONAL_GIBBS.md'],['MT-001_COERCIVITY_CONNECTIVITY.md','source_context/UF_INTAKE_USER.md'],
 'Conditional count Gibbs, not full Gibbs, has a quartic soft factor times a Gaussian in all ten transverse directions; tails and Stirling lattice are stated explicitly.'),
('MT-004','POLYNOMIAL_LOCAL_MIXING_WINDOW_PROOF_SUBMITTED',['MT-004_REFLECTED_MIXING_WINDOW.md'],['MT-001_COERCIVITY_CONNECTIVITY.md','MT-002_IID_TRAPPING.md','MT-003_CONDITIONAL_GIBBS.md'],
 'Canonical-path energy ceiling O(log N/N) implies gap>=poly(N)^-1, allowing t_loc=poly(N) and parent TV closeness to pi_D up to exp(cN). Critical upper gap <=C N^-1/2.'),
('MT-005','HANDOFF_COMPLETE_REVIEW_REQUIRED',['HANDOFF.md','NEXT_OPEN_ATOM.md','NONCONCLUSIONS.md','RESEARCH_SUMMARY_PL.md'],['MT-001_COERCIVITY_CONNECTIVITY.md','MT-002_IID_TRAPPING.md','MT-003_CONDITIONAL_GIBBS.md','MT-004_REFLECTED_MIXING_WINDOW.md'],
 'Full evidence and dependency classification; unresolved independence of review distinguished from author-submitted proof.')]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
claims={
 'MT-001':{'class':'M','status':'analytic_proof_submitted','objects':['local_F_cylinder','connected_count_sublevel','iid_start_cover','boundary_height'], 'numerics':'rational_geometry_conditional_on_inherited_exact_spectral_premise'},
 'MT-002':{'class':'M','status':'analytic_proof_submitted','objects':['finite_N_weight_domination','actual_exit_compensator','one_sided_exponential_residence'], 'numerics':'no_exit_lifetime_measurement'},
 'MT-003':{'class':'M','status':'analytic_proof_submitted','objects':['conditional_count_Gibbs_limit','quartic_soft_Gaussian_fast','moment_UI'], 'numerics':'nominal_point_probes_not_moment_proof'},
 'MT-004':{'class':'M','status':'analytic_proof_submitted_review_sensitive','objects':['local_root_paths','congestion_lower_gap','parent_TV_window'], 'numerics':'no_certified_K_or_N0; discrete_path_audit_required'},
 'MT-005':{'class':'S','status':'handoff_ready_for_independent_intake','objects':['reproducible_package','claim_scope','next_audit_atom'], 'numerics':'small_checks_only'} }

(ROOT/'CLAIM_REGISTER.json').write_text(json.dumps({'campaign':'MT-001–MT-005', 'scientific_approval':'not_independently_accepted', 'claims':claims,
 'physical_claims_admitted':[], 'future_decision':'architect_review', 'source_kernel_choice':'supplied_not_derived'},indent=2)+'\n')
versions={'python':'3','numpy':None,'scipy':None}
try:
 import numpy,scipy
 versions['numpy']=numpy.__version__;versions['scipy']=scipy.__version__
except ImportError:pass
inputs=['source_context/MT_PLAN_USER.md','source_context/UF_INTAKE_USER.md','source_context/UF_PREDECESSOR.zip','source_context/theory_core.py']
(ROOT/'PROVENANCE.json').write_text(json.dumps({'environment':versions,'repo_checkout':'not_claimed','inherited_review_spectrum':'gamma>0.0186; exact gc; see UF intake','inputs':[{'path':p,'sha256':sha(ROOT/p),'bytes':(ROOT/p).stat().st_size} for p in inputs],
 'external_reference':'canonical paths finite reversible chains; https://link.springer.com/article/10.1007/s00440-023-01257-x',
 'source_capture_scope':'user-uploaded plan and intake, original UF source archive and inherited source-context core; no original MT review'},indent=2)+'\n')
for stage,status,docs,deps,abstract in stages:
 p=ROOT/stage;p.mkdir(exist_ok=True)
 (p/'REPORT.md').write_text('# '+stage+' — '+status+'\n\n'+abstract+'\n\nPrimary argument: '+', '.join('../'+x for x in docs)+'.\n\nThis is a model-theorem **submission** under inherited input premises, not a certified physical observation, external peer review or proof-assistant check.\n')
 (p/'PROOF_OR_OPEN_LEMMA.md').write_text('# '+stage+' proof status\n\n'+abstract+'\n\n'+('**Required independent audit:** exact all-type descent/central-corridor energy ceiling; loss of this statement invalidates the polynomial gap/window.\n' if stage in ('MT-001','MT-004') else '**Remaining status:** separate peer/intake review, no physical input or numerical N0.\n'))
 (p/'NONCONCLUSIONS.md').write_text('# '+stage+' nonconclusions\n\nNo independent source of A7, physical time, measured metastable lifetime, global mixing, capacity/prefactor, exact finite-N parity bit, selected sign, autonomous particle, locality, QM/GR/ToE. Source generator and gain sequence remain modelling premises.\n')
 (p/'REPLAY.md').write_text('# '+stage+' replay\n\nRun `python ../mt_analysis.py && python ../independent_verify.py` from this task directory or run both from campaign root. Code corroborates a finite list of rational/numerical assertions only. Consult the written analytical proof in '+', '.join('../'+x for x in docs)+'.\n')
 (p/'RESULTS.json').write_text(json.dumps({'task':stage,'execution_status':'COMPLETED_AUTHOR_SUBMITTED','scientific_verdict':status,
 'supplied_objects':['q12 exact A7 kernel','g_c exact spectral root','LOO HB count generator','iid uniform initial law','fixed g0 MET/BAR controls','accepted relative gap gamma>0.0186'],
 'derived_objects':docs, 'fitted_objects':[], 'proof_scope':{'N':'sufficiently large multiples of 12','lambda':'bounded fixed interval',
 'domain':'|m|<1/100, ||w||<1/10000, F<=1e-10','time':'model time up to exp(cN) c<5e-11'},
 'arithmetic_status':'analytic inequalities + exact Fraction geometry conditional on inherited spectral premise + binary64 diagnostics; no finite-N N0 or gap coefficient certified',
 'input_dependencies':deps,'discrepancies':[], 'independent_review':'OPEN','next_action':'Independent MT discrete-path proof audit; architect decision after intake'},indent=2)+'\n')
 hlines=[]
 for x in deps:
  f=ROOT/x;assert f.is_file(),f
  hlines.append(sha(f)+'  ../'+x)
 (p/'INPUTS.sha256').write_text('\n'.join(hlines)+'\n')
(ROOT/'REPORT.md').write_text('# MT 001–005 report\n\n'+ '\n\n'.join(f'**{stage}:** {abstract}' for stage,status,docs,deps,abstract in stages) + '\n\n**Final status:** LOCAL_CRITICAL_WINDOW_PROOF_SUBMITTED_FOR_INTAKE. See `NEXT_OPEN_ATOM.md` for the discrete canonical-path lemma that must survive independent audit.\n')
print('Authored task metadata for',len(stages),'tasks')
