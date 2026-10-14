#!/usr/bin/env python3
"""Author research status artifacts and a byte manifest; no scientific results are fabricated."""
import json,hashlib,sys,platform,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
sh=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=[root/'source_context'/'SF_PLAN_USER.md',root/'source_context'/'TH_INTAKE_USER.md',root/'source_context'/'FIN_TH_001_TH_005_FULL_HANDOFF_20261014.zip',root/'source_context'/'TH005_PREDECESSOR_STOCHASTIC_FOLD_LIMIT.md',root/'source_context'/'TH004_PREDECESSOR_SOFT_MODE_BRIDGE.md',root/'source_context'/'FOLD_SOURCE_BINDING.json']
(root/'SOURCE_INPUTS.sha256').write_text(''.join(sh(p)+'  '+str(p.relative_to(root))+'\n' for p in sources))
refs={'SF-001':['source_context/SF_PLAN_USER.md','source_context/TH_INTAKE_USER.md','source_context/FOLD_SOURCE_BINDING.json','SF-001/FAST_CHART_PROOF.md'],
 'SF-002':['source_context/SF_PLAN_USER.md','SF-001/FAST_CHART_PROOF.md','SF-002/GENERATOR_ENERGY_PROOF.md'],
 'SF-003':['source_context/SF_PLAN_USER.md','SF-002/GENERATOR_ENERGY_PROOF.md','SF-003/STOPPED_CONVERGENCE_THEOREM.md'],
 'SF-004':['SF-003/STOPPED_CONVERGENCE_THEOREM.md','SF-004/FIRST_PASSAGE_PREDICTION.md','SF-004/EXIT_SIDE_CERTIFICATE.json'],
 'SF-005':['source_context/SF_PLAN_USER.md','SF-005/HANDOFF.md','SF-005/NEXT_OPEN_ATOM.md']}
statuses={
 'SF-001':('ANALYTIC_CHART_PROOF','One soft and ten stable occupation directions; positive derived P; prepared actual rounded types','exact symbolic identity under inherited validated fold root/inertia'),
 'SF-002':('LOCAL_ENERGY_PROOF','Actual LOO generator local Lyapunov bound, O(N^-2/3) stopped occupation energy and O(N^-1/6) rare fast exit','analytic asymptotic uniform-on-frozen-compact proof; no numerical N0'),
 'SF-003':('LOCAL_WEAK_LIMIT_PROOF_SUBMITTED','C-tight stopped Z_N to nonlinear diffusion plus boundary/exit mapping by outer-localization','analytic weak-convergence proof submitted for independent review; inherited fold interval premise'),
 'SF-004':('ANALYTIC_EXIT_LAW_CONDITIONAL_ON_SF003','Canonical dimensionless exit-side scale integral and complete finite-horizon BVP/Laplace law; rationally bounded eventual upper exit','exact rational interval for canonical event, no certified finite-horizon numeric values'),
 'SF-005':('SCOPED_HANDOFF','Report all paid and unpaid claims, unchanged physical boundaries, and architect checkpoint','reporting/meta; no new scientific computation')}
scope='Fixed lambda in bounded compact, soft R/S finite, q12, prepared rounded actual counts; local stopped process; inherited exact isolated fold.'
for code,(status,conclusion,arith) in statuses.items():
    d=root/code
    obj={'task':code,'execution_status':'COMPLETED','scientific_verdict':status,'supplied_objects':['fixed exact-decimal X/A, N-copy leave-one-out jump generator, per-copy model clock','inherited exact isolated simple-fold root/inertia and fixed sign'],
         'derived_objects':[conclusion],'fitted_objects':[],'proof_scope':scope,'arithmetic_status':arith,
         'initial_preparation':'deterministic balanced rounded counts from pf+N^(-1/3)z0 r; compact z0 strictly inside (-R,R); optional uniformly prepared random extension',
         'stop_scope':'first |Z|>=R, or ||w||>=N^(-1/4), or scaled horizon S; fast exit negligible under conditions',
         'dependencies':refs[code],'discrepancies':[],'next_action':'submit proof and exact-fold-source dependence for independent mathematical review; architect selects subsequent direction' if code in ('SF-003','SF-005') else 'follow next numbered task without altering parent inputs',
         'is_empirical_result':False,'is_physical_source_result':False}
    (d/'RESULTS.json').write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
    (d/'REPORT.md').write_text('# '+code+' — report\n\n**'+status+'**. '+conclusion+'\n\nScope: '+scope+'\n\nRead the detailed analytic object(s) listed in RESULTS.json. Finite floating diagnostics do not certify the asymptotic proof.\n')
    (d/'PROOF_OR_OPEN_LEMMA.md').write_text('# '+code+' — proof status\n\n'+('Proof is given in `FAST_CHART_PROOF.md`, using the accepted simple-fold/inertia result.' if code=='SF-001' else
            'Proof with all jump terms and stopped estimates is in `GENERATOR_ENERGY_PROOF.md`.' if code=='SF-002' else
            'Proof is in `STOPPED_CONVERGENCE_THEOREM.md` and `BOUNDARY_AND_TOPOLOGY.md`. Requires independent expert review and inherited exact fold certificate.' if code=='SF-003' else
            'Exact exit probability integral and canonicalized finite-horizon laws in `FIRST_PASSAGE_PREDICTION.md`; rational series certificate is independently checked. Finite-horizon numeric enclosure not supplied.' if code=='SF-004' else
            'No new mathematics in reporting task. Architect review remains open.')+'\n')
    (d/'NONCONCLUSIONS.md').write_text('# '+code+' — excluded interpretations\n\nNo independent A7 source or new kinetic law, no physical clock, no global mixing/transition/prefactor, no empirical validation, no spatial/QM/GR/SM/ToE derivation. Uncertified binary64 diagnostic values are not exact-root enclosures.\n')
    (d/'REPLAY.md').write_text('# '+code+' — replay\n\nFrom package root, follow `REPLAY.md` (primary). Detailed proof is human-reviewable Markdown; scripts check numerical/algebraic prerequisites but cannot replace an asymptotic theorem proof.\n')
    (d/'INPUTS.sha256').write_text(''.join(sh(root/x)+'  '+('../'+x if not x.startswith(code+'/') else x[len(code)+1:])+'\n' for x in refs[code]))
claim=[
 {'id':'SF-001-CHART','type':'exact analytical','status':'PROOF_SUBMITTED','premise':'certified H7 simple fold','claim':'full occupation tangent one soft and ten stable, oblique Q, positive Lyapunov P'},
 {'id':'SF-001-PREP','type':'exact analytical','status':'PROOF_SUBMITTED','premise':'rounded true N-copy types, z0 interior','claim':'w0=O(1/N), Z0->z0'},
 {'id':'SF-002-ENERGY','type':'analytic asymptotic','status':'PROOF_SUBMITTED','premise':'exact LOO rates + stopped prepared class','claim':'stopped energy O(N^-2/3) model time and fast tube exit O(N^-1/6)'},
 {'id':'SF-003-FOLD','type':'weak-convergence theorem','status':'PROOF_SUBMITTED_PENDING_INDEPENDENT_REVIEW','premise':'SF-001 and SF-002 plus inherited interval fold hypothesis','claim':'compact-stopped process Z_N => nonlinear diffusion in Skorokhod J1'},
 {'id':'SF-004-CANON','type':'exact analytic','status':'PROOF_SUBMITTED','premise':'SF-003 conditional model limit','claim':'canonical dY=(kappa lambda -Y^2) dt +dW and finite-horizon BVP laws'},
 {'id':'SF-004-EXIT-SIDE','type':'exact rational interval','status':'CHECKED','premise':'canonical lambda=0,start0,bounds±1','claim':'P_upper∈[0.416855651698382283,0.416855651698382284] at eventual exit'},
 {'id':'NUMERICAL-CHART','type':'binary64 nominal','status':'CHECKED','premise':'nominal fold and reconstructed kernel','claim':'numerical eigengap/P and model identities match within floating tolerance'},
 {'id':'PHYSICAL_SOURCE','type':'physical premise','status':'NO_NEW_SOURCE','claim':'No independent physical A7/clock/units or experimental data'}]
(root/'SF-005'/'CLAIM_REGISTER.json').write_text(json.dumps(claim,indent=2)+'\n')
(root/'SF-005'/'DEPENDENCY_DAG.json').write_text(json.dumps({'nodes':['supplied A7/LOO rates','accepted exact fold/inertia','SF-001 chart/starts','SF-002 generator energy','SF-003 stochastic limit','SF-004 first passage','SF-005 architect checkpoint'], 'edges':[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6]], 'no_automatic_next_campaign':True},indent=2)+'\n')
(root/'ENVIRONMENT.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform()},indent=2)+'\n')
manifest=[]
for p in sorted(root.rglob('*')):
    if p.is_file() and p.name!='MANIFEST.sha256' and not '__pycache__' in p.parts:
        manifest.append(sh(p)+'  '+str(p.relative_to(root)))
(root/'MANIFEST.sha256').write_text('\n'.join(manifest)+'\n')
print('manifest entries',len(manifest))
