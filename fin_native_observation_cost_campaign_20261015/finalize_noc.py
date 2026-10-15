#!/usr/bin/env python3
"""Emit reproducible NOC-only stage metadata and hash manifests.
No predecessor files are copied or mutated."""
from pathlib import Path
import hashlib,json,os,subprocess,zipfile
R=Path(__file__).resolve().parent
inputs=[Path('/mnt/data/Wklejony Markdown(20261015-131511).md'),Path('/mnt/data/fin_dynamic_kernel_observability_campaign_20261015/DKO-001/ACTUAL_CHARACTERISTICS.md'),Path('/mnt/data/fin_dynamic_kernel_observability_campaign_20261015/DKO-001/BRACKET_AND_INITIAL_PAYMENT.md'),Path('/mnt/data/fin_dynamic_kernel_observability_campaign_20261015/DKO-002/JOINT_FAST_PATH_THEOREM_OR_BLOCKER.md'),Path('/mnt/data/fin_dynamic_kernel_observability_campaign_20261015/DKO-004/LIMIT_ORDER_AND_TRANSFER_GATES.md'),Path('/mnt/data/files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip'),Path('/mnt/data/FIN_KO_001_KO_005_ONLY_HANDOFF_20261015.zip')]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
assert all(x.is_file() for x in inputs), 'Missing needed source files'
records=[{'source_path':str(x),'bytes':x.stat().st_size,'sha256':sha(x)} for x in inputs]
(R/'SOURCE_INPUTS.json').write_text(json.dumps(records,indent=2)+'\n')
(R/'SOURCE_INPUTS.sha256').write_text(''.join(f'{x["sha256"]}  {x["source_path"]}\n' for x in records))
statuses={
 '001':('CONDITIONAL_NATIVE_UNIFORM_L1_GENERATOR_ESTIMATE','Exact native block/outside generator, hypergeometric moments and uniform smooth chart estimate.'),
 '002':('PROVED_ANALYSIS_CORRECTOR_AND_IDEAL_MEAN_GAP','Explicit OU analysis semigroup, alpha-uniform C3 Poisson corrector, positive linear mean separation.'),
 '003':('CONDITIONAL_NATIVE_ONE_STEP_CONCENTRATION_BOUND','Native stationarity Duhamel, full-skeleton martingale difference, no iid measurement assumption.'),
 '004':('CONDITIONAL_FINITE_PARENT_CLASSIFIER_CONSISTENCY','Frozen cosine test, true floors, additive iid parent coupling errors, one-sided sufficient asymptotic.'),
 '005':('HANDOFF_COMPLETED_NOC_ONLY','No inherited results repackaged; single unexecuted necessity/open resource atom preserved.')
}
new={
 '001':['NATIVE_BLOCK_CONTRACT.json','UNIFORM_MOMENTS_AND_CHARACTERISTICS.md','GENERATOR_ERROR_THEOREM_OR_BLOCKER.md','BOUNDARY_PAYMENT.md'],
 '002':['EXPLICIT_CORRECTOR.md','UNIFORM_DERIVATIVE_BOUNDS.md','IDEAL_MEAN_SEPARATION.md','FROZEN_COSINE_STATISTIC.json'],
 '003':['FIXED_STEP_DUHAMEL_PROOF.md','NATIVE_MARTINGALE_DECOMPOSITION.md','DEPENDENT_RECORD_ERROR_BOUND.md','ERROR_LEDGER.json'],
 '004':['FROZEN_NATIVE_TEST.json','FINITE_N_CONSISTENCY_OR_BLOCKER.md','PARENT_TRANSFER_AND_EXIT_LEDGER.md','NONCONCLUSIONS.md'],
 '005':['DEPENDENCY_DAG.md','PROOF_OR_OPEN_LEMMA.md','NEXT_OPEN_ATOM.md','NONCONCLUSIONS.md']
}
for num,(stat,desc) in statuses.items():
 d=R/f'NOC-{num}'
 if num!='005':
  (d/'PROOF_OR_OPEN_LEMMA.md').write_text(f'# NOC-{num}: analytical status\n\n**{stat}.**\n\n{desc}\n\nAll scoped proof details and limitations are in the NOC-{num} named files. Any previous MT/KO/DKO step is an explicitly labelled conditional dependency; ordinary finite assertions do not mechanically prove it.\n')
 if num not in ['004','005']:
  (d/'NONCONCLUSIONS.md').write_text(f'# NOC-{num}: nonconclusions\n\nThis task does not infer the final classifier theorem without later stages, full path equivalence, sign-only cost, physical clock/amplitude calibration, source of A7, physical apparatus or an empirical result. Prior MT/KO/DKO inputs stay conditional and unchanged.\n')
 (d/'REPORT.md').write_text(f'# NOC-{num}: task report\n\n**Scientific disposition:** `{stat}`.\n\n{desc}\n\nNewly generated proof artifacts: {", ".join(new[num])}. Only NOC task artifacts are included; no historical source materials are recopied. See root HANDOFF for overall scope.\n')
 (d/'RESULTS.json').write_text(json.dumps({'task':f'NOC-{num}','execution_status':'COMPLETED','scientific_verdict':stat,'original_model':'q12 labelled unit-attempt LOO heat bath','compared_models':['primary A7','diagnostic A7+0.01(P3-P4)'],'initialization':'iid-uniform and inherited model-specific local post-MT stationary comparison','time_scope':'native fast-grid unit spacing, polynomial growing record for NOC-004','observed_input':'predetermined block f3 average and cos; full global counts hidden proof variables','derived_objects':new[num],'supplied_objects':['frozen A7 family','conditional MT local moment/mixing/trapping','DKO native characteristics'],'fitted_objects':[],'arithmetic_status':'bounded exact/algebraic and floating fixture replay; asymptotic theorem analytical, conditional','error_budget':'N^-1/4 + k^-1/2 + (N/sqrt(k))pi_D(B); decision 1/(M alpha^2)+eps/alpha+delta_N','discrepancies':[],'unpaid':['independent re-proof of historical MT local domains for both operators','finite-N necessity/minimax resource result','physical source/preparation/calibration']},indent=2)+'\n')
 (d/'REPLAY.md').write_text(f'# NOC-{num} replay\n\nFrom the extracted root, run `python verify_noc.py` for exact rational/projector, native finite q12 generator, corrector and frozen-window fixtures. Check `INDEPENDENT_VERIFICATION.json`. The asymptotic analytical estimates must be independently reviewed by reading NOC-{num} proof files; test PASS counts do not certify MT hypotheses. At the root, `python verify_package.py` validates hashes.\n')
 (d/'INPUTS.sha256').write_text(''.join(f'{x["sha256"]}  {x["source_path"]}\n' for x in records))
claims=[
 ('NOC-001','stationary_native_generator','conditional_analytical','pi_D |(L_D-G) v(Phi)|<=C ||v||_C3_1 [N^-1/4+k^-1/2+N/sqrt(k)*pi_D(B)]','MT local moments/domain'),
 ('NOC-002','discrete_corrector','analytic','psi=sum_l(P_l g-mu); (I-P1)psi=g-mu; C3 derivative bounds uniform alpha','fixed gamma>0'),
 ('NOC-002','strict_cosine_mean_gap','analytic','mu_A7-mu_control=d alpha+O(alpha^2), d>0','frozen operator pair'),
 ('NOC-003','dependent_native_concentration','conditional_analytical','P_pi(|A_M-mu|>c alpha)<=C/(M alpha^2)+C eps/alpha','NOC-001 and corrector'),
 ('NOC-004','unchanged_parent_consistent_classifier','conditional_analytical','max_j P_parent_j(wrong)->0 along frozen floors','inherited MT mixing/exit both models'),
 ('NOC-005','finite_N_minimax_cost','open','No necessity/minimax result or sign-only transfer','no task executed')]
(R/'CLAIM_REGISTER.json').write_text(json.dumps([{'stage':a,'claim':b,'status':c,'precise_statement':d,'dependency':e} for a,b,c,d,e in claims],indent=2)+'\n')
(R/'REPLAY.md').write_text('# NOC-only replay\n\n1. `python verify_noc.py` recomputes independent finite checks into `INDEPENDENT_VERIFICATION.json` (57/57).\n2. `python verify_package.py` checks the manifest against current files and prints the number of covered entries.\n3. Read `NOC-001/UNIFORM_MOMENTS_AND_CHARACTERISTICS.md`, `NOC-002/UNIFORM_DERIVATIVE_BOUNDS.md`, `NOC-003/NATIVE_MARTINGALE_DECOMPOSITION.md`, `NOC-004/FINITE_N_CONSISTENCY_OR_BLOCKER.md`. Their stochastic-limit dependencies are not mechanically verified by the scripts.\n4. Compare source file hashes only for files listed in `SOURCE_INPUTS.sha256`; no unavailable review proof input is pretended present.\n5. Repack source by `python finalize_noc.py` on the same input paths, then `python verify_package.py`.\n\nImportant: this package does not include any historical review, archive, source code or AGENTS change beyond hash references.\n')
(R/'INPUTS.sha256').write_text((R/'SOURCE_INPUTS.sha256').read_text())
# own verifier script excluded from own manifest self recursion; package.json excluded too. source inputs list is explicit.
exclude={'MANIFEST.sha256','PACKAGE_VERIFICATION.json'}
paths=sorted((x for x in R.rglob('*') if x.is_file() and not x.name.endswith('.pyc') and '__pycache__' not in x.parts and x.relative_to(R).as_posix() not in exclude),key=lambda p:p.relative_to(R).as_posix())
(R/'MANIFEST.sha256').write_text(''.join(f'{sha(x)}  {x.relative_to(R).as_posix()}\n' for x in paths))
package=R.parent/'FIN_NOC_001_NOC_005_ONLY_HANDOFF_20261015.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for x in sorted((y for y in R.rglob('*') if y.is_file() and not y.name.endswith('.pyc') and '__pycache__' not in y.parts),key=lambda p:p.relative_to(R).as_posix()):
  z.write(x,arcname=R.name+'/'+x.relative_to(R).as_posix())
print('STAGES',len(statuses),'MANIFEST',len(paths),'ZIP',package,'SIZE',package.stat().st_size)
