#!/usr/bin/env python3
from pathlib import Path
import os, subprocess, sys, json, time
ROOT=Path(__file__).resolve().parent
files=[
 'tests/test_B.py','tests/test_C_face.py','tests/test_D_coexistence.py','tests/test_D_fold.py',
 'tests/test_E_edge_reflection.py','tests/test_E_stationary_counterexample.py','tests/test_F_boundary_ising.py',
 'tests/test_G_boundary_cover.py','tests/test_G_checker.py','tests/test_H_intraparity.py','tests/test_H_intraparity_closure.py',
 'tests/test_I_local_cone.py','tests/test_I_off_face.py','tests/test_I_tail_partial.py','tests/test_K_phase_cumulants.py',
 'tests/test_M_global_frontier.py','tests/test_N_passive_finiteN.py','tests/test_O_quantum_bridge.py','tests/test_schema_and_resume.py']
env=os.environ.copy(); env['PYTHONDONTWRITEBYTECODE']='1'; env['PYTHONPATH']='.:..:src'
rows=[]
for f in files:
 t=time.time(); p=subprocess.run([sys.executable,'-m','pytest','-q',f],cwd=ROOT,env=env,text=True,capture_output=True)
 rows.append({'file':f,'returncode':p.returncode,'seconds':time.time()-t,'stdout':p.stdout,'stderr':p.stderr})
 print(f, 'PASS' if p.returncode==0 else 'FAIL', flush=True)
 if p.returncode:
  print(p.stdout); print(p.stderr,file=sys.stderr)
summary={'rows':rows,'passed_files':sum(r['returncode']==0 for r in rows),'failed_files':sum(r['returncode']!=0 for r in rows)}
(ROOT/'logs/sharded_replay_latest.json').write_text(json.dumps(summary,indent=2)+'\n')
raise SystemExit(0 if summary['failed_files']==0 else 1)
