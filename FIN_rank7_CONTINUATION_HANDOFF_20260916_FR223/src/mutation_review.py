from __future__ import annotations
import subprocess,sys,time,json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CASES=[
 ('boundary_tree_mutations','tests/test_G_checker.py','altered spectral interval, corrupt split, flipped reason, removed leaf plus valid replay'),
 ('schema_and_input_hash_mutations','tests/test_schema_and_resume.py','bad exact endpoint, unresolved global pass, input hash mutation plus resume'),
 ('local_cone_negative_radius','tests/test_I_local_cone.py','accepted rho=1/8192 and rejection of rho=1/4096 by the same checker')]

def main():
 env=os.environ.copy();env['PYTHONPATH']='.:..:src';env['PYTHONDONTWRITEBYTECODE']='1'; rows=[]; logs=[]
 for name,test,meaning in CASES:
  t=time.time();cp=subprocess.run([sys.executable,'-m','pytest','-q',test],cwd=ROOT,env=env,text=True,capture_output=True,timeout=60);sec=time.time()-t
  rows.append({'name':name,'test':test,'returncode':cp.returncode,'seconds':sec,'meaning':meaning});logs.append(f'== {name} ==\n{cp.stdout}{cp.stderr}')
 out={'status':'PASS' if all(x['returncode']==0 for x in rows) else 'FAIL','cases':rows,
      'independence_level':{'boundary_checker':'independent replay of frozen proof spec/tree and mutation rejection; shares the mathematical polynomial specification, so not a proof-assistant-independent formalization','schema':'independent structural validation and copied-package input mutation','local_cone':'same interval checker with an explicit larger-radius negative control'}}
 (ROOT/'results/R7P-122_mutation_review.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/'logs/R7P-122_mutation_review.log').write_text('\n'.join(logs)+'\n');print(json.dumps(out,indent=2));return 0 if out['status']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
