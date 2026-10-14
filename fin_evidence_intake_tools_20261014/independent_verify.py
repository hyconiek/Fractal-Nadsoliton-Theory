#!/usr/bin/env python3
"""Small independent verifier for READY-001 artifacts.
Does not import strict_records/strict_decision/strict_process_gate.
"""
from pathlib import Path
import hashlib,json,subprocess,sys,math
import mpmath as mp
ROOT=Path(__file__).resolve().parent
checks=[]
def add(name,ok,detail=None):checks.append({'name':name,'pass':bool(ok),'detail':detail})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Source binding hashes.
b=json.loads((ROOT/'PROCESS_BINDING.json').read_text())
for fn,h in b['source_hashes'].items():
    p=ROOT/'source_context'/fn
    add('hash_'+fn,p.exists() and sha(p)==h,sha(p) if p.exists() else 'missing')
# Unit test result exists and is complete.
u=json.loads((ROOT/'UNIT_TESTS.json').read_text())
add('unit_tests_status',u['status']=='PASS',(u['passed'],u['total']))
add('unit_tests_all_pass',u['passed']==u['total'] and u['total']>=40,u['failed'])
# Exact high-precision threshold recheck, separate formula implementation.
mp.mp.dps=100
rho=mp.mpf('0.00119114'); eps=rho/3; alpha=mp.mpf('0.025'); k=13
x=mp.log((2**k-2)/alpha)/(2*eps**2)
n0=int(mp.ceil(x))
add('threshold_40278676',n0==40278676,str(x))
def rad(n):return mp.sqrt(mp.log((2**k-2)/alpha)/(2*n))
add('threshold_minimality',rad(n0)<=eps and rad(n0-1)>eps,{'r_n':str(rad(n0)),'r_prev':str(rad(n0-1)),'eps':str(eps)})
rho2=mp.mpf('0.00039114'); eps2=rho2/3
x2=mp.log((2**k-2)/alpha)/(2*eps2**2);n2=int(mp.ceil(x2))
add('illustrative_threshold_373539331',n2==373539331,str(x2))
# CLI behavior checks without importing modules.
def run(args):
    return subprocess.run([sys.executable,*args],cwd=ROOT,text=True,capture_output=True,timeout=60)
r=run(['strict_records.py','fixtures/valid.csv'])
add('cli_valid_records_pass',r.returncode==0 and 'PASS_SYNTAX_INTEGRITY_ONLY' in r.stdout,(r.returncode,r.stdout[-300:],r.stderr))
r=run(['strict_records.py','fixtures/nan_elapsed.csv'])
add('cli_nan_records_fail_closed',r.returncode==2 and 'FAIL_CLOSED' in r.stdout,(r.returncode,r.stdout[-300:],r.stderr))
r=run(['strict_records.py','fixtures/role_collision.csv'])
add('cli_role_collision_fail_closed',r.returncode==2 and 'FAIL_CLOSED' in r.stdout,(r.returncode,r.stdout[-300:],r.stderr))
for model in ['FIN_A7','FULL_POTTS_TRACE','FLAT_P7_TRACE']:
  for kin in ['heat_bath','metropolis','barker']:
    args=['strict_process_gate.py','--model',model,'--kinetic',kin,'--gain','3']
    if kin=='heat_bath':args+=['--shortcut']
    r=run(args);add(f'process_cli_{model}_{kin}',r.returncode==0 and 'PASS_SOFTWARE_CONTRACT_ONLY' in r.stdout,r.returncode)
# Invalid decision n through CLI.
(ROOT/'fixtures/decision_bad_n.json').write_text(json.dumps({'n':'NaN','counts':[0]*13})+'\n')
r=run(['strict_decision.py','fixtures/decision_bad_n.json'])
add('decision_cli_nan_n_fail_closed',r.returncode==2 and 'FAIL_CLOSED' in r.stdout,(r.returncode,r.stdout))
status='PASS' if all(c['pass'] for c in checks) else 'FAIL'
out={'status':status,'passed':sum(c['pass'] for c in checks),'total':len(checks),'checks':checks,'scope':'independent software arithmetic/replay checks only; no physical calibration/provenance certification'}
(ROOT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':status,'passed':out['passed'],'total':out['total'],'failed':[c['name'] for c in checks if not c['pass']]},indent=2))
raise SystemExit(0 if status=='PASS' else 1)
