"""Numerical integrity and replay checks, not a theorem verifier."""
from pathlib import Path
import hashlib,json,subprocess,sys
import numpy as np
H=Path(__file__).resolve().parent
for filename,expected in json.loads((H/'INPUTS.json').read_text()).items():
    assert hashlib.sha256((H.parent/filename).read_bytes()).hexdigest()==expected,filename
for filename,expected in json.loads((H/'MANIFEST.json').read_text()).items():
    assert hashlib.sha256((H/filename).read_bytes()).hexdigest()==expected,filename
previous=(H/'comparison.json').read_bytes()
subprocess.run([sys.executable,str(H/'analyze.py')],check=True)
assert (H/'comparison.json').read_bytes()==previous
out=json.loads(previous);assert len(out['runs'])==3
for r in out['runs']:
    a=np.loadtxt(H/r['file'].replace('.json','.csv'),delimiter=',',skiprows=1)
    assert np.all(np.isfinite(a)) and np.all(a[:,1:]>0)
    assert len(a)==4001 and abs(a[-1,0]-400)<1e-10
    assert np.all(np.diff(a[:,0])>0) and np.max(np.diff(a[:,0]))<.101
    assert r['L']>=r['T']+30 and r['relative_global_energy_drift']<1e-5
    assert max(r['max_relative_local_balance_residual'])<1e-5
assert len(out['exploratory_extrapolations'])==16
assert max(out['max_trace_difference_relative_initial_action'])<.001
print('PASS: dependency/output hashes, 3 runs, matched times, energy/flux controls, trace refinement and deterministic analysis replay')
