from pathlib import Path
import hashlib,json
import numpy as np
from compare import fit,coefficients
H=Path(__file__).resolve().parent
for fn,expected in json.loads((H/'INPUTS.json').read_text()).items():
    assert hashlib.sha256((H.parent/fn).read_bytes()).hexdigest()==expected,fn
for fn,expected in json.loads((H/'MANIFEST.json').read_text()).items():
    assert hashlib.sha256((H/fn).read_bytes()).hexdigest()==expected,fn
out=json.loads((H/'comparison.json').read_text());assert len(out['runs'])==4
for r in out['runs']:
    p=H/f"eps{r['epsilon']}_phase{r['phase']:.6f}.csv"
    a=np.loadtxt(p,delimiter=',',skiprows=1)
    assert np.all(np.isfinite(a)) and np.all(a[:,1:]>0)
    assert np.all(np.diff(a[:,0])>0) and np.max(np.diff(a[:,0]))<.101
    assert abs(a[-1,0]*r['epsilon']**2-1)<1e-12
    assert r['L']>=r['T']+30 and r['relative_global_energy_drift']<1e-5
    assert max(r['max_relative_local_balance_residual'])<1e-5
    for start in [.25,.5]:
        for stride in [1,2,5,10]:
            assert np.allclose(fit(a,start*r['T'],stride),r['fits'][str(start)][str(stride)],rtol=1e-12,atol=1e-14)
# Synthetic exact common mixed decay: I=J=u, ka=kb is not assumed.
# Instead verify integration convergence and regression recovery for each
# action's exactly known effective slope in the constant-ratio u'= -u^2 law.
t=np.linspace(0,10,10001);u=.1/(1+.1*t);a=np.column_stack([t,u,u])
expected=1-coefficients(1.5)[:2]
assert np.allclose(fit(a,0),expected,rtol=1e-7)
print('PASS: hashes, four matched-time runs, physical energy controls, fit replay and synthetic estimator control')
