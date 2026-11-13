"""Matched slow-time comparison and explicitly model-dependent extrapolations."""
from pathlib import Path
import sys,json
import numpy as np
H=Path(__file__).resolve().parent
P=H.parent/'fin_mixed_scaling_phase_20261113'
sys.path.insert(0,str(P))
from compare import fit,coefficients

def main():
    rows=[]
    for p in sorted(H.glob('phase*.json')):
        r=json.loads(p.read_text());a=np.loadtxt(p.with_suffix('.csv'),delimiter=',',skiprows=1)
        assert np.all(np.isfinite(a)) and np.all(a[:,1:]>0)
        assert np.all(np.diff(a[:,0])>0) and abs(a[-1,0]*.05**2-1)<1e-12
        r['fits']={str(s):{str(k):fit(a,s*r['T'],k) for k in [1,2,5,10]} for s in [.25,.5]}
        r['dense_integration_sampling']={str(s):{str(k):fit(a,s*r['T'],1,k) for k in [1,5,10]} for s in [.25,.5]}
        r['file']=p.name;rows.append(r)
    assert len(rows)==3
    old=json.loads((P/'comparison.json').read_text())['runs']
    extrap=[]
    for phase in [0.,np.pi/2]:
        selected=sorted([r for r in old+rows if r['phase']==phase and r['dx']==.1],key=lambda r:r['epsilon'])
        assert len(selected)==3
        for start in ['0.25','0.5']:
            for action in [0,1]:
                eps=np.array([r['epsilon'] for r in selected]);y=np.array([r['fits'][start]['1'][action] for r in selected])
                for power in [1,2]:
                    slope,intercept=np.polyfit(eps**power,y,1)
                    extrap.append({'phase':phase,'start_tau':float(start),'action':['I','J'][action],'assumed_power':power,'intercept':float(intercept),'relative_intercept_deviation':float(intercept/coefficients(1.5)[2]-1),'max_fit_residual':float(max(abs(y-(intercept+slope*eps**power))))})
    coarse=next(r for r in rows if r['phase']==0 and r['dx']==.1)
    fine=next(r for r in rows if r['dx']==.05)
    refinement={s:(np.array(fine['fits'][s]['1'])-np.array(coarse['fits'][s]['1'])).tolist() for s in ['0.25','0.5']}
    ac=np.loadtxt(H/coarse['file'].replace('.json','.csv'),delimiter=',',skiprows=1)
    af=np.loadtxt(H/fine['file'].replace('.json','.csv'),delimiter=',',skiprows=1)
    assert np.allclose(ac[:,0],af[:,0],rtol=0,atol=1e-12)
    trace_error=(np.max(abs(ac[:,1:]-af[:,1:]),axis=0)/af[0,1:]).tolist()
    out={'theory':coefficients(1.5).tolist(),'runs':rows,'refinement_fit_delta':refinement,'max_trace_difference_relative_initial_action':trace_error,'exploratory_extrapolations':extrap,'extrapolation_warning':'Only three amplitudes; assumed power is not derived. No confidence interval, certified limit or zero-discrepancy claim.'}
    (H/'comparison.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in rows:print(r['file'],r['fits']['0.25']['1'],r['fits']['0.5']['1'])
    print('refinement',refinement,'trace',trace_error)
if __name__=='__main__':main()
