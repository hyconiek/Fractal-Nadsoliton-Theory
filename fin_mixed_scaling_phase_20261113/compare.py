"""Separate integration-grid sensitivity from regression sampling on fixed nodes."""
from pathlib import Path
import sys,json
import numpy as np
H=Path(__file__).resolve().parent
sys.path.insert(0,str(H.parent/'fin_pde_two_mode_20261108'))
from analyze import coefficients

def fit(a,start,stride=1,node_stride=10):
    a=a[a[:,0]>=start-1e-9];idx=np.unique(np.r_[np.arange(0,len(a),stride),len(a)-1])
    t,I,J=a[idx].T
    def integral(v):return np.r_[0,np.cumsum(np.diff(t)*(v[1:]+v[:-1])/2)]
    ka,kb,kx=coefficients(1.5)
    X=integral(I*J);YI=I[0]-I-ka*integral(I*I);YJ=J[0]-J-kb*integral(J*J)
    # Identical regression nodes for every integration-grid comparison.
    nodes=np.unique(np.r_[np.arange(0,len(a),node_stride),len(a)-1]);tt=a[nodes,0]
    xx=np.interp(tt,t,X)
    return [float(np.polyfit(xx,np.interp(tt,t,y),1)[0]) for y in [YI,YJ]]

def main():
    rows=[]
    for p in sorted(H.glob('eps*.json')):
        r=json.loads(p.read_text());a=np.loadtxt(p.with_suffix('.csv'),delimiter=',',skiprows=1)
        assert np.all(np.isfinite(a)) and np.all(np.diff(a[:,0])>0)
        assert abs(a[-1,0]*r['epsilon']**2-1)<1e-12
        assert np.all(a[:,1:]>0)
        row={**r,'sample_count':len(a),'fits':{}}
        for start in [.25,.5]:
            row['fits'][str(start)]={str(k):fit(a,start*r['T'],k) for k in [1,2,5,10]}
        row['dense_integration_sampling']={str(start):{str(k):fit(a,start*r['T'],1,k) for k in [1,5,10]} for start in [.25,.5]}
        rows.append(row)
    assert len(rows)==4
    out={'theory':coefficients(1.5).tolist(),'runs':rows}
    (H/'comparison.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in rows:print(r['epsilon'],r['phase'],r['fits']['0.25']['1'],r['fits']['0.5']['1'],r['relative_global_energy_drift'])
if __name__=='__main__':main()
