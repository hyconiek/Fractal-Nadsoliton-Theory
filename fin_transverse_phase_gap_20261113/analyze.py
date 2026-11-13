from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent
r=json.loads((H/'results.json').read_text());assert len(r)==9
for v in r:
 a=np.loadtxt(H/(v['name']+'.csv'),delimiter=',',skiprows=1)
 assert a.shape==(201,5) and np.isfinite(a).all() and np.all(a[:,1]>0)
 assert np.max(abs(a[0,1:]-1))<1e-12
 assert abs(max(a[:,2])-v['max_energy_gain'])<1e-12
assert r[0]['static_max_deviation']<1e-8
pairs=[]
for w in [1.25,1.5]:
 for p in [0.,np.pi/2]:
  z=[v for v in r if v['driven'] and v['w']==w and v['phase']==p];assert len(z)==2
  c,f=sorted(z,key=lambda v:-v['dx']);a=np.loadtxt(H/(c['name']+'.csv'),delimiter=',',skiprows=1);b=np.loadtxt(H/(f['name']+'.csv'),delimiter=',',skiprows=1)
  err=float(max(abs(a[:,2]-b[:,2])));assert err<1e-3
  pairs.append(dict(w=w,phase=p,max_gain=f['max_energy_gain'],final_gain=f['final_max_energy_gain'],max_gain_trace_refinement_difference=err))
(H/'comparison.json').write_text(json.dumps(pairs,indent=2)+'\n');print(json.dumps(pairs,indent=2))
