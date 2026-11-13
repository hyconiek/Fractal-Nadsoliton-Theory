from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent
rows=json.loads((H/'results.json').read_text());assert len(rows)==3
out=[];tr=[]
for r in rows:
 a=np.loadtxt(H/(r['name']+'.csv'),delimiter=',',skiprows=1);assert a.shape==(401,9) and np.isfinite(a).all()
 residual=a[:,3:5]-1-a[:,5:7]-a[:,7:9]
 err=float(np.max(abs(residual)));assert err<1e-4
 old=np.loadtxt(H.parent/'fin_transverse_horizon_20261113'/(r['name']+'.csv'),delimiter=',',skiprows=1)
 assert np.max(abs(a[:,:5]-old))<1e-10
 snap=[]
 for t in [100,116.5,150,200]:
  v=a[round(2*t)];snap.append(dict(t=t,gain=v[3:5].tolist(),phi_work=v[5:7].tolist(),psi_work=v[7:9].tolist()))
 out.append(dict(name=r['name'],max_normalized_balance_residual=err,samples=snap));tr.append(a)
assert np.max(abs(tr[0][:,5:]))<1e-12
ref=float(np.max(abs(tr[1][:,5:]-tr[2][:,5:])));assert ref<1e-3
(H/'comparison.json').write_text(json.dumps({'runs':out,'work_trace_refinement_difference':ref},indent=2)+'\n')
print(json.dumps({'runs':out,'work_trace_refinement_difference':ref},indent=2))
