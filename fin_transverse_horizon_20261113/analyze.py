from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent
rows=json.loads((H/'results.json').read_text());assert len(rows)==3
tr=[]
for r in rows:
 a=np.loadtxt(H/(r['name']+'.csv'),delimiter=',',skiprows=1);assert a.shape==(401,5) and np.isfinite(a).all() and np.all(a[:,1]>0);tr.append(a)
assert rows[0]['static_max_deviation']<1e-8
c,f=tr[1:];err=float(np.max(abs(c[:,1:]-f[:,1:])));assert err<1e-3
old=H.parent/'fin_transverse_phase_gap_20261113'/('driven_w1.25_phase1.570796_dx0.1_dt0.02.csv')
o=np.loadtxt(old,delimiter=',',skiprows=1);overlap=float(np.max(abs(o[:,1:]-c[:201,1:])));assert overlap<1e-6
snap=[]
for t in [0,50,100,150,200]:
 a=f[round(t*2)];snap.append(dict(t=t,min_gain=a[1],max_gain=a[2],fixed_position_gain=a[3],fixed_velocity_gain=a[4]))
out={'samples':snap,'max_gain':float(max(f[:,2])),'time_of_sampled_max':float(f[np.argmax(f[:,2]),0]),'time_refinement_max_absolute_gain_difference':err,'domain_overlap_max_absolute_gain_difference':overlap,'static_control_drift':rows[0]['static_max_deviation']}
(H/'comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
