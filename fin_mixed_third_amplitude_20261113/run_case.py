"""Third amplitude and joint space/time refinement; frozen solver dependency."""
from pathlib import Path
import sys,json,csv
import numpy as np
H=Path(__file__).resolve().parent
sys.path.insert(0,str(H.parent/'fin_mixed_scaling_phase_20261113'))
from campaign import run

def job(index):
    phase,dx,dt=[(0.,.1,.02),(np.pi/2,.1,.02),(0.,.05,.01)][index]
    name=f'phase{phase:.6f}_dx{dx}'
    if (H/(name+'.json')).exists():return
    r=run(1.5,dx,dt,.05,phase);samples=r.pop('samples')
    with (H/(name+'.csv')).open('w') as fp:
        w=csv.writer(fp);w.writerow(['t','I','J']);w.writerows([[v['t'],v['I'],v['J']] for v in samples])
    (H/(name+'.json')).write_text(json.dumps(r,indent=2)+'\n')
    print(name,r['relative_global_energy_drift'],flush=True)
if __name__=='__main__':job(int(sys.argv[1]))
