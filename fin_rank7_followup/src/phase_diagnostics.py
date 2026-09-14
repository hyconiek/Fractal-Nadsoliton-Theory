"""R7P-085 fixed-seed phase C2 diagnostics."""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from .phase_cumulants import full_phase_value_grad_hess,k4_phase_value_grad_hess

FIXTURE={'r3':0.1131879146,'r4':0.1698528641,'r5':0.2269339093,'z6':-0.3380663037}

def run(power=16,seed=85085):
    sob=qmc.Sobol(d=3,scramble=True,seed=seed)
    unit=sob.random_base2(power)
    phases=2*math.pi*unit
    dv=np.empty(len(phases)); dg=np.empty(len(phases)); dh=np.empty(len(phases))
    args=(FIXTURE['r3'],FIXTURE['r4'],FIXTURE['r5'],FIXTURE['z6'])
    for i,ph in enumerate(phases):
        vf,gf,Hf=full_phase_value_grad_hess(*args,ph)
        v4,g4,H4=k4_phase_value_grad_hess(*args,ph)
        dv[i]=abs(vf-v4); dg[i]=np.linalg.norm(gf-g4); dh[i]=np.linalg.norm(Hf-H4,2)
    summary={
      'seed':seed,'sobol_power':power,'samples':len(phases),'fixture':FIXTURE,
      'max_abs_value_diff':float(dv.max()),'argmax_value':int(dv.argmax()),
      'max_grad_diff':float(dg.max()),'argmax_grad':int(dg.argmax()),
      'max_hessian_op_diff':float(dh.max()),'argmax_hessian':int(dh.argmax()),
    }
    return summary,unit,phases,dv,dg,dh

def main(outdir):
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    summary,unit,phases,dv,dg,dh=run()
    np.savez_compressed(out/'R7P-085_sobol_samples.npz',unit=unit,phases=phases,value_diff=dv,grad_diff=dg,hessian_op_diff=dh)
    (out/'R7P-085_phase_diagnostics.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':
    import sys; main(sys.argv[1] if len(sys.argv)>1 else '.')
