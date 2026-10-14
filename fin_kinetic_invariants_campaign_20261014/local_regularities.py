#!/usr/bin/env python3
"""Tiny deterministic fold-line rate checks; diagnostic only, no trajectories."""
from pathlib import Path
import json, numpy as np
import sys
base=Path(__file__).resolve().parent
sys.path.insert(0,str(base/'source_context'))
import theory_core as tc
P=json.loads((base/'KI-002'/'DERIVED_SOFT_CLOCKS.json').read_text());f=tc.solve_fold();g=f['g'];p=f['p'];X=tc.X;v=f['v7'];r=(np.diag(p)-np.outer(p,p))@X@v;A=tc.A
rows=[]
for kin in ['HB','MET','BAR']:
    ell=np.array(P['models'][kin]['ell']);m=P['models'][kin]['mobility_m_nominal']
    for x in [.002,.0005,-.002,-.0005]:
        q=p+x*r;h=g*(A@q);d=h[None,:]-h[:,None]
        if kin=='HB':rate=np.tile(tc.softmax_vec(h),(12,1))
        elif kin=='MET':rate=np.exp(np.minimum(0.,d))/12
        else:rate=2/(12*(1+np.exp(-d)))
        C=q[:,None]*rate
        b=C.sum(axis=0)-C.sum(axis=1)
        projected=ell@b;pred=-m*f['b']*x*x/2
        rows.append({'kinetic':kin,'x':x,'true_projected_drift':float(projected),'quadratic_prediction':float(pred),'absolute_residual':float(abs(projected-pred)),'residual_over_absx3':float(abs(projected-pred)/abs(x)**3)})
(base/'KI-003'/'LOCAL_DRIFT_DIAGNOSTICS.json').write_text(json.dumps({'rows':rows,'scope':'binary64 local finite dimensional drift check; not a stochastic convergence theorem'},indent=2)+'\n')
print(json.dumps({'max_cubic_scaled_residual':max(t['residual_over_absx3'] for t in rows),'first_rows':rows[:3]},indent=2))
