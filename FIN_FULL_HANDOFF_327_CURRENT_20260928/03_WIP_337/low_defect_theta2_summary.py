import json, numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'02_RESEARCH_331_336'/'fin336'
meta=json.loads((ROOT/'DEFECT_RESPONSE_336.json').read_text())
tail6={int(r['N']):float(r['theta2_tail_Dgt6']) for r in meta['rows']}
out={'status':'WIP_NOT_A_COMPLETED_RESEARCH_TASK','theta':2.0,'rows':[]}
for N in range(7,11):
    z=np.load(ROOT/f'DEFECT_RESPONSE_N{N}_336.npz')
    st=z['states']; pi=z['pi']; D=N-st[:,0]
    # theta=2 => pi_theta proportional pi * exp(theta*n0)
    lw=np.log(pi)+2.0*st[:,0]
    lw-=lw.max(); w=np.exp(lw); w/=w.sum()
    mass_cond=float(w[D<=2].sum())
    q6=tail6[N]
    mass_full=(1-q6)*mass_cond
    out['rows'].append({
        'N':N,
        'Dle6_library_states':int(len(st)),
        'Dle2_states':int(np.sum(D<=2)),
        'mass_Dle2_conditional_on_Dle6':mass_cond,
        'certified_tail_Dgt6':q6,
        'certified_full_mass_Dle2':mass_full,
        'certified_full_tail_Dgt2':1.0-mass_full,
    })
print(json.dumps(out,indent=2))
(Path(__file__).with_name('LOW_DEFECT_THETA2_SUMMARY_337_WIP.json')).write_text(json.dumps(out,indent=2)+'\n')
