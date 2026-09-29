import json, time, numpy as np
from fin_core import *
A7,L=A7_default(); Ap,lbar=projector_countermodel(); G=G_FROZEN
ge_fin=3.7183448981203777; ge_p=g_eq(Ap,lo=2.0,hi=8.0,n=40); g_p=G*ge_p/ge_fin
out={}
t=time.time()
rf,Rf,asym,_=slow_modes(5,G,A7,"barker"); rp,Rp,_,_=slow_modes(5,g_p,Ap,"barker")
out["Rk_N5_barker"]={"R_FIN":Rf,"R_proj":Rp,"rel_diff_pct":{k:100*(Rf[k]-Rp[k])/Rp[k] for k in Rf if k!=4},"DB_asym":asym}
print("N5 barker FINvsProj%:",{k:round(v,1) for k,v in out["Rk_N5_barker"]["rel_diff_pct"].items()},f"[{time.time()-t:.0f}s]",flush=True)
out["static_S"]={}
for N in (6,8):
    for g in (2.0,3.0):
        f=static_S(N,g,A7); p=static_S(N,g,Ap)
        out["static_S"][f"N{N}|g{g}"]={"FIN":f,"proj_same_g":p,"gauss_FIN":{k:(1/(1-g*float(L[k])/12) if k>=3 else 1.0) for k in range(1,7)}}
        print(f"static N={N} g={g}: FIN "+" ".join(f"k{k}:{f[k]:.3f}" for k in range(1,7))+" | proj "+" ".join(f"k{k}:{p[k]:.3f}" for k in range(1,7)),flush=True)
json.dump(out,open("seed_tail_results.json","w"),indent=1,default=float)
