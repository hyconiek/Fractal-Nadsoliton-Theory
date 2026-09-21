#!/usr/bin/env python3
import os
import json, pathlib
from mpmath import iv, mp
mp.dps=80; iv.dps=60
ROOT=pathlib.Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup'))
c=json.load(open(ROOT/'certificates/R7P-026_equal_energy_event.json'))

def rat(s):
    if '/' in s:
        a,b=s.split('/'); return mp.mpf(a)/mp.mpf(b)
    return mp.mpf(s)
def pair(x): return [float(x.a),float(x.b)]
# d=(lambda3/6,lambda4/6,lambda5/6,lambda6/12)
d=[]
for k,den in [('3',6),('4',6),('5',6),('6',12)]:
    a,b=c['spectral_intervals'][k]
    d.append(iv.mpf([str(rat(a)/den),str(rat(b)/den)]))
s=[iv.mpf([a,b]) for a,b in c['root_box'][:4]]; g=iv.mpf(c['root_box'][4])
J=[s[i]*iv.sqrt(d[i]) for i in range(4)]
rt3=iv.sqrt(3); half=iv.mpf('0.5')
c3=[1,0,-1,0]*3
c4=[1,-half,-half]*4
c5=[1,-rt3/2,half,0,-half,rt3/2,-1,rt3/2,-half,0,half,-rt3/2]
alt=[1 if j%2==0 else -1 for j in range(12)]
F=[[iv.mpf(c3[j]),iv.mpf(c4[j]),iv.mpf(c5[j]),iv.mpf(alt[j])] for j in range(12)]
h=[sum(F[j][i]*J[i] for i in range(4)) for j in range(12)]; w=[iv.exp(x) for x in h]; Z=sum(w,iv.mpf(0)); p=[x/Z for x in w]
mu=[sum(p[j]*F[j][i] for j in range(12)) for i in range(4)]
C=[[sum(p[j]*(F[j][i]-mu[i])*(F[j][k]-mu[k]) for j in range(12)) for k in range(4)] for i in range(4)]
AJ=[[ (iv.mpf(1) if i==k else iv.mpf(0))-g*d[i]*C[i][k] for k in range(4)] for i in range(4)]
Inv=iv.matrix(AJ)**-1
S=[[Inv[i,j]*J[j] for j in range(4)] for i in range(4)]
row_abs=[]
for i in range(4):
    r=0.0
    for j in range(4): r+=max(abs(float(S[i][j].a)),abs(float(S[i][j].b)))
    row_abs.append(r)
# relative J sensitivity: d log J_i / d log d_j = S_ij/J_i
Rels=[[S[i][j]/J[i] for j in range(4)] for i in range(4)]
rel_rows=[]
for i in range(4):
    rel_rows.append(sum(max(abs(float(Rels[i][j].a)),abs(float(Rels[i][j].b))) for j in range(4)))
out={
 'task':'MP7-021 quantitative local moving-root sensitivity',
 'scientific_state':'PROVED_INTERVAL_ASSISTED_LOCAL_SENSITIVITY',
 'gain_box':c['root_box'][4],
 'J_box':[pair(x) for x in J],
 'dJ_dloglambda_matrix':[[pair(S[i][j]) for j in range(4)] for i in range(4)],
 'row_abs_sensitivity_upper':row_abs,
 'relative_dlogJ_dloglambda_row_sum_upper':rel_rows,
 'max_relative_inf_condition_upper':max(rel_rows),
 'interpretation':'For an infinitesimal vector alpha=d log lambda, ||d log J||_inf <= kappa ||alpha||_inf with displayed local interval upper kappa. This is derivative sensitivity at the certified root, not a finite perturbation radius by itself.'
}
print(json.dumps(out,indent=2))
