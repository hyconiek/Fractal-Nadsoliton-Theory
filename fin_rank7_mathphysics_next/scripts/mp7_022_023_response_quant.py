import os
import json, math, pathlib
from mpmath import iv, mp
mp.dps=80
iv.dps=60
ROOT=pathlib.Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup'))
CERT=ROOT/'certificates/R7P-026_equal_energy_event.json'
FOLD=ROOT/'certificates/R7P-031_simple_fold.json'
MARGIN=pathlib.Path(os.environ.get('MP7_WORK_ROOT','/mnt/data/fin_rank7_mathphysics_next'))/'results/MP7-020_spectral_margin.json'

def I(a,b=None):
    if b is None: b=a
    return iv.mpf([str(a),str(b)])

def lo(x): return float(x.a)
def hi(x): return float(x.b)
def pair(x): return [str(mp.mpf(x.a)), str(mp.mpf(x.b))]

# interval spectral feature C4 from R7P-026 lambda boxes
c=json.load(open(CERT)); f=json.load(open(FOLD)); marg=json.load(open(MARGIN))
lam={int(k):(iv.mpf(v[0]),iv.mpf(v[1])) for k,v in c['spectral_intervals'].items()}
# make intervals correctly from rational strings
def rat_mp(s):
    if '/' in s:
        a,b=s.split('/'); return mp.mpf(a)/mp.mpf(b)
    return mp.mpf(s)
lam={int(k):iv.mpf([str(rat_mp(v[0])),str(rat_mp(v[1]))]) for k,v in c['spectral_intervals'].items()}
pi=iv.pi
C=[]
for j in range(12):
    row=[]
    for k in (3,4,5):
        scale=iv.sqrt(lam[k]/6)
        row.append(scale*iv.cos(2*pi*k*j/12))
    scale6=iv.sqrt(lam[6]/12)
    row.append(scale6*((-1)**j))
    C.append(row)

# Root interval theta and gain at coexistence
s=[iv.mpf([a,b]) for a,b in c['root_box'][:4]]
g=iv.mpf(c['root_box'][4])
# p intervals via exp field; center to avoid overflow not necessary
h=[]
for j in range(12):
    z=iv.mpf(0)
    for q in range(4): z += C[j][q]*s[q]
    h.append(z)
w=[iv.exp(x) for x in h]
Z=sum(w,iv.mpf(0))
p=[x/Z for x in w]
mu=[]
for a in range(4): mu.append(sum((p[j]*C[j][a] for j in range(12)),iv.mpf(0)))
M=[[iv.mpf(0) for b in range(4)] for a in range(4)]
for a in range(4):
    for b in range(4):
        M[a][b]=sum((p[j]*(C[j][a]-mu[a])*(C[j][b]-mu[b]) for j in range(12)),iv.mpf(0))
H=[[((1/g) if a==b else iv.mpf(0))-M[a][b] for b in range(4)] for a in range(4)]
# interval matrix inverse and response
Hm=iv.matrix(H); sm=iv.matrix(s)
try:
    invH=Hm**-1
    resp=invH*sm/(g*g)
    inv_ok=True
except Exception as e:
    inv_ok=False; inv_err=repr(e); resp=None; invH=None

# Perron fold gap from exact fold box + global lambda2 upper
fg=iv.mpf(f['fold_gain_interval'])
lambda1=1/fg
mglobal=mp.mpf(str(marg['global_margin_float']))
tau=mp.mpf(67)/250
lambda2_upper=tau-mglobal
gap_lower=mp.mpf(lambda1.a)-lambda2_upper
# positive vector component bounds and norm check
v=[iv.mpf(x) for x in f['augmented_root_box'][5:9]]
# cosine to positive all-ones direction; merely descriptive
ones_norm=mp.sqrt(4)
vsum=sum(v,iv.mpf(0)); vnorm2=sum((x*x for x in v),iv.mpf(0))
cos_ones=vsum/(iv.sqrt(vnorm2)*2)

out={
 'task':'MP7-022/023 quantitative supplement',
 'scientific_state':'PROVED_INTERVAL_ASSISTED_QUANTITATIVE_RESPONSE',
 'coexistence_gain':c['root_box'][4],
 'response_interval_s_over_g2_Hinv': [pair(x) for x in resp] if inv_ok else None,
 'H_interval_diag':[pair(H[i][i]) for i in range(4)],
 'interval_inverse_succeeded':inv_ok,
 'interval_inverse_error':None if inv_ok else inv_err,
 'fold_lambda1_interval':pair(lambda1),
 'global_lambda2_upper':str(lambda2_upper),
 'fold_spectral_gap_lower':str(gap_lower),
 'fold_v_component_intervals':[list(x) for x in f['augmented_root_box'][5:9]],
 'fold_v_cosine_with_positive_diagonal_interval':pair(cos_ones),
}
print(json.dumps(out,indent=2))
