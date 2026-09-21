from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time
ROOT=Path(__file__).resolve().parents[1]
R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src')); import physical_centered_moment as pc
seq=[int(x) for x in sys.argv[1].split(',')]
probe=json.load(open(ROOT/'results/depth6_axis_probe.json'))['rows']
def split(cell,axis):
 c=[(F(a),F(b)) for a,b in cell];lo,hi=c[axis];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[axis]=(lo,m);R[axis]=(m,hi);return L,R
def rec(cell,d):
 z=pc.certify(cell)
 if z['ok']: return 1,0,1
 if d==len(seq): return 0,1,1
 L,R=split(cell,seq[d]);a=rec(L,d+1);b=rec(R,d+1);return a[0]+b[0],a[1]+b[1],a[2]+b[2]+1
st=time.time();rows=[]
for r in probe:
 cell=[(F(a),F(b)) for a,b in r['cell']];s,u,n=rec(cell,0);rows.append({'original_index':r['original_index'],'path':r['path'],'safe_terminal':s,'unresolved_terminal':u,'fully_closed':u==0,'nodes':n})
out={'sequence':seq,'rows':rows,'closed':sum(x['fully_closed'] for x in rows),'unresolved_leaves':sum(x['unresolved_terminal'] for x in rows),'elapsed':time.time()-st}
print(json.dumps(out,indent=2))
