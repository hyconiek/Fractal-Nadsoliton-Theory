from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time
ROOT=Path(__file__).resolve().parents[1];R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src')); import physical_centered_moment as pc
seq=[int(x) for x in sys.argv[1].split(',')]
rows=json.load(open(ROOT/'results/depth6_axis_probe.json'))['rows']
rows=[r for r in rows if r['original_index'] in (332,357) and r['path'] in ('RRRLR','RRRRR')]
def split(c,ax):
 c=[(F(a),F(b)) for a,b in c];lo,hi=c[ax];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[ax]=(lo,m);R[ax]=(m,hi);return L,R
def rec(c,d):
 z=pc.certify(c)
 if z['ok']:return 1,0
 if d==len(seq):return 0,1
 L,R=split(c,seq[d]);a=rec(L,d+1);b=rec(R,d+1);return a[0]+b[0],a[1]+b[1]
st=time.time();out=[]
for r in rows:
 s,u=rec([(F(a),F(b)) for a,b in r['cell']],0);out.append((r['original_index'],r['path'],s,u))
print(json.dumps({'seq':seq,'rows':out,'closed':sum(u==0 for _,_,s,u in out),'unresolved':sum(u for *_,u in out),'elapsed':time.time()-st},indent=2))
