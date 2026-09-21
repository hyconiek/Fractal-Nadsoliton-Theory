from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time
ROOT=Path(__file__).resolve().parents[1]
R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src')); import physical_centered_moment as pc
MAX=int(sys.argv[1]) if len(sys.argv)>1 else 4
probe=json.load(open(ROOT/'results/depth6_axis_probe.json'))['rows']
def split(cell,axis):
 c=[(F(a),F(b)) for a,b in cell];lo,hi=c[axis];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[axis]=(lo,m);R[axis]=(m,hi);return L,R,m
def cert(cell): return pc.certify(cell)
def rec(cell,d,stats):
 z=cert(cell); stats['cert_calls']+=1
 if z['ok']: return 1,0,1
 if d==MAX: return 0,1,1
 candidates=[]
 for ax in (0,1,2):
  L,R,m=split(cell,ax); zl=cert(L);zr=cert(R);stats['cert_calls']+=2
  score=int(zl['ok'])+int(zr['ok'])
  # tie-break fewer aggregate widths then fixed axis
  width=(zl['entry_width_max'] if not zl['ok'] else 0)+(zr['entry_width_max'] if not zr['ok'] else 0)
  candidates.append((-score,width,ax,L,R,zl,zr))
 candidates.sort(key=lambda q:(q[0],q[1],q[2])); _,_,ax,L,R,zl,zr=candidates[0]
 stats['axis_counts'][ax]+=1
 s=u=n=0
 for C,Z in ((L,zl),(R,zr)):
  if Z['ok']: s+=1;n+=1
  else:
   a=rec(C,d+1,stats);s+=a[0];u+=a[1];n+=a[2]+1
 return s,u,n
st=time.time();rows=[]
for r in probe:
 cell=[(F(a),F(b)) for a,b in r['cell']];stats={'cert_calls':0,'axis_counts':[0,0,0]};s,u,n=rec(cell,0,stats);rows.append({'original_index':r['original_index'],'path':r['path'],'safe_terminal':s,'unresolved_terminal':u,'fully_closed':u==0,'nodes':n,**stats})
out={'max_extra_depth':MAX,'rows':rows,'closed':sum(x['fully_closed'] for x in rows),'unresolved_leaves':sum(x['unresolved_terminal'] for x in rows),'elapsed':time.time()-st,'axis_counts_total':[sum(r['axis_counts'][a] for r in rows) for a in range(3)]}
(ROOT/'results/dynamic_extension_probe.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'max_extra_depth':MAX,'closed':out['closed'],'unresolved_leaves':out['unresolved_leaves'],'elapsed':out['elapsed'],'axis_counts_total':out['axis_counts_total'],'failed':[(r['original_index'],r['path'],r['unresolved_terminal']) for r in rows if not r['fully_closed']]},indent=2))
