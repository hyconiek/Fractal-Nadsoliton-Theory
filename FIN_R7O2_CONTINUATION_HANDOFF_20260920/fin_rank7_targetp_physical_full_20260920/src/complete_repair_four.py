from pathlib import Path
from fractions import Fraction as F
import json,sys,math,hashlib
ROOT=Path(__file__).resolve().parents[1];R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src'));import physical_centered_moment as pc
src=json.load(open(ROOT/'results/repair_four_parents.json'))
def V(c):
 z=F(1)
 for a,b in c:z*=b-a
 return z
def split(c,ax):
 c=[(F(a),F(b)) for a,b in c];lo,hi=c[ax];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[ax]=(lo,m);R[ax]=(m,hi);return L,R,m
safe=list(src['safe_terminal_leaves']); new=[]
for u in src['unresolved_terminal_leaves']:
 c=[(F(a),F(b)) for a,b in u['cell']];L,R,m=split(c,0)
 for side,C in [('L',L),('R',R)]:
  z=pc.certify(C);assert z['ok']
  new.append({'original_index':u['original_index'],'depth':u['depth']+1,'path':u['path']+side,'stage':'final-r-repair','axis':0,'cell':[[str(a),str(b)] for a,b in C],'ok':True,'reason':z['reason'],'basis_den':z['basis_den'],'basis_num':z['basis_num'],'rank_det':z['rank_det'],'d1':z['d1'],'d2':z['d2'],'d3':z['d3'],'gersh_lower':z['gersh_lower'],'entry_width_max':z['entry_width_max']})
 # exact replacement
 assert V(L)+V(R)==V(c)
safe+=new
# exact coverage each parent using terminal leaves
base=json.load(open(ROOT/'inputs/target_p_residual_5432.json'))['refined_failed'];parents=[]
for idx in (332,338,357,363):
 pcells=[x for x in safe if x['original_index']==idx]; pv=V([(F(a),F(b)) for a,b in base[idx]['cell']]);sv=sum((V([(F(a),F(b)) for a,b in x['cell']]) for x in pcells),F(0));assert sv==pv
 parents.append({'index':idx,'safe_terminal_count':len(pcells),'coverage_exact':True})
out={'task':'R7O2-complete-repair-four-parents','source':'repair_four_parents.json','parents':parents,'safe_terminal_leaves':safe,'unresolved_terminal_leaves':[],'coverage_exact':True,'scientific_status':'INTERVAL_CERTIFIED_LOCAL_REPAIR'}
p=ROOT/'results/repair_four_parents_complete.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'parents':parents,'safe_leaves':len(safe),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
