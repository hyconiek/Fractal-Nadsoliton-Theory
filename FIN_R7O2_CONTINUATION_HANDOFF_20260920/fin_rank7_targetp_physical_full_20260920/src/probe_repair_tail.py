from pathlib import Path
from fractions import Fraction as F
import json,sys,math
ROOT=Path(__file__).resolve().parents[1];R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src'));import physical_centered_moment as pc
d=json.load(open(ROOT/'results/repair_four_parents.json'))
def split(c,ax):
 c=[(F(a),F(b)) for a,b in c];lo,hi=c[ax];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[ax]=(lo,m);R[ax]=(m,hi);return L,R,m
rows=[]
for u in d['unresolved_terminal_leaves']:
 c=[(F(a),F(b)) for a,b in u['cell']]; axs=[]
 for ax in range(4):
  L,R,m=split(c,ax);zl,zr=pc.certify(L),pc.certify(R)
  axs.append({'axis':ax,'split':str(m),'pass_count':int(zl['ok'])+int(zr['ok']),'left_ok':zl['ok'],'right_ok':zr['ok'],'left_d3':zl['d3'],'right_d3':zr['d3'],'left_g':zl['gersh_lower'],'right_g':zr['gersh_lower']})
 rows.append({'original_index':u['original_index'],'path':u['path'],'cell':u['cell'],'axes':axs})
out={'task':'R7O2-probe-repair-tail','rows':rows};(ROOT/'results/repair_tail_axis_probe.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
