from pathlib import Path
from fractions import Fraction as F
import json,sys,math
ROOT=Path(__file__).resolve().parents[1];R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src')); import physical_centered_moment as pc
rows=json.load(open(ROOT/'results/hard2_axis_probe.json'))['rows']
def split(c,ax):
 c=[(F(a),F(b)) for a,b in c];lo,hi=c[ax];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[ax]=(lo,m);R[ax]=(m,hi);return L,R,m
def V(c):
 z=F(1)
 for a,b in c:z*=b-a
 return z
outs=[]
for row in rows:
 c=[(F(a),F(b)) for a,b in row['cell']]; rootv=V(c); safe=[];steps=[]; final=None
 for depth in range(16):
  z0=pc.certify(c)
  if z0['ok']:
   final={'ok':True,'cell':[[str(a),str(b)] for a,b in c],'cert':{k:z0[k] for k in ['reason','d1','d2','d3','gersh_lower','basis_den','basis_num','rank_det']}};break
  cand=[]
  for ax in (0,1,2):
   L,R,m=split(c,ax); zl,zr=pc.certify(L),pc.certify(R)
   kids=[(L,zl,'L'),(R,zr,'R')]; passes=[k for k in kids if k[1]['ok']]; fails=[k for k in kids if not k[1]['ok']]
   # score: maximize immediate pass count, then best d3 lower of remaining failure, then gersh lower
   if fails:
    bestfail=max(fails,key=lambda k:(k[1]['d3'][0],k[1]['gersh_lower']))
    quality=(len(passes),bestfail[1]['d3'][0],bestfail[1]['gersh_lower'])
   else:
    bestfail=None;quality=(2,1e9,1e9)
   cand.append((quality,ax,m,kids,bestfail))
  cand.sort(key=lambda x:x[0],reverse=True);quality,ax,m,kids,bf=cand[0]
  safe_k=[k for k in kids if k[1]['ok']]
  for C,Z,side in safe_k:
   safe.append({'depth':depth+1,'axis':ax,'side':side,'cell':[[str(a),str(b)] for a,b in C],'reason':Z['reason'],'d3':Z['d3'],'gersh_lower':Z['gersh_lower'],'basis_den':Z['basis_den'],'basis_num':Z['basis_num'],'rank_det':Z['rank_det']})
  steps.append({'depth':depth+1,'axis':ax,'split':str(m),'pass_count':len(safe_k),'chosen_quality':quality,'children':[{'side':side,'ok':Z['ok'],'d3':Z['d3'],'g':Z['gersh_lower']} for C,Z,side in kids]})
  if bf is None:
   final={'ok':True,'cell':None,'cert':'both children safe'}; c=None; break
  c=bf[0]
 else:
  z=pc.certify(c); final={'ok':z['ok'],'cell':[[str(a),str(b)] for a,b in c],'cert':{k:z[k] for k in ['reason','d1','d2','d3','gersh_lower','basis_den','basis_num','rank_det']}}
 sv=sum((V([(F(a),F(b)) for a,b in s['cell']]) for s in safe),F(0)); fv=F(0) if c is None or final['ok'] else V(c)
 outs.append({'original_index':row['original_index'],'source_cell':row['cell'],'steps':steps,'safe_siblings':safe,'final':final,'covered_fraction':str(sv/rootv),'remaining_fraction':str(fv/rootv)})
res={'task':'R7O2-hard2-greedy-corner-chase','rows':outs}
(ROOT/'results/hard2_corner_chase.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps([{'index':x['original_index'],'steps':len(x['steps']),'final_ok':x['final']['ok'],'covered_fraction':x['covered_fraction'],'remaining_fraction':x['remaining_fraction'],'axes':[s['axis'] for s in x['steps']]} for x in outs],indent=2))
