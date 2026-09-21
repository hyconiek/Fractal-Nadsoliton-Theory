from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time,multiprocessing as mp,hashlib
ROOT=Path(__file__).resolve().parents[1];R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
base=json.load(open(ROOT/'inputs/target_p_residual_5432.json'));D=base['refined_failed'];root=[tuple(map(F,p)) for p in base['root_hull']]
TARGET=[332,338,357,363];fixed=[2,0,1,2,0];MAX_EXTRA=6
def V(c):
 z=F(1)
 for a,b in c:z*=F(b)-F(a)
 return z
def _init():
 global pc
 sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
 from verify_continuation import bounded_rationals
 bounded_rationals(9)
 sys.path.insert(0,str(ROOT/'src'));import physical_centered_moment as _pc;pc=_pc
def split(c,ax):
 c=[(F(a),F(b)) for a,b in c];lo,hi=c[ax];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[ax]=(lo,m);R[ax]=(m,hi);return L,R,m
def pack(idx,c,depth,path,z,stage,axis=None):
 return {'original_index':idx,'depth':depth,'path':path,'stage':stage,'axis':axis,'cell':[[str(a),str(b)] for a,b in c],'ok':z['ok'],'reason':z['reason'],'basis_den':z['basis_den'],'basis_num':z['basis_num'],'rank_det':z['rank_det'],'d1':z['d1'],'d2':z['d2'],'d3':z['d3'],'gersh_lower':z['gersh_lower'],'entry_width_max':z['entry_width_max']}
def dynamic(c,extra,path,idx,safe,unres,depth0):
 z=pc.certify(c)
 if z['ok']:safe.append(pack(idx,c,depth0+extra,path,z,'dynamic'));return
 if extra==MAX_EXTRA:unres.append(pack(idx,c,depth0+extra,path,z,'dynamic'));return
 cand=[]
 for ax in (0,1,2):
  L,R,m=split(c,ax);zl,zr=pc.certify(L),pc.certify(R);kids=[(L,zl,'L'),(R,zr,'R')];passes=sum(k[1]['ok'] for k in kids);fails=[k for k in kids if not k[1]['ok']]
  if fails:
   best=max(fails,key=lambda k:(k[1]['d3'][0],k[1]['gersh_lower']));quality=(passes,best[1]['d3'][0],best[1]['gersh_lower'],-ax)
  else:quality=(2,1e9,1e9,-ax)
  cand.append((quality,ax,kids))
 cand.sort(key=lambda x:x[0],reverse=True);_,ax,kids=cand[0]
 for C,Z,side in kids:
  if Z['ok']:safe.append(pack(idx,C,depth0+extra+1,path+side,Z,'dynamic',ax))
  else:dynamic(C,extra+1,path+side,idx,safe,unres,depth0)
def fixedrec(c,d,path,idx,safe,unres):
 z=pc.certify(c)
 if z['ok']:safe.append(pack(idx,c,d,path,z,'fixed'));return
 if d==len(fixed):dynamic(c,0,path,idx,safe,unres,d);return
 ax=fixed[d];L,R,m=split(c,ax);fixedrec(L,d+1,path+'L',idx,safe,unres);fixedrec(R,d+1,path+'R',idx,safe,unres)
def work(idx):
 c=tuple((F(a),F(b)) for a,b in D[idx]['cell']);safe=[];unres=[];fixedrec(c,0,'',idx,safe,unres)
 assert sum((V(tuple((F(a),F(b)) for a,b in x['cell'])) for x in safe+unres),F(0))==V(c)
 return {'index':idx,'parent_cell':D[idx]['cell'],'safe_terminal_count':len(safe),'unresolved_terminal_count':len(unres),'fully_closed':len(unres)==0},safe,unres
if __name__=='__main__':
 st=time.time()
 with mp.Pool(4,initializer=_init) as pool:rs=pool.map(work,TARGET)
 out={'task':'R7O2-repair-four-depth5-bottlenecks','policy':{'fixed_axes':['t','r','s','t','r'],'dynamic_extra_depth':MAX_EXTRA,'dynamic_axes':['r','s','t'],'selection':'maximize immediate PASS count, then failed-child d3.lo, then Gershgorin lower, tie r<s<t'},'parents':[],'safe_terminal_leaves':[],'unresolved_terminal_leaves':[]}
 for p,s,u in rs:out['parents'].append(p);out['safe_terminal_leaves']+=s;out['unresolved_terminal_leaves']+=u
 out['stats']={'parents':len(out['parents']),'closed':sum(x['fully_closed'] for x in out['parents']),'safe_leaves':len(out['safe_terminal_leaves']),'unresolved_leaves':len(out['unresolved_terminal_leaves']),'elapsed':time.time()-st}
 p=ROOT/'results/repair_four_parents.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['stats'],indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
