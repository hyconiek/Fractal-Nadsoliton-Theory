from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time,hashlib,multiprocessing as mp
ROOT=Path(__file__).resolve().parents[1]
R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
base=json.load(open(ROOT/'inputs/target_p_residual_5432.json'));D=base['refined_failed'];root=[tuple(map(F,p)) for p in base['root_hull']]
DIR=ROOT/'checkpoints/full_cover_chunks';IDX=DIR/'index.json'
N=int(sys.argv[1]) if len(sys.argv)>1 else 100;PROCS=int(sys.argv[2]) if len(sys.argv)>2 else 5
fixed=[2,0,1,2,0];MAX_EXTRA=7

def V(c):
 z=F(1)
 for a,b in c:z*=F(b)-F(a)
 return z
rootvol=V(root)

def _init():
 global pc
 sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
 from verify_continuation import bounded_rationals
 bounded_rationals(9)
 sys.path.insert(0,str(ROOT/'src'));import physical_centered_moment as _pc;pc=_pc

def split(c,ax):
 c=[(F(a),F(b)) for a,b in c];lo,hi=c[ax];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[ax]=(lo,m);R[ax]=(m,hi);return L,R

def pack(idx,c,depth,path,z,stage,axis=None):
 return {'original_index':idx,'depth':depth,'path':path,'stage':stage,'axis':axis,'cell':[[str(a),str(b)] for a,b in c],'ok':z['ok'],'reason':z['reason'],'basis_den':z['basis_den'],'basis_num':z['basis_num'],'rank_det':z['rank_det'],'d1':z['d1'],'d2':z['d2'],'d3':z['d3'],'gersh_lower':z['gersh_lower'],'entry_width_max':z['entry_width_max']}

def dynamic(c,extra,path,idx,safe,unres,depth0):
 z=pc.certify(c)
 if z['ok']:safe.append(pack(idx,c,depth0+extra,path,z,'dynamic'));return
 if extra==MAX_EXTRA:unres.append(pack(idx,c,depth0+extra,path,z,'dynamic'));return
 cand=[]
 for ax in (0,1,2):
  L,R=split(c,ax);zl,zr=pc.certify(L),pc.certify(R);kids=[(L,zl,'L'),(R,zr,'R')];passes=sum(k[1]['ok'] for k in kids);fails=[k for k in kids if not k[1]['ok']]
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
 ax=fixed[d];L,R=split(c,ax);fixedrec(L,d+1,path+'L',idx,safe,unres);fixedrec(R,d+1,path+'R',idx,safe,unres)

def work(idx):
 c=tuple((F(a),F(b)) for a,b in D[idx]['cell']);safe=[];unres=[];fixedrec(c,0,'',idx,safe,unres)
 assert sum((V(tuple((F(a),F(b)) for a,b in x['cell'])) for x in safe+unres),F(0))==V(c)
 sv=sum((V(tuple((F(a),F(b)) for a,b in x['cell'])) for x in safe),F(0))
 return {'index':idx,'parent_cell':D[idx]['cell'],'parent_path':D[idx].get('path'),'parent_volume_fraction_root':str(V(c)/rootvol),'safe_terminal_count':len(safe),'unresolved_terminal_count':len(unres),'safe_volume_fraction_of_parent':str(sv/V(c)),'fully_closed':len(unres)==0},safe,unres

if __name__=='__main__':
 idxd=json.load(open(IDX));done=set(idxd['processed_indices']);sel=[i for i in range(len(D)) if i not in done][:N]
 st=time.time()
 with mp.Pool(PROCS,initializer=_init) as pool:rs=pool.map(work,sel)
 policy={'version':'v2','backend':'QI outward rounded 1e-9','fixed_axis_sequence':['t','r','s','t','r'],'dynamic_extra_depth':MAX_EXTRA,'dynamic_axes':['r','s','t'],'dynamic_selection':'maximize immediate PASS; then failed-child d3.lo; then Gershgorin lower; fixed tie r<s<t'}
 out={'task':'R7O2-full-cover-chunk-v2','policy':policy,'processes':PROCS,'processed_parents':[],'safe_terminal_leaves':[],'unresolved_terminal_leaves':[]}
 for p,s,u in rs:out['processed_parents'].append(p);out['safe_terminal_leaves']+=s;out['unresolved_terminal_leaves']+=u
 out['stats']={'processed_parent_count':len(out['processed_parents']),'fully_closed_parent_count':sum(x['fully_closed'] for x in out['processed_parents']),'safe_terminal_leaf_count':len(out['safe_terminal_leaves']),'unresolved_terminal_leaf_count':len(out['unresolved_terminal_leaves']),'elapsed_seconds':time.time()-st}
 no=len(idxd['chunks']);name=f'chunk_{no:04d}_v2.json';p=DIR/name;p.write_text(json.dumps(out,indent=2)+'\n');h=hashlib.sha256(p.read_bytes()).hexdigest()
 idxd['chunks'].append({'file':name,'sha256':h,'policy_version':'v2','processed_count':len(out['processed_parents']),'closed_count':sum(x['fully_closed'] for x in out['processed_parents']),'unresolved_terminal_count':len(out['unresolved_terminal_leaves'])})
 idxd['processed_indices'].extend(x['index'] for x in out['processed_parents']);idxd['total_processed']=len(idxd['processed_indices']);idxd['total_closed']+=sum(x['fully_closed'] for x in out['processed_parents']);newu=[x['index'] for x in out['processed_parents'] if not x['fully_closed']];idxd['unresolved_parent_indices'].extend(newu);idxd['effective_closed_count']=idxd['total_closed']+idxd.get('repaired_closed_count',0);IDX.write_text(json.dumps(idxd,indent=2)+'\n')
 print(json.dumps({'chunk':name,'processed':len(out['processed_parents']),'closed':sum(x['fully_closed'] for x in out['processed_parents']),'unresolved_leaves':len(out['unresolved_terminal_leaves']),'unresolved_parents':newu,'elapsed':out['stats']['elapsed_seconds'],'total_processed':idxd['total_processed'],'effective_closed':idxd['effective_closed_count']},indent=2))
