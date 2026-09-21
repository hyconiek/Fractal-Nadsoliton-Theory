from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time,hashlib,multiprocessing as mp, os
ROOT=Path(__file__).resolve().parents[1]
R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
base=json.load(open(ROOT/'inputs/target_p_residual_5432.json'));D=base['refined_failed'];root=[tuple(map(F,p)) for p in base['root_hull']]
OLDIDX=ROOT/'checkpoints/full_cover_chunks/index.json'
DIR=ROOT/'checkpoints/stream_v2'; PDIR=DIR/'parents'; PDIR.mkdir(parents=True,exist_ok=True); IDX=DIR/'index.json'
PROCS=int(sys.argv[1]) if len(sys.argv)>1 else 5; BUDGET=float(sys.argv[2]) if len(sys.argv)>2 else 35.0; MAXN=int(sys.argv[3]) if len(sys.argv)>3 else 10000
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
 c=tuple((F(a),F(b)) for a,b in D[idx]['cell']);safe=[];unres=[];st=time.time();fixedrec(c,0,'',idx,safe,unres)
 assert sum((V(tuple((F(a),F(b)) for a,b in x['cell'])) for x in safe+unres),F(0))==V(c)
 return {'task':'R7O2-stream-parent-v2','index':idx,'parent_cell':D[idx]['cell'],'parent_path':D[idx].get('path'),'parent_volume_fraction_root':str(V(c)/rootvol),'safe_terminal_leaves':safe,'unresolved_terminal_leaves':unres,'fully_closed':len(unres)==0,'elapsed_seconds':time.time()-st,'policy':{'fixed_axes':['t','r','s','t','r'],'dynamic_extra_depth':MAX_EXTRA,'dynamic_axes':['r','s','t']}}
def load_index():
 old=json.load(open(OLDIDX)); old_done=set(old['processed_indices'])
 files=sorted(PDIR.glob('parent_*.json')); rows=[]
 for p in files:
  try:d=json.load(open(p)); rows.append({'index':d['index'],'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'fully_closed':d['fully_closed'],'unresolved_terminal_count':len(d['unresolved_terminal_leaves'])})
  except Exception: pass
 return old,old_done,rows
def write_index(old,rows):
 d={'task':'R7O2-stream-v2-index','base_chunk_index':'../full_cover_chunks/index.json','base_processed_count':old['total_processed'],'base_effective_closed_count':old.get('effective_closed_count',old['total_closed']),'stream_rows':rows,'stream_processed_count':len(rows),'stream_closed_count':sum(r['fully_closed'] for r in rows),'stream_unresolved_parent_indices':[r['index'] for r in rows if not r['fully_closed']],'combined_processed_count':old['total_processed']+len(rows),'combined_effective_closed_count':old.get('effective_closed_count',old['total_closed'])+sum(r['fully_closed'] for r in rows)}
 tmp=IDX.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');os.replace(tmp,IDX);return d
if __name__=='__main__':
 old,old_done,rows=load_index();stream_done={r['index'] for r in rows};sel=[i for i in range(len(D)) if i not in old_done and i not in stream_done][:MAXN]
 st=time.time();new=0
 pool=mp.Pool(PROCS,initializer=_init)
 try:
  it=pool.imap_unordered(work,sel,chunksize=1)
  for d in it:
   p=PDIR/f"parent_{d['index']:04d}.json";tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');os.replace(tmp,p);new+=1
   rows.append({'index':d['index'],'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'fully_closed':d['fully_closed'],'unresolved_terminal_count':len(d['unresolved_terminal_leaves'])})
   if new%5==0 or not d['fully_closed']: write_index(old,rows)
   if time.time()-st>=BUDGET or not d['fully_closed']: break
 finally:
  pool.terminate();pool.join()
 idx=write_index(old,rows)
 print(json.dumps({'new_this_call':new,'elapsed':time.time()-st,'stream_processed':idx['stream_processed_count'],'stream_closed':idx['stream_closed_count'],'stream_unresolved':idx['stream_unresolved_parent_indices'][-10:],'combined_processed':idx['combined_processed_count'],'combined_effective_closed':idx['combined_effective_closed_count'],'remaining':len(D)-idx['combined_processed_count']},indent=2))
