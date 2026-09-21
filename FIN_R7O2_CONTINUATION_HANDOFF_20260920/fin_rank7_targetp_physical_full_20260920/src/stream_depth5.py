from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time,hashlib,multiprocessing as mp, os
ROOT=Path(__file__).resolve().parents[1];R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
base=json.load(open(ROOT/'inputs/target_p_residual_5432.json'));D=base['refined_failed'];root=[tuple(map(F,p)) for p in base['root_hull']]
BASEIDX=ROOT/'checkpoints/full_cover_chunks/index.json';V2DIR=ROOT/'checkpoints/stream_v2/parents';DIR=ROOT/'checkpoints/depth5_stream';PDIR=DIR/'parents';PDIR.mkdir(parents=True,exist_ok=True);IDX=DIR/'index.json'
PROCS=int(sys.argv[1]) if len(sys.argv)>1 else 5;BUDGET=float(sys.argv[2]) if len(sys.argv)>2 else 20;MAXN=int(sys.argv[3]) if len(sys.argv)>3 else 10000
axes=[2,0,1,2,0]
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
def pack(idx,c,d,path,z):
 return {'original_index':idx,'depth':d,'path':path,'cell':[[str(a),str(b)] for a,b in c],'ok':z['ok'],'reason':z['reason'],'basis_den':z['basis_den'],'basis_num':z['basis_num'],'rank_det':z['rank_det'],'d1':z['d1'],'d2':z['d2'],'d3':z['d3'],'gersh_lower':z['gersh_lower'],'entry_width_max':z['entry_width_max']}
def rec(c,d,path,idx,safe,unres):
 z=pc.certify(c)
 if z['ok']:safe.append(pack(idx,c,d,path,z));return
 if d==len(axes):unres.append(pack(idx,c,d,path,z));return
 L,R=split(c,axes[d]);rec(L,d+1,path+'L',idx,safe,unres);rec(R,d+1,path+'R',idx,safe,unres)
def work(idx):
 c=tuple((F(a),F(b)) for a,b in D[idx]['cell']);safe=[];unres=[];st=time.time();rec(c,0,'',idx,safe,unres)
 assert sum((V(tuple((F(a),F(b)) for a,b in x['cell'])) for x in safe+unres),F(0))==V(c)
 return {'task':'R7O2-depth5-parent','index':idx,'parent_cell':D[idx]['cell'],'parent_path':D[idx].get('path'),'parent_volume_fraction_root':str(V(c)/rootvol),'safe_terminal_leaves':safe,'unresolved_terminal_leaves':unres,'fully_closed':len(unres)==0,'elapsed_seconds':time.time()-st,'policy':{'axis_sequence':['t','r','s','t','r'],'max_depth':5}}
def scan_files():
 rows=[]
 for p in PDIR.glob('parent_*.json'):
  try:d=json.load(open(p));rows.append({'index':d['index'],'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'fully_closed':d['fully_closed'],'unresolved_terminal_count':len(d['unresolved_terminal_leaves'])})
  except:pass
 return rows
def write_idx(rows,done_base):
 d={'task':'R7O2-depth5-stream-index','base_done_count':len(done_base),'rows':rows,'processed_count':len(rows),'closed_count':sum(r['fully_closed'] for r in rows),'unresolved_parent_indices':[r['index'] for r in rows if not r['fully_closed']]}
 tmp=IDX.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');os.replace(tmp,IDX);return d
if __name__=='__main__':
 baseidx=json.load(open(BASEIDX));done_base=set(baseidx['processed_indices']);done_base|={int(p.stem.split('_')[1]) for p in V2DIR.glob('parent_*.json')}
 rows=scan_files();done={r['index'] for r in rows};sel=[i for i in range(len(D)) if i not in done_base and i not in done][:MAXN]
 st=time.time();new=0;pool=mp.Pool(PROCS,initializer=_init)
 try:
  for d in pool.imap_unordered(work,sel,chunksize=1):
   p=PDIR/f"parent_{d['index']:04d}.json";tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');os.replace(tmp,p);new+=1
   rows.append({'index':d['index'],'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'fully_closed':d['fully_closed'],'unresolved_terminal_count':len(d['unresolved_terminal_leaves'])})
   if new%20==0 or not d['fully_closed']:write_idx(rows,done_base)
   if time.time()-st>=BUDGET:break
 finally:pool.terminate();pool.join()
 idx=write_idx(rows,done_base)
 print(json.dumps({'new_this_call':new,'elapsed':time.time()-st,'processed':idx['processed_count'],'closed':idx['closed_count'],'unresolved_count':len(idx['unresolved_parent_indices']),'recent_unresolved':idx['unresolved_parent_indices'][-10:],'base_done_count':len(done_base),'remaining':len(D)-len(done_base)-idx['processed_count']},indent=2))
