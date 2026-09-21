from pathlib import Path
from fractions import Fraction as F
import json,sys,math,time,hashlib,multiprocessing as mp, os
ROOT=Path(__file__).resolve().parents[1]
R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920'); H=R7N/'inputs/FR223_20260916'
base=json.load(open(ROOT/'inputs/target_p_residual_5432.json')); D=base['refined_failed']; root=[tuple(map(F,p)) for p in base['root_hull']]
DIR=ROOT/'checkpoints/full_cover_chunks'; DIR.mkdir(parents=True,exist_ok=True); IDX=DIR/'index.json'
N=int(sys.argv[1]) if len(sys.argv)>1 else 100; PROCS=int(sys.argv[2]) if len(sys.argv)>2 else 4
axes=[2,0,1,2,0]
def V(cell):
    z=F(1)
    for a,b in cell:z*=F(b)-F(a)
    return z
rootvol=V(root)
def _init():
    global pc
    sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
    from verify_continuation import bounded_rationals
    bounded_rationals(9)
    sys.path.insert(0,str(ROOT/'src')); import physical_centered_moment as _pc; pc=_pc
def split(cell,axis):
    c=[(F(a),F(b)) for a,b in cell];lo,hi=c[axis];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
    if not lo<m<hi:m=(lo+hi)/2
    L=list(c);R=list(c);L[axis]=(lo,m);R[axis]=(m,hi);return L,R
def cert_rec(cell,depth,path,idx,safe,unres):
    z=pc.certify(cell)
    rec={'original_index':idx,'depth':depth,'path':path,'cell':[[str(a),str(b)] for a,b in cell],'ok':z['ok'],'reason':z['reason'],'basis_den':z['basis_den'],'basis_num':z['basis_num'],'rank_det':z['rank_det'],'d1':z['d1'],'d2':z['d2'],'d3':z['d3'],'gersh_lower':z['gersh_lower'],'entry_width_max':z['entry_width_max']}
    if z['ok']:safe.append(rec);return
    if depth==len(axes):unres.append(rec);return
    L,R=split(cell,axes[depth]);cert_rec(L,depth+1,path+'L',idx,safe,unres);cert_rec(R,depth+1,path+'R',idx,safe,unres)
def work(idx):
    cell=tuple((F(a),F(b)) for a,b in D[idx]['cell']);safe=[];unres=[];cert_rec(cell,0,'',idx,safe,unres)
    total=sum((V(tuple((F(a),F(b)) for a,b in x['cell'])) for x in safe+unres),F(0)); assert total==V(cell)
    sv=sum((V(tuple((F(a),F(b)) for a,b in x['cell'])) for x in safe),F(0))
    parent={'index':idx,'parent_cell':D[idx]['cell'],'parent_path':D[idx].get('path'),'parent_volume_fraction_root':str(V(cell)/rootvol),'safe_terminal_count':len(safe),'unresolved_terminal_count':len(unres),'safe_volume_fraction_of_parent':str(sv/V(cell)),'fully_closed':len(unres)==0}
    return parent,safe,unres
if __name__=='__main__':
    if IDX.exists(): idxd=json.load(open(IDX))
    else: idxd={'task':'R7O2-full-5432-physical-coupled-cover','source_count':len(D),'policy':{'backend':'QI outward rounded 1e-9','axis_sequence':['t','r','s','t','r'],'max_split_depth':5,'split':'rationalized geometric midpoint'},'chunks':[],'processed_indices':[],'total_processed':0,'total_closed':0,'unresolved_parent_indices':[]}
    done=set(idxd['processed_indices']); sel=[i for i in range(len(D)) if i not in done][:N]
    st=time.time()
    with mp.Pool(PROCS,initializer=_init) as pool: results=pool.map(work,sel)
    out={'task':'R7O2-full-cover-chunk','policy':idxd['policy'],'processes':PROCS,'processed_parents':[],'safe_terminal_leaves':[],'unresolved_terminal_leaves':[]}
    for parent,safe,unres in results: out['processed_parents'].append(parent);out['safe_terminal_leaves'].extend(safe);out['unresolved_terminal_leaves'].extend(unres)
    out['stats']={'processed_parent_count':len(out['processed_parents']),'fully_closed_parent_count':sum(x['fully_closed'] for x in out['processed_parents']),'safe_terminal_leaf_count':len(out['safe_terminal_leaves']),'unresolved_terminal_leaf_count':len(out['unresolved_terminal_leaves']),'elapsed_seconds':time.time()-st}
    no=len(idxd['chunks']); name=f'chunk_{no:04d}.json'; p=DIR/name; p.write_text(json.dumps(out,indent=2)+'\n'); h=hashlib.sha256(p.read_bytes()).hexdigest()
    idxd['chunks'].append({'file':name,'sha256':h,'processed_count':len(out['processed_parents']),'closed_count':sum(x['fully_closed'] for x in out['processed_parents']),'unresolved_terminal_count':len(out['unresolved_terminal_leaves'])})
    idxd['processed_indices'].extend(x['index'] for x in out['processed_parents']); idxd['total_processed']=len(idxd['processed_indices']); idxd['total_closed']+=sum(x['fully_closed'] for x in out['processed_parents']); idxd['unresolved_parent_indices'].extend(x['index'] for x in out['processed_parents'] if not x['fully_closed']); IDX.write_text(json.dumps(idxd,indent=2)+'\n')
    print(json.dumps({'chunk':name,'processed':len(out['processed_parents']),'closed':sum(x['fully_closed'] for x in out['processed_parents']),'unresolved_leaves':len(out['unresolved_terminal_leaves']),'unresolved_parents':[x['index'] for x in out['processed_parents'] if not x['fully_closed']],'elapsed':out['stats']['elapsed_seconds'],'total_processed':idxd['total_processed']},indent=2))
