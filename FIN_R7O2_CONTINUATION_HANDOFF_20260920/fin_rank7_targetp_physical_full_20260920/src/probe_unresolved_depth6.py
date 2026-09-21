from pathlib import Path
from fractions import Fraction as F
import json,sys,math,numpy as np
ROOT=Path(__file__).resolve().parents[1]
R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920'); H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(R7N/'src'))
import physical_centered_moment as pc
import off_face, compression_interval_probe as old
TAU=67/250
# collect unresolved from all chunks
DIR=ROOT/'checkpoints/full_cover_chunks'
idx=json.load(open(DIR/'index.json'))
U=[]
for c in idx['chunks']:
    d=json.load(open(DIR/c['file']))
    U.extend(d['unresolved_terminal_leaves'])
# numeric center lambda2
obs=np.array([[float((old.OBS[i][r].lo+old.OBS[i][r].hi)/2) for r in range(4)] for i in range(7)])
def lam2(cell):
    x=np.array([float((F(a)+F(b))/2) for a,b in cell])
    p=off_face.p_from_aligned_compact(*x)
    pa=np.array([p[0],p[4]+p[8],p[6],p[2]+p[10],p[3]+p[9],p[5]+p[7],p[1]+p[11]])
    mu=pa@obs; M=(obs-mu).T@((obs-mu)*pa[:,None]); e=np.linalg.eigvalsh(M); return e.tolist(), float(TAU-e[-2])
def split(cell,axis):
    c=[(F(a),F(b)) for a,b in cell];lo,hi=c[axis];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
    if not lo<m<hi:m=(lo+hi)/2
    L=list(c);R=list(c);L[axis]=(lo,m);R[axis]=(m,hi);return L,R,m
rows=[]
for j,u in enumerate(U):
    cell=[(F(a),F(b)) for a,b in u['cell']]
    eig,gap=lam2(cell); axes=[]
    for ax in range(4):
        L,R,m=split(cell,ax); zl=pc.certify(L); zr=pc.certify(R)
        axes.append({'axis':ax,'split':str(m),'left_ok':zl['ok'],'right_ok':zr['ok'],'pass_count':int(zl['ok'])+int(zr['ok']),'left_d3':zl['d3'],'right_d3':zr['d3'],'left_g':zl['gersh_lower'],'right_g':zr['gersh_lower']})
    rows.append({'id':j,'original_index':u['original_index'],'path':u['path'],'cell':u['cell'],'center_eigenvalues':eig,'tau_gap_center':gap,'axes':axes})
out={'task':'R7O2-probe-depth6-unresolved','count':len(rows),'rows':rows}
(ROOT/'results/depth6_axis_probe.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'count':len(rows),'by_parent':{str(i):sum(r['original_index']==i for r in rows) for i in sorted({r['original_index'] for r in rows})},'axis_full_pass_counts':{str(ax):sum(r['axes'][ax]['pass_count']==2 for r in rows) for ax in range(4)},'min_center_gap':min(r['tau_gap_center'] for r in rows)},indent=2))
