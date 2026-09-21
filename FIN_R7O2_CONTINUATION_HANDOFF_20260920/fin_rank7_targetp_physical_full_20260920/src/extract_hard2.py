from pathlib import Path
from fractions import Fraction as F
import json,sys,math,numpy as np
ROOT=Path(__file__).resolve().parents[1];R7N=Path('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920');H=R7N/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(R7N/'inputs/intake_review_20260919'))
from verify_continuation import bounded_rationals
bounded_rationals(9)
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(R7N/'src'))
import physical_centered_moment as pc, off_face, compression_interval_probe as old
seq=[1,0,1,0,1]
rows=json.load(open(ROOT/'results/depth6_axis_probe.json'))['rows']; rows=[r for r in rows if r['original_index'] in (332,357) and r['path']=='RRRRR']
def split(c,ax):
 c=[(F(a),F(b)) for a,b in c];lo,hi=c[ax];m=F(format(math.sqrt(float(lo)*float(hi)),'.16g'))
 if not lo<m<hi:m=(lo+hi)/2
 L=list(c);R=list(c);L[ax]=(lo,m);R[ax]=(m,hi);return L,R,m
def rec(c,d,U):
 z=pc.certify(c)
 if z['ok']:return
 if d==len(seq):U.append(c);return
 L,R,m=split(c,seq[d]);rec(L,d+1,U);rec(R,d+1,U)
obs=np.array([[float((old.OBS[i][r].lo+old.OBS[i][r].hi)/2) for r in range(4)] for i in range(7)])
def eigs(c):
 x=np.array([float((a+b)/2) for a,b in c]);p=off_face.p_from_aligned_compact(*x);pa=np.array([p[0],p[4]+p[8],p[6],p[2]+p[10],p[3]+p[9],p[5]+p[7],p[1]+p[11]]);mu=pa@obs;M=(obs-mu).T@((obs-mu)*pa[:,None]);return np.linalg.eigvalsh(M)
out=[]
for r in rows:
 U=[];rec([(F(a),F(b)) for a,b in r['cell']],0,U);assert len(U)==1
 c=U[0]; e=eigs(c); axs=[]
 for ax in range(4):
  L,R,m=split(c,ax);zl=pc.certify(L);zr=pc.certify(R);axs.append({'axis':ax,'split':str(m),'pass_count':int(zl['ok'])+int(zr['ok']),'left_ok':zl['ok'],'right_ok':zr['ok'],'left_d3':zl['d3'],'right_d3':zr['d3'],'left_g':zl['gersh_lower'],'right_g':zr['gersh_lower']})
 out.append({'original_index':r['original_index'],'cell':[[str(a),str(b)] for a,b in c],'center_eigs':e.tolist(),'tau_gap':float(F(67,250))-float(e[-2]),'axes':axs})
(ROOT/'results/hard2_axis_probe.json').write_text(json.dumps({'rows':out},indent=2)+'\n')
print(json.dumps(out,indent=2))
