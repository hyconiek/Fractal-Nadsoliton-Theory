from pathlib import Path
from fractions import Fraction as F
import sys,json,itertools,numpy as np
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); H=ROOT/'inputs/FR223_20260916'; sys.path[:0]=[str(ROOT/'src'),str(H/'src')]
import target_p_trace_cover as tc
from intervals import QI, sqrt_interval
import off_face
TAU=F(67,250)
# rigorous 4D feature intervals for representative states
S=[sqrt_interval(tc.L[3]/6,40),sqrt_interval(tc.L[4]/6,40),sqrt_interval(tc.L[5]/6,40),sqrt_interval(tc.L[6]/12,40)]
OBS=[]
for idx in tc.REPS:
    OBS.append([tc.CV[k][idx]*S[k] for k in range(4)])

def pair_prob_bounds(wb,i,j):
    Dhi=sum(b for a,b in wb)
    lo=wb[i][0]*wb[j][0]/(Dhi*Dhi)
    hi=tc.pair_prob_product_upper(wb,i,j)
    return QI(lo,hi)

def rank_minor(Bn):
    # Bn integer 4x3
    from itertools import combinations
    best=(None,0)
    for rows in combinations(range(4),3):
        A=np.array([Bn[r] for r in rows],dtype=object)
        det=(A[0,0]*(A[1,1]*A[2,2]-A[1,2]*A[2,1])
             -A[0,1]*(A[1,0]*A[2,2]-A[1,2]*A[2,0])
             +A[0,2]*(A[1,0]*A[2,1]-A[1,1]*A[2,0]))
        if abs(det)>abs(best[1]): best=(rows,int(det))
    return best

def center_basis(cell,den=100000):
    x=np.array([float((a+b)/2) for a,b in cell])
    p=off_face.p_from_aligned_compact(*x)
    _,_,_,_,C4=off_face.constants(); mu=p@C4; Y=C4-mu; M=Y.T@(p[:,None]*Y)
    w,V=np.linalg.eigh(M)
    B=V[:,:3]
    Bn=np.rint(B*den).astype(np.int64)
    rows,det=rank_minor(Bn.tolist())
    return Bn,den,w,rows,det

def certify(cell,den=100000):
    Bn,den,eigs,rows,det=center_basis(cell,den)
    if det==0:return {'ok':False,'reason':'RANK'}
    B=[[F(int(Bn[i,j]),den) for j in range(3)] for i in range(4)]
    # exact a B^T B
    K=[[QI(TAU*sum(B[r][i]*B[r][j] for r in range(4))) for j in range(3)] for i in range(3)]
    wb=tc.weight_bounds(cell)
    for i,j in itertools.combinations(range(7),2):
        q=pair_prob_bounds(wb,i,j)
        dz=[]
        for k in range(3):
            z=QI(0)
            for r in range(4):z += B[r][k]*(OBS[i][r]-OBS[j][r])
            dz.append(z)
        for a in range(3):
            for b in range(a,3):
                K[a][b]=K[a][b]-q*dz[a]*dz[b]
                if a!=b:K[b][a]=K[a][b]
    d1=K[0][0]
    d2=K[0][0]*K[1][1]-K[0][1]*K[0][1]
    d3=(K[0][0]*K[1][1]*K[2][2]+2*K[0][1]*K[0][2]*K[1][2]
        -K[0][0]*K[1][2]*K[1][2]-K[1][1]*K[0][2]*K[0][2]-K[2][2]*K[0][1]*K[0][1])
    sylv=d1.lo>0 and d2.lo>0 and d3.lo>0
    g=[]
    for a in range(3):
        off=sum(max(abs(K[a][b].lo),abs(K[a][b].hi)) for b in range(3) if b!=a)
        g.append(K[a][a].lo-off)
    gersh=min(g)>0
    ok=sylv or gersh
    return {'ok':ok,'reason':'SYLVESTER_PD' if sylv else ('GERSHGORIN_PD' if gersh else 'FAILED'),'rank_rows':rows,'rank_det':det,
            'basis_den':den,'basis_num':Bn.tolist(),'center_eigs':eigs.tolist(),
            'gersh_lower':float(min(g)),'d1':[float(d1.lo),float(d1.hi)],'d2':[float(d2.lo),float(d2.hi)],'d3':[float(d3.lo),float(d3.hi)]}

if __name__=='__main__':
 d=json.load(open(ROOT/'checkpoints/R7N-020_trace_e2_hybrid_v2.json'))
 arr=d['e2_failed']
 n=int(sys.argv[1]) if len(sys.argv)>1 else 50
 out=[]
 for x in arr[:n]:
    cell=tuple(tuple(map(F,p)) for p in x['cell']); r=certify(cell); out.append({'path':x['path'],**r})
 print(json.dumps({'n':n,'pass':sum(x['ok'] for x in out),'results':out},indent=2))
