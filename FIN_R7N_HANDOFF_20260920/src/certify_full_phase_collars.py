from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,numpy as np,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); IR=ROOT/'inputs/intake_review_20260919'; sys.path.insert(0,str(IR))
import scientific_rechecks as sr
CAT=ROOT/'inputs/FR223_20260916/certificates/R7P-092_full_phase_roots.json'
RADII=['0.002','0.0015','0.001','0.00075','0.0005','0.0003','0.0002','0.00015','0.0001','0.000075','0.00005']

def qbound(center,r):
    rr=F(r); c=[sr.I(F(str(x))) for x in center]
    box=[sr.iv.mpf([sr.I(F(str(x))-rr),sr.I(F(str(x))+rr)]) for x in center]
    _,H0=sr.phase_FH(c,'full'); _,HB=sr.phase_FH(box,'full')
    Hm=np.array([[sr.mid(x) for x in row] for row in H0],float); A=np.linalg.inv(Hm)
    Ar=[[F(repr(float(A[i,j]))) for j in range(3)] for i in range(3)]
    # exact nonsingularity
    det=(Ar[0][0]*(Ar[1][1]*Ar[2][2]-Ar[1][2]*Ar[2][1])-Ar[0][1]*(Ar[1][0]*Ar[2][2]-Ar[1][2]*Ar[2][0])+Ar[0][2]*(Ar[1][0]*Ar[2][1]-Ar[1][1]*Ar[2][0]))
    assert det!=0
    rows=[]
    for i in range(3):
        s=F(0)
        for j in range(3):
            e=sr.I(int(i==j))-sum((sr.I(Ar[i][k])*HB[k][j] for k in range(3)),sr.I(0))
            lo,hi=sr.bounds(e); s+=max(abs(lo),abs(hi))
        rows.append(s)
    return max(rows),Ar

def run():
    roots=json.load(open(CAT))['roots']; outrows=[]; st=time.time()
    for idx,rec in enumerate(roots):
        best=None; bestA=None; tests=[]
        for r in RADII:
            q,A=qbound(rec['phase'],r); tests.append({'radius':r,'q_inf':str(q),'pass':q<1})
            if q<1 and best is None: best=r; bestA=A
        assert best is not None
        outrows.append({'id':idx,'quartic_id':rec['quartic_id'],'center':[str(x) for x in rec['phase']],
                        'certified_radius':best,'tests':tests,
                        'preconditioner':[[str(x) for x in row] for row in bestA]})
        if (idx+1)%10==0: print('done',idx+1,flush=True)
    hist={}
    for x in outrows: hist[x['certified_radius']]=hist.get(x['certified_radius'],0)+1
    out={'task':'R7N-041','method':'full log-mgf uniform infinity-norm injectivity sup ||I-A H_full(phi)||_inf < 1; existence inherited from certified radius-1e-7 full root box',
         'count':60,'certified_count':60,'radius_histogram':hist,'min_radius':min(F(x['certified_radius']) for x in outrows),
         'max_radius':max(F(x['certified_radius']) for x in outrows),'elapsed_seconds':time.time()-st,'roots':outrows}
    (ROOT/'results/R7N-041_full_uniqueness_collars.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps({k:(str(v) if isinstance(v,F) else v) for k,v in out.items() if k!='roots'},indent=2))
if __name__=='__main__': run()
