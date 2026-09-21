from pathlib import Path
from fractions import Fraction as F
import sys,json,math,numpy as np,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR))
import scientific_rechecks as sr
CAT=ROOT/'inputs/FR223_20260916/results/R7P-089_quartic_roots.json'
# exact K4 phase coefficients
def terms():
 rt3=sr.iv.sqrt(sr.I(3));r3,r4,r5,z6=map(sr.I,['0.1131879146','0.1698528641','0.2269339093','-0.3380663037']);a,b,c,d=[x/(2*rt3) for x in [r3,r4,r5,z6]]
 cub=[(6*a*a*d,(2,0,0)),(12*a*b*c,(1,1,1)),(2*b**3,(0,3,0))]
 qua=[(2*a**4,(4,0,0)),(24*a*b*b*c,(1,-2,1)),(48*a*b*c*d,(-1,1,1)),(8*a*c**3,(1,0,-3)),(24*b*c*c*d,(0,1,-2))]
 return [(coef/sr.I(6),n) for coef,n in cub]+[(coef/sr.I(24),n) for coef,n in qua]
TERMS=terms()
def Ibox(c,r):
 c=F(str(c));r=F(str(r));return sr.iv.mpf([sr.I(c-r),sr.I(c+r)])
def hess(phi):
 H=[[sr.I(0) for _ in range(3)] for __ in range(3)]
 for coef,n in TERMS:
  ang=sum((sr.I(n[k])*phi[k] for k in range(3)),sr.I(0));co=sr.iv.cos(ang)
  for i in range(3):
   for j in range(3):H[i][j]+=-coef*sr.I(n[i]*n[j])*co
 return H
def qbound(center,r):
 c=[sr.I(F(str(x))) for x in center]; H0=hess(c);Hm=np.array([[sr.mid(x) for x in row] for row in H0],float);A=np.linalg.inv(Hm)
 # rationalized float A; exact nonsingularity check
 Ar=[[F(repr(float(A[i,j]))) for j in range(3)] for i in range(3)]
 det=(Ar[0][0]*(Ar[1][1]*Ar[2][2]-Ar[1][2]*Ar[2][1])-Ar[0][1]*(Ar[1][0]*Ar[2][2]-Ar[1][2]*Ar[2][0])+Ar[0][2]*(Ar[1][0]*Ar[2][1]-Ar[1][1]*Ar[2][0]))
 assert det!=0
 X=[Ibox(x,r) for x in center];HB=hess(X); rows=[]
 for i in range(3):
  s=F(0)
  for j in range(3):
   e=sr.I(int(i==j))-sum((sr.I(Ar[i][k])*HB[k][j] for k in range(3)),sr.I(0));lo,hi=sr.bounds(e);s+=max(abs(lo),abs(hi))
  rows.append(s)
 return max(rows),Ar

def run():
 roots=json.load(open(CAT))['roots']; radii=['0.05','0.03','0.02','0.015','0.01','0.0075','0.005','0.003','0.002','0.001']
 rows=[];st=time.time()
 for rec in roots:
  best=None;tests=[];bestA=None
  for r in radii:
   q,A=qbound(rec['phase'],r);tests.append({'radius':r,'q_inf':str(q),'pass':q<1})
   if q<1 and best is None:best=r;bestA=A
  rows.append({'id':rec['id'],'center':[str(x) for x in rec['phase']],'certified_radius':best,'tests':tests,'preconditioner':[[str(x) for x in row] for row in bestA] if bestA else None})
 out={'task':'R7N-035','method':'uniform infinity-norm injectivity: sup ||I-A H(phi)||_inf < 1 on lifted convex phase box; existence inherited from certified 1e-7 root box','count':len(rows),'certified_count':sum(x['certified_radius'] is not None for x in rows),'radius_histogram':{},'elapsed_seconds':time.time()-st,'roots':rows}
 for x in rows:out['radius_histogram'][str(x['certified_radius'])]=out['radius_histogram'].get(str(x['certified_radius']),0)+1
 (ROOT/'results/R7N-035_quartic_uniqueness_collars.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='roots'},indent=2))
if __name__=='__main__':run()
