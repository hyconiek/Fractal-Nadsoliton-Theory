from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import json, math, hashlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
import sys
sys.path.insert(0, str(HERE))
import intervals
from intervals import QI, sqrt_interval

N=4
TAU=F(67,250)
AL=(19,11); AH=(26,15)
_ROUNDING=None
_ORIG_QI_INIT=intervals.QI.__init__

def enable_outward_rounding(digits:int=9):
    global _ROUNDING
    if _ROUNDING is not None:
        if _ROUNDING != digits:
            raise RuntimeError(f'rounding already initialized at {_ROUNDING}')
        return
    q=10**digits
    original=_ORIG_QI_INIT
    def initialize(self,lo,hi=None):
        original(self,lo,hi)
        self.lo=F((self.lo*q).__floor__(),q)
        self.hi=F((self.hi*q).__ceil__(),q)
    intervals.QI.__init__=initialize
    _ROUNDING=digits
    for a,b in [(F(1,3),F(2,3)),(F(-2,7),F(-1,7)),(F(0),F(0))]:
        v=QI(a,b); assert v.lo<=a<=b<=v.hi

def _as_jet(x): return x if isinstance(x,Jet4) else Jet4(x)
class Jet4:
    def __init__(self,v,g=None,H=None):
        self.v=QI.cast(v)
        self.g=[QI(0) for _ in range(N)] if g is None else g
        self.H=[[QI(0) for _ in range(N)] for __ in range(N)] if H is None else H
    @staticmethod
    def var(v,i):
        g=[QI(0) for _ in range(N)]; g[i]=QI(1); return Jet4(v,g)
    def __add__(self,o):
        o=_as_jet(o); return Jet4(self.v+o.v,[self.g[i]+o.g[i] for i in range(N)],[[self.H[i][j]+o.H[i][j] for j in range(N)] for i in range(N)])
    __radd__=__add__
    def __neg__(self): return Jet4(-self.v,[-x for x in self.g],[[-x for x in row] for row in self.H])
    def __sub__(self,o): return self+(-_as_jet(o))
    def __rsub__(self,o): return _as_jet(o)-self
    def __mul__(self,o):
        o=_as_jet(o)
        g=[self.g[i]*o.v+self.v*o.g[i] for i in range(N)]
        H=[[self.H[i][j]*o.v+self.v*o.H[i][j]+self.g[i]*o.g[j]+self.g[j]*o.g[i] for j in range(N)] for i in range(N)]
        return Jet4(self.v*o.v,g,H)
    __rmul__=__mul__
    def inv(self):
        v=self.v; iv=QI(1)/v
        g=[-self.g[i]/(v*v) for i in range(N)]
        H=[[2*self.g[i]*self.g[j]/(v*v*v)-self.H[i][j]/(v*v) for j in range(N)] for i in range(N)]
        return Jet4(iv,g,H)
    def __truediv__(self,o): return self*_as_jet(o).inv()
    def __rtruediv__(self,o): return _as_jet(o)*self.inv()
    def __pow__(self,n):
        if n<0:return (self**(-n)).inv()
        out=Jet4(1)
        for _ in range(n): out=out*self
        return out

@lru_cache(maxsize=200000)
def pow_frac_point(x:F,p:int,q:int,digits:int=18):
    x=F(x)
    if x<0: raise ValueError('negative base')
    if x==0:
        if p>0:return QI(0)
        raise ValueError('zero negative power')
    if p<0:
        I=pow_frac_point(x,-p,q,digits); return QI(1/I.hi,1/I.lo)
    Q=10**digits; Nn,D=x.numerator,x.denominator
    targetN=pow(Nn,p)*pow(Q,q); targetD=pow(D,p)
    approx=float(x)**(p/q); n=max(0,int(math.floor(approx*Q)))
    def le(k): return pow(k,q)*targetD <= targetN
    while n>0 and not le(n): n-=1
    while le(n+1): n+=1
    if pow(n,q)*targetD==targetN:return QI(F(n,Q))
    return QI(F(n,Q),F(n+1,Q))

def _z_of_v(v:Jet4,vlo:F,vhi:F,rt3:QI):
    tlo=F(1)-F(vhi); thi=F(1)-F(vlo)
    if not (0<tlo<=thi<=1): raise ValueError('v box outside 0<=v<1')
    zv=QI(pow_frac_point(tlo,*AH).lo, pow_frac_point(thi,*AL).hi)
    d1=QI(rt3.lo*tlo,rt3.hi)
    g=[QI(0) for _ in range(N)]; g[2]=-d1
    af=rt3*(rt3-QI(1)); H=[[QI(0) for _ in range(N)] for __ in range(N)]; H[2][2]=af*QI(1,F(1)/tlo)
    return Jet4(zv,g,H)

def _absmax(I): return max(abs(I.lo),abs(I.hi))
def _range_from_jets(full,cen,rads):
    lin=cen.v
    for i,r in enumerate(rads): lin=lin+cen.g[i]*QI(-r,r)
    rem=F(0)
    for i in range(N):
        for j in range(N): rem += F(1,2)*_absmax(full.H[i][j])*rads[i]*rads[j]
    return QI(lin.lo-rem,lin.hi+rem)

def _internal_bounds(cell):
    cell=tuple((F(lo),F(hi)) for lo,hi in cell)
    (rlo,rhi),(slo,shi),(tlo,thi),(ylo,yhi)=cell
    A=sqrt_interval(QI(rlo,rhi),60)
    return ((A.lo,A.hi),(F(1)-shi,F(1)-slo),(F(1)-thi,F(1)-tlo),(ylo,yhi))
def _midpoint_bounds(bounds):
    mids=[]; rads=[]
    for lo,hi in bounds:
        m=(lo+hi)/2; mids.append((m,m)); rads.append((hi-lo)/2)
    return tuple(mids),rads

@lru_cache(maxsize=1)
def _load_obs():
    d=json.loads((ROOT/'inherited/spectral_obs_rounded9.json').read_text())
    return [[QI(F(lo),F(hi)) for lo,hi in row] for row in d['obs']]

def _weights_jet(bounds):
    (Alo,Ahi),(ulo,uhi),(vlo,vhi),(ylo,yhi)=bounds
    A=Jet4.var(QI(Alo,Ahi),0);u=Jet4.var(QI(ulo,uhi),1);v=Jet4.var(QI(vlo,vhi),2);y=Jet4.var(QI(ylo,yhi),3)
    s=Jet4(1)-u; t=Jet4(1)-v
    rt3=sqrt_interval(QI(3),40); z=_z_of_v(v,vlo,vhi,rt3); r=A*A
    return [Jet4(1),2*s*(t**3),r*(t**4),2*r*s*t,2*A*(t**2)*y,2*A*s*(t**2)*y/z,2*A*s*(t**2)*z*y]

def _moment_jets(bounds,B,c,OBS):
    w=_weights_jet(bounds); D=sum(w,Jet4(0)); invD=D.inv(); Z=[]
    for i in range(7):
        row=[]
        for k in range(3):
            zz=QI(0)
            for r in range(4): zz += B[r][k]*OBS[i][r]
            row.append(zz-QI(c[k]))
        Z.append(row)
    E=[[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(a,3):
            Nn=Jet4(0)
            for i in range(7): Nn += w[i]*(Z[i][a]*Z[i][b])
            E[a][b]=Nn*invD;E[b][a]=E[a][b]
    return E

def _pd_result(K):
    d1=K[0][0]
    d2=K[0][0]*K[1][1]-K[0][1]*K[0][1]
    d3=(K[0][0]*K[1][1]*K[2][2]+2*K[0][1]*K[0][2]*K[1][2]-K[0][0]*K[1][2]*K[1][2]-K[1][1]*K[0][2]*K[0][2]-K[2][2]*K[0][1]*K[0][1])
    sylv=d1.lo>0 and d2.lo>0 and d3.lo>0
    g=[]
    for a in range(3):
        offsum=sum(max(abs(K[a][b].lo),abs(K[a][b].hi)) for b in range(3) if b!=a)
        g.append(K[a][a].lo-offsum)
    gersh=min(g)>0
    return sylv,gersh,d1,d2,d3,g

def certify_fixed(cert, digits=9):
    enable_outward_rounding(digits)
    if cert.get('threshold')!='67/250': return {'ok':False,'reason':'WRONG_THRESHOLD'}
    cell=tuple((F(a),F(b)) for a,b in cert['cell'])
    den=int(cert['basis_den']); Bn=cert['basis_num']; B=[[F(int(Bn[i][j]),den) for j in range(3)] for i in range(4)]; c=[F(x) for x in cert['center_c']]
    # exact rank evidence and chosen 3x3 minor
    rows=tuple(cert['rank_minor_rows']); A=[[Bn[r][j] for j in range(3)] for r in rows]
    det=(A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))
    detq=F(det,den**3)
    if detq==0:return {'ok':False,'reason':'RANK'}
    if str(detq)!=cert['rank_minor_exact']:return {'ok':False,'reason':'RANK_EVIDENCE_MISMATCH','computed':str(detq)}
    OBS=_load_obs(); bounds=_internal_bounds(cell);cb,rads=_midpoint_bounds(bounds)
    full=_moment_jets(bounds,B,c,OBS);cen=_moment_jets(cb,B,c,OBS);Er=[[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(a,3):Er[a][b]=_range_from_jets(full[a][b],cen[a][b],rads);Er[b][a]=Er[a][b]
    Gram=[[sum(B[r][a]*B[r][b] for r in range(4)) for b in range(3)] for a in range(3)]
    K=[[QI(TAU*Gram[a][b])-Er[a][b] for b in range(3)] for a in range(3)]
    sylv,gersh,d1,d2,d3,g=_pd_result(K);ok=sylv or gersh
    return {'ok':bool(ok),'reason':'PHYSICAL_CENTERED_SYLVESTER_PD' if sylv else ('PHYSICAL_CENTERED_GERSHGORIN_PD' if gersh else 'FAILED'),
            'rank_minor_exact':str(detq),'gram_matrix_exact':[[str(x) for x in row] for row in Gram],
            'internal_chart_bounds':[[str(a),str(b)] for a,b in bounds],
            'moment_entry_enclosures':[[[str(Er[a][b].lo),str(Er[a][b].hi)] for b in range(3)] for a in range(3)],
            'pd_bounds_exact':{'d1':[str(d1.lo),str(d1.hi)],'d2':[str(d2.lo),str(d2.hi)],'d3':[str(d3.lo),str(d3.hi)],'gersh_lower_bounds':[str(x) for x in g]}}

def exact_match(cert,out):
    keys=['gram_matrix_exact','internal_chart_bounds','moment_entry_enclosures','pd_bounds_exact']
    return {k:(cert.get(k)==out.get(k)) for k in keys}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('certificate_json');ap.add_argument('--digits',type=int,default=9);args=ap.parse_args()
    cert=json.loads(Path(args.certificate_json).read_text());o=certify_fixed(cert,args.digits);o['exact_match']=exact_match(cert,o);print(json.dumps(o,indent=2))
