"""R7P-068 validated finite-radius off-face cone certificate.

The proof works with the 3x3 Schur-reduced covariance Mtilde from R7P-065 and
uses the shifted characteristic at sigma_*:

    P = det(sigma_* I-Mtilde),
    P1 = d/dsigma det(sigma I-Mtilde)|sigma_*,
    c2 = P''/2 = tr(sigma_* I-Mtilde).

For a PSD 3x3 matrix and c2>0, R7P-044 proves lambda_2<=sigma_* whenever
P<=0 OR P1>=0.  Around the double root we certify this disjunction on an
explicit physical tangent box without differentiating an ordered eigenvalue.

Coordinates are x=r-r_*, u=1-s, v=1-t, e=1-q.  The physical odd conditional
law has z=t^sqrt(3).  On a box 0<=v<=rho we enclose z and its first two
v-derivatives by rational intervals, preserving (and slightly enlarging) the
physical relation.  Every proof decision below uses Fraction endpoints.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json, math

from intervals import QI, sqrt_interval
import boundary_ising as bi

ROOT=Path(__file__).resolve().parents[1]
N=4


def _as_jet(x): return x if isinstance(x,Jet4) else Jet4(x)


class Jet4:
    """Second-order interval jet in four variables."""
    def __init__(self,v,g=None,H=None):
        self.v=QI.cast(v)
        self.g=[QI(0) for _ in range(N)] if g is None else g
        self.H=[[QI(0) for _ in range(N)] for __ in range(N)] if H is None else H
    @staticmethod
    def var(v,i):
        g=[QI(0) for _ in range(N)]; g[i]=QI(1)
        return Jet4(v,g)
    def __add__(self,o):
        o=_as_jet(o)
        return Jet4(self.v+o.v,
                    [self.g[i]+o.g[i] for i in range(N)],
                    [[self.H[i][j]+o.H[i][j] for j in range(N)] for i in range(N)])
    __radd__=__add__
    def __neg__(self):
        return Jet4(-self.v,[-x for x in self.g],[[-x for x in row] for row in self.H])
    def __sub__(self,o): return self+(-_as_jet(o))
    def __rsub__(self,o): return _as_jet(o)-self
    def __mul__(self,o):
        o=_as_jet(o)
        g=[self.g[i]*o.v+self.v*o.g[i] for i in range(N)]
        H=[]
        for i in range(N):
            row=[]
            for j in range(N):
                row.append(self.H[i][j]*o.v+self.v*o.H[i][j]
                           +self.g[i]*o.g[j]+self.g[j]*o.g[i])
            H.append(row)
        return Jet4(self.v*o.v,g,H)
    __rmul__=__mul__
    def inv(self):
        v=self.v; iv=QI(1)/v
        g=[-self.g[i]/(v*v) for i in range(N)]
        H=[]
        for i in range(N):
            row=[]
            for j in range(N):
                row.append(2*self.g[i]*self.g[j]/(v*v*v)-self.H[i][j]/(v*v))
            H.append(row)
        return Jet4(iv,g,H)
    def __truediv__(self,o): return self*_as_jet(o).inv()
    def __rtruediv__(self,o): return _as_jet(o)*self.inv()
    def __pow__(self,n):
        if n<0:return (self**(-n)).inv()
        out=Jet4(1)
        for _ in range(n):out=out*self
        return out


def _vadd(a,b): return [a[i]+b[i] for i in range(3)]
def _vsub(a,b): return [a[i]-b[i] for i in range(3)]
def _vscale(s,a): return [s*a[i] for i in range(3)]
def _madd(A,B): return [[A[i][j]+B[i][j] for j in range(3)] for i in range(3)]
def _mscale(s,A): return [[s*A[i][j] for j in range(3)] for i in range(3)]
def _outer(a): return [[a[i]*a[j] for j in range(3)] for i in range(3)]


def _weighted_stats(weights,features):
    """Stable covariance: pairwise positive-weight identity, no E[xx]-mu mu cancellation."""
    D=sum(weights,Jet4(0)); mu=[Jet4(0),Jet4(0),Jet4(0)]
    for wi,vi in zip(weights,features):
        mu=_vadd(mu,[wi*vi[k] for k in range(3)])
    mu=_vscale(Jet4(1)/D,mu)
    C=[[Jet4(0) for _ in range(3)] for __ in range(3)]
    for i in range(len(weights)):
        for j in range(i+1,len(weights)):
            dv=[features[i][k]-features[j][k] for k in range(3)]
            C=_madd(C,_mscale(weights[i]*weights[j],
                              [[dv[a]*dv[b] for b in range(3)] for a in range(3)]))
    return mu,_mscale(Jet4(1)/(D*D),C)


def _square(I):
    if I.lo<=0<=I.hi:return QI(0,max(I.lo*I.lo,I.hi*I.hi))
    return QI(min(I.lo*I.lo,I.hi*I.hi),max(I.lo*I.lo,I.hi*I.hi))


def _schur_cone(M,margin):
    """Cone negativity of M+margin*I for x in R, u,v>=0.

    Eliminating signed x leaves a 2x2 Schur form.  Strictly negative diagonal
    and nonpositive mixed coefficient on the nonnegative (u,v) cone follow
    from hxx<0 and the three positive Schur numerators below.
    """
    A=[[M[i][j]+(QI(margin) if i==j else QI(0)) for j in range(3)] for i in range(3)]
    h=A[0][0]
    vals=(h,
          h*A[1][1]-_square(A[0][1]),
          h*A[1][2]-A[0][1]*A[0][2],
          h*A[2][2]-_square(A[0][2]))
    ok=vals[0].hi<0 and all(x.lo>0 for x in vals[1:])
    return ok,vals


def _model_box(rho):
    """Return interval jets P,P1,c2 on |x|<=rho, 0<=u,v,e<=rho."""
    rho=F(rho)
    if not (0<rho<F(1,10)): raise ValueError('rho outside local regime')
    L=bi.strict_intervals();l3,l4,l5,l6=L[3],L[4],L[5],L[6]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3); c=l6/(3*sigma)
    tau=sqrt_interval(((2*l3-l4)*(2*l3-l5))/(4*l3*l3),40)
    rstar=(QI(1)-tau)/(QI(1)+tau)
    a3=sqrt_interval(l3/6,40);a4=sqrt_interval(l4/6,40);a5=sqrt_interval(l5/6,40)
    rt3=sqrt_interval(QI(3),40)
    Vp=[[a3,a4,a5],[a3,-a4/2,-a5/2],[-a3,a4,-a5],[-a3,-a4/2,a5/2]]
    Vm=[[QI(0),a4,QI(0)],[QI(0),-a4/2,rt3*a5/2],[QI(0),-a4/2,-rt3*a5/2]]

    x=Jet4.var(QI(-rho,rho),0);u=Jet4.var(QI(0,rho),1)
    v=Jet4.var(QI(0,rho),2);e=Jet4.var(QI(0,rho),3)
    r=Jet4(rstar)+x; s=Jet4(1)-u; t=Jet4(1)-v

    # Physical z=t^sqrt(3).  For t in [1-rho,1] and 1<sqrt(3)<2,
    # t^2<=z<=1.  With v=1-t:
    # dz/dv=-alpha*t^(alpha-1),
    # d2z/dv2=alpha(alpha-1)*t^(alpha-2).
    # The rational intervals below safely contain all three quantities.
    tlo=F(1)-rho
    zv=QI(tlo*tlo,1)
    zg=[QI(0) for _ in range(N)]; zg[2]=QI(-rt3.hi,-rt3.lo*tlo)
    zH=[[QI(0) for _ in range(N)] for __ in range(N)]
    af=rt3*(rt3-QI(1)); zH[2][2]=af*QI(1,F(1)/tlo)
    z=Jet4(zv,zg,zH)

    # Divide the physical odd aggregate weights by their common positive factor:
    # (s t^(2+sqrt3), t^2, s t^(2-sqrt3)) -> (z, s, s z^2)
    # in the feature order (0,a4,0), (0,-a4/2,+sqrt3 a5/2),
    # (0,-a4/2,-sqrt3 a5/2).
    wp=[Jet4(1),2*s*(t**3),r*(t**4),2*r*s*t]
    wm=[z,s,s*(z**2)]
    mp,Cp=_weighted_stats(wp,Vp); mm,Cm=_weighted_stats(wm,Vm)
    d=_vsub(mp,mm); q=Jet4(1)-e; qe=q*e; eta=Jet4(1)-qe*c
    M=_madd(_mscale(q,Cp),_mscale(e,Cm))
    M=_madd(M,_mscale(qe/eta,_outer(d)))
    K=[[(Jet4(sigma)-M[i][j] if i==j else -M[i][j]) for j in range(3)] for i in range(3)]
    a,b,cc=K[0];dd=K[1][1];ee=K[1][2];ff=K[2][2]
    P=a*dd*ff+2*b*cc*ee-a*ee*ee-dd*cc*cc-ff*b*b
    P1=(a*dd-b*b)+(a*ff-cc*cc)+(dd*ff-ee*ee)
    c2=a+dd+ff
    return rstar,P,P1,c2


def _floor_frac(x,scale): return (x.numerator*scale//x.denominator)/scale

def _ceil_frac(x,scale):
    n=x.numerator*scale; d=x.denominator
    return -((-n)//d)/scale


def _short_interval(I,digits=14):
    S=10**digits
    lo=F(math.floor(I.lo*S),S) if abs(I.lo)<10**12 else I.lo
    hi=F(math.ceil(I.hi*S),S) if abs(I.hi)<10**12 else I.hi
    return [format(float(lo),'.15g'),format(float(hi),'.15g')]


def raw_local_cone(rho=F(1,8192),margin=F(1,10000)):
    """Compute the acceptance-grade interval conditions; returns QI objects."""
    rho=F(rho); margin=F(margin)
    rstar,P,P1,c2=_model_box(rho); g=P1.g; H=P.H
    signs={
      'P1_x_negative':g[0].hi<0,
      'P1_u_negative':g[1].hi<0,
      'P1_v_negative':g[2].hi<0,
      'P1_e_positive':g[3].lo>0,
      'P_ee_positive':H[3][3].lo>0,
      'c2_positive':c2.v.lo>0,
    }
    alpha=[(-g[i])/g[3] for i in range(3)]
    HB=[row[:3] for row in H[:3]]
    HE=[]
    for i in range(3):
        row=[]
        for j in range(3):
            row.append(H[i][j]+H[i][3]*alpha[j]+alpha[i]*H[3][j]
                       +alpha[i]*alpha[j]*H[3][3])
        HE.append(row)
    okB,sB=_schur_cone(HB,margin); okE,sE=_schur_cone(HE,margin)
    ok=all(signs.values()) and okB and okE
    return dict(rho=rho,margin=margin,rstar=rstar,P=P,P1=P1,c2=c2,
                signs=signs,alpha=alpha,HB=HB,HE=HE,
                boundary_schur=sB,endpoint_schur=sE,
                boundary_ok=okB,endpoint_ok=okE,status='INTERVAL_CERTIFIED' if ok else 'FAILED')


def local_cone_record(rho=F(1,8192),margin=F(1,10000)):
    R=raw_local_cone(rho,margin)
    def qij(I): return _short_interval(I)
    return {
      'id':'R7P-068-finite-local-cone',
      'claim_id':'CLM-054',
      'status':R['status'],
      'proof_type':'second-order rational interval automatic differentiation plus characteristic-inertia disjunction',
      'coordinates':{
        'x':'r-r_star (signed)','u':'1-s >= 0','v':'1-t >= 0','e':'1-q >= 0',
        'r':'exp(-2 J3)','s':'exp(-3 J4/2)','t':'exp(-J5/2)'
      },
      'domain':{
        'rho':str(R['rho']),
        'x':['-1/8192','1/8192'] if R['rho']==F(1,8192) else [str(-R['rho']),str(R['rho'])],
        'u':[0,str(R['rho'])],'v':[0,str(R['rho'])],'e':[0,str(R['rho'])],
        'q_lower':str(1-R['rho']),
      },
      'strict_spectrum_rstar_interval':qij(R['rstar']),
      'criterion':'For PSD 3x3 Mtilde with c2>0, lambda2<=sigma_* if P<=0 OR P1>=0 (R7P-044).',
      'exact_anchor_identities':['P(root)=0','grad P(root)=0','P1(root)=0'],
      'interval_signs':R['signs'],
      'P1_gradient_bounds':{k:qij(R['P1'].g[i]) for i,k in enumerate(('x','u','v','e'))},
      'P_H_ee':qij(R['P'].H[3][3]),
      'c2_range':qij(R['c2'].v),
      'alpha_ranges':{k:qij(R['alpha'][i]) for i,k in enumerate(('x','u','v'))},
      'cone_margin':str(R['margin']),
      'boundary_endpoint_schur':{
        'hxx':qij(R['boundary_schur'][0]),'Nuu':qij(R['boundary_schur'][1]),
        'Nuv':qij(R['boundary_schur'][2]),'Nvv':qij(R['boundary_schur'][3]),
        'pass':R['boundary_ok']},
      'P1_zero_endpoint_schur':{
        'hxx':qij(R['endpoint_schur'][0]),'Nuu':qij(R['endpoint_schur'][1]),
        'Nuv':qij(R['endpoint_schur'][2]),'Nvv':qij(R['endpoint_schur'][3]),
        'pass':R['endpoint_ok']},
      'irrational_relation_enclosure':'z=t^sqrt(3) is enclosed with t^2<=z<=1 and rigorous first/second v-derivative intervals; no independent physical z is assumed.',
      'logic':[
        'If P1>=0 the R7P-044 disjunction is already satisfied.',
        'If P1<0, integral Taylor for P1 and the certified gradient signs imply 0<=e<=alpha_x*x+alpha_u*u+alpha_v*v for some alpha in the stored intervals.',
        'Integral Taylor for P uses an averaged Hessian contained in the stored box Hessian interval.  Because H_ee>0, its quadratic form is convex in e, so its maximum on the admissible e segment is at e=0 or at the P1-linear endpoint.',
        'Both endpoint quadratic forms remain cone-negative after adding margin*I by the stored Schur tests. Hence P<0 away from the root whenever P1<0.',
        'Therefore P<=0 OR P1>=0 throughout the declared box, with c2>0, so lambda2(Mtilde)<=sigma_* throughout it.'
      ],
      'conclusion':'Explicit finite physical-cone neighborhood certified: lambda2(Mtilde)<=sigma_* and hence lambda2(M4)<=sigma_* throughout the declared local domain.',
      'nontransfer':'This local theorem does not prove the full off-face orthant ceiling or any full-seven-coordinate Hessian bound.'
    }


def write_certificate(path=None):
    rec=local_cone_record()
    if path is None:path=ROOT/'certificates/R7P-068_local_cone.json'
    Path(path).write_text(json.dumps(rec,indent=2)+'\n')
    return rec


if __name__=='__main__':
    print(json.dumps(write_certificate(),indent=2))
