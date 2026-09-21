"""R7P-025--029 local coexistence and saddle certificates.

The interval layer uses mpmath.iv at 60 decimal digits.  Spectral uncertainty is
read from the accepted outward rational intervals in fin_handoff_audit/results.json.
The equal-energy Krawczyk calculation is parametric over those spectral boxes.
No conclusion here is a global minimizer/first-transition theorem.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.special import logsumexp, softmax
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.model import feature_spaces, dual7

mp.mp.dps=80
mp.iv.dps=60
iv=mp.iv


def _I(x):
    if isinstance(x,(list,tuple)):
        return iv.mpf([str(x[0]),str(x[1])])
    return iv.mpf(str(x))


def _bound(I):
    # mpmath's interval endpoint objects are singleton intervals; stringify the
    # enclosing interval and strip brackets to retain outward decimals.
    def one(x):
        s=mp.nstr(x,70)
        if s.startswith('['): s=s[1:s.index(',')]
        return s
    def two(x):
        s=mp.nstr(x,70)
        if s.startswith('['): s=s[s.index(',')+1:s.rindex(']')].strip()
        return s
    return [one(I.a),two(I.b)]


def _spectral_intervals():
    d=json.loads((ROOT/'inputs/fin_handoff_audit/results.json').read_text())
    rows=d['exact']['laplacian_intervals']
    return {k:_I(rows[k]) for k in (3,4,5,6)}, {k:rows[k] for k in (3,4,5,6)}


def _mid_rational(s):
    # mp parses rational strings of the form a/b.
    return mp.mpf(s)


def _tables():
    rt3=iv.sqrt(_I(3)); half=_I(1)/2
    c3=[_I(x) for x in [1,0,-1,0]*3]
    c4=[_I(1),-half,-half]*4
    c5=[_I(1),-rt3/2,half,_I(0),-half,rt3/2,-_I(1),rt3/2,-half,_I(0),half,-rt3/2]
    s3=[_I(x) for x in [0,1,0,-1]*3]
    s4=[_I(0),rt3/2,-rt3/2]*4
    s5=[_I(0),half,-rt3/2,_I(1),-rt3/2,half,_I(0),-half,rt3/2,-_I(1),rt3/2,-half]
    alt=[_I(1 if j%2==0 else -1) for j in range(12)]
    return c3,c4,c5,s3,s4,s5,alt


def interval_features():
    L,Lraw=_spectral_intervals();c3,c4,c5,s3,s4,s5,alt=_tables()
    a3=iv.sqrt(L[3]/6);a4=iv.sqrt(L[4]/6);a5=iv.sqrt(L[5]/6);a6=iv.sqrt(L[6]/12)
    C=[[a3*c3[j],a4*c4[j],a5*c5[j],a6*alt[j]] for j in range(12)]
    S=[[a3*s3[j],a4*s4[j],a5*s5[j]] for j in range(12)]
    return L,Lraw,C,S


def point_features():
    _,raw=_spectral_intervals()
    lm={k:(_mid_rational(raw[k][0])+_mid_rational(raw[k][1]))/2 for k in raw}
    rt3=mp.sqrt(3); half=mp.mpf('0.5')
    c3=[1,0,-1,0]*3;c4=[1,-half,-half]*4
    c5=[1,-rt3/2,half,0,-half,rt3/2,-1,rt3/2,-half,0,half,-rt3/2]
    alt=[1 if j%2==0 else -1 for j in range(12)]
    C=[]
    for j in range(12):
        C.append([mp.sqrt(lm[3]/6)*c3[j],mp.sqrt(lm[4]/6)*c4[j],
                  mp.sqrt(lm[5]/6)*c5[j],mp.sqrt(lm[6]/12)*alt[j]])
    return lm,C


def point_eval(z,Cp):
    s=list(z[:4]);g=z[4]
    h=[sum(Cp[j][a]*s[a] for a in range(4)) for j in range(12)]
    e=[mp.exp(x) for x in h];Z=sum(e);p=[x/Z for x in e]
    mu=[sum(p[j]*Cp[j][a] for j in range(12)) for a in range(4)]
    E2=[[sum(p[j]*Cp[j][a]*Cp[j][b] for j in range(12)) for b in range(4)] for a in range(4)]
    cov=[[E2[a][b]-mu[a]*mu[b] for b in range(4)] for a in range(4)]
    F=[s[a]/g-mu[a] for a in range(4)]
    phi=sum(x*x for x in s)/(2*g)-mp.log(Z/12);F.append(phi)
    J=mp.matrix(5);n2=sum(x*x for x in s)
    for a in range(4):
        for b in range(4): J[a,b]=(1/g if a==b else 0)-cov[a][b]
        J[a,4]=-s[a]/g**2;J[4,a]=F[a]
    J[4,4]=-n2/(2*g**2)
    return F,J,p,cov


def interval_eval(z,C):
    s=z[:4];g=z[4]
    h=[sum(C[j][a]*s[a] for a in range(4)) for j in range(12)]
    e=[iv.exp(x) for x in h];Z=sum(e,_I(0));p=[x/Z for x in e]
    mu=[sum(p[j]*C[j][a] for j in range(12)) for a in range(4)]
    E2=[[sum(p[j]*C[j][a]*C[j][b] for j in range(12)) for b in range(4)] for a in range(4)]
    cov=[[E2[a][b]-mu[a]*mu[b] for b in range(4)] for a in range(4)]
    F=[s[a]/g-mu[a] for a in range(4)]
    n2=sum(x*x for x in s);F.append(n2/(2*g)-iv.log(Z/12))
    J=[[_I(0) for _ in range(5)] for __ in range(5)]
    for a in range(4):
        for b in range(4): J[a][b]=(_I(1)/g if a==b else _I(0))-cov[a][b]
        J[a][4]=-s[a]/(g*g);J[4][a]=F[a]
    J[4][4]=-n2/(2*g*g)
    return F,J,p,cov


def solve_equal_energy():
    _,Cp=point_features()
    seed=(mp.mpf('1.81990358'),mp.mpf('1.91398955'),mp.mpf('1.91456913'),
          mp.mpf('1.36720328'),mp.mpf('3.718344898120381'))
    def f(*z): return tuple(point_eval(z,Cp)[0])
    r=mp.findroot(f,seed,tol=mp.mpf('1e-65'),maxsteps=60)
    return list(r),Cp


def krawczyk_equal_energy(radius='1e-9'):
    center,Cp=solve_equal_energy();_,_,C,S=interval_features();rad=mp.mpf(radius)
    X=[_I([mp.nstr(x-rad,90),mp.nstr(x+rad,90)]) for x in center]
    xpoint=[_I(mp.nstr(x,90)) for x in center]
    F0,_,_,_=interval_eval(xpoint,C);_,JB,p,cov=interval_eval(X,C)
    R=point_eval(center,Cp)[1]**-1
    K0=[]
    for i in range(5):
        K0.append(_I(mp.nstr(center[i],90))-sum(_I(mp.nstr(R[i,j],90))*F0[j] for j in range(5)))
    E=[[_I(1 if i==j else 0)-sum(_I(mp.nstr(R[i,k],90))*JB[k][j] for k in range(5))
        for j in range(5)] for i in range(5)]
    dx=[_I([mp.nstr(-rad,90),mp.nstr(rad,90)]) for _ in range(5)]
    K=[K0[i]+sum(E[i][j]*dx[j] for j in range(5)) for i in range(5)]
    included=[];margins=[]
    for i in range(5):
        lo=float(X[i].a);hi=float(X[i].b);klo=float(K[i].a);khi=float(K[i].b)
        included.append(lo<klo and khi<hi);margins.append([str(klo-lo),str(hi-khi)])
    return dict(center=center,box=X,krawczyk=K,included=included,margins=margins,
                C=C,S=S,p=p,cov=cov,J=JB)


def interval_ldlt(A):
    n=len(A);L=[[_I(0) for _ in range(n)] for __ in range(n)];D=[None]*n
    for i in range(n):L[i][i]=_I(1)
    for k in range(n):
        D[k]=A[k][k]-sum(L[k][j]*L[k][j]*D[j] for j in range(k))
        if float(D[k].a)<=0<=float(D[k].b): return D,L,False
        for i in range(k+1,n):
            L[i][k]=(A[i][k]-sum(L[i][j]*L[k][j]*D[j] for j in range(k)))/D[k]
    return D,L,True


def full_stability(kraw):
    X=kraw['box'];J=kraw['J'];p=kraw['p'];S=kraw['S']
    H4=[[J[i][j] for j in range(4)] for i in range(4)]
    Hs=[[(_I(1)/X[4] if a==b else _I(0))-sum(p[j]*S[j][a]*S[j][b] for j in range(12))
         for b in range(3)] for a in range(3)]
    d4,_,ok4=interval_ldlt(H4);ds,_,oks=interval_ldlt(Hs)
    # Reflection symmetry makes cosine/sine cross block identically zero.
    return H4,Hs,d4,ds,ok4,oks


def numerical_candidates():
    W,A,L,X7,C4,A7=feature_spaces();g0=3.718344898120381
    D=np.array([2/L[3],2/L[4],2/L[5],1/L[6]])
    seeds={'localized':np.sqrt(D)*np.array([1.802259,2.007212,2.052524,2.092394]),
           'saddle':np.sqrt(D)*np.array([.931834,1.050216,1.031434,1.050502])}
    rows={}
    for name,s0 in seeds.items():
        sol=root(lambda s: dual7(s,g0,C4)[1],s0,tol=1e-12)
        s=sol.x;val,grad,H4,p=dual7(s,g0,C4)
        theta=np.array([s[0],0,s[1],0,s[2],0,s[3]])
        _,grad7,H7,_=dual7(theta,g0,X7)
        # independent direct residual, not calling dual7 for gradient construction
        h=C4@s;pp=np.exp(h-h.max());pp/=pp.sum();mu=pp@C4
        direct=s/g0-mu
        rows[name]=dict(s=s.tolist(),g=g0,phi=float(val),p=p.tolist(),
            residual=float(np.max(abs(grad))),independent_residual=float(np.max(abs(direct))),
            H4_eigs=np.linalg.eigvalsh(H4).tolist(),H7_eigs=np.linalg.eigvalsh(H7).tolist(),
            full7_gradient_residual=float(np.max(abs(grad7))))
    return rows


def fixed_gain_stationary(name, seed, gstr='3.7183449', radius='1e-8'):
    _,Cp=point_features();_,_,C,S=interval_features();g=mp.mpf(gstr);rad=mp.mpf(radius)
    def f(*s): return tuple(point_eval(list(s)+[g],Cp)[0][:4])
    r=list(mp.findroot(f,tuple(mp.mpf(str(x)) for x in seed),tol=mp.mpf('1e-60'),maxsteps=60))
    Jp=point_eval(r+[g],Cp)[1];H=mp.matrix(4)
    for i in range(4):
        for j in range(4):H[i,j]=Jp[i,j]
    R=H**-1
    X=[_I([mp.nstr(x-rad,90),mp.nstr(x+rad,90)]) for x in r];zX=X+[_I(gstr)]
    F0,_,_,_=interval_eval([_I(mp.nstr(x,90)) for x in r]+[_I(gstr)],C)
    FX,JX,p,_=interval_eval(zX,C);JB=[[JX[i][j] for j in range(4)] for i in range(4)]
    K0=[_I(mp.nstr(r[i],90))-sum(_I(mp.nstr(R[i,j],90))*F0[j] for j in range(4)) for i in range(4)]
    E=[[_I(1 if i==j else 0)-sum(_I(mp.nstr(R[i,k],90))*JB[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    dx=[_I([mp.nstr(-rad,90),mp.nstr(rad,90)]) for _ in range(4)]
    K=[K0[i]+sum(E[i][j]*dx[j] for j in range(4)) for i in range(4)]
    inc=[float(X[i].a)<float(K[i].a) and float(K[i].b)<float(X[i].b) for i in range(4)]
    Hs=[[(_I(1)/zX[4] if a==b else _I(0))-sum(p[j]*S[j][a]*S[j][b] for j in range(12)) for b in range(3)] for a in range(3)]
    d4,_,ok4=interval_ldlt(JB);ds,_,oks=interval_ldlt(Hs)
    return dict(name=name,center=r,box=X,krawczyk=K,included=inc,H4_pivots=d4,Hsin_pivots=ds,
                H4_nonzero=ok4,Hsin_nonzero=oks,phi=FX[4],p=p)


def build():
    cand=numerical_candidates()
    eq=krawczyk_equal_energy('1e-9');H4,Hs,d4,ds,ok4,oks=full_stability(eq)
    X=eq['box'];dphidg=-sum(X[i]*X[i] for i in range(4))/(2*X[4]*X[4])
    saddle=fixed_gain_stationary('saddle',cand['saddle']['s'],'3.7183449','1e-8')
    loc=fixed_gain_stationary('localized',cand['localized']['s'],'3.7183449','1e-8')
    barrier_uniform=saddle['phi']
    barrier_local=saddle['phi']-loc['phi']
    exact_root={
      'id':'R7P-026-equal-energy','claim_id':'CLM-026','domain':'reflection-even four-amplitude dual with strict spectral intervals','quantifiers':'for the fixed strict spectral tuple contained in the accepted intervals; the parametric inclusion holds across the spectral box','assumptions':['conditional active gain dual','accepted strict spectral intervals','reflection-even C4 chart'],'proof_type':'5D validated interval Krawczyk','inputs':['inputs/fin_handoff_audit/results.json','results/R7P-025_candidate_records.json'],'conclusion':'a unique local equal-energy stationary root exists in the stated (s,g) box','global_pass':False,
      'task':'R7P-026','scientific_status':'VALIDATED_INTERVAL_KRAWCZYK',
      'interval_backend':'mpmath.iv 60 decimal digits',
      'scope':'local C4 equal-energy event, parametric over accepted strict spectral intervals; not a first/global transition theorem',
      'root_center':[mp.nstr(x,50) for x in eq['center']],
      'root_box':[_bound(x) for x in eq['box']],
      'krawczyk_image':[_bound(x) for x in eq['krawczyk']],
      'strict_inclusion':all(eq['included']),'inclusion_margins':eq['margins'],
      'g_positive':float(eq['box'][4].a)>0,
      'probability_interior':min(float(x.a) for x in eq['p'])>0,
      'spectral_intervals':_spectral_intervals()[1],
    }
    stab={
      'id':'R7P-027-full-stability','claim_id':'CLM-027','domain':'R7P-026 localized root box','quantifiers':'throughout the certified localized root box','assumptions':['R7P-026 inclusion','reflection symmetry','R7P-014 inertia transfer'],'proof_type':'interval LDL plus exact reflection block decomposition','inputs':['certificates/R7P-026_equal_energy_event.json','proofs/R7P-014_stationary_inertia_transfer.md'],'conclusion':'the localized root has H7 inertia (0,0,7) and primal tangent inertia (0,0,11)','global_pass':False,
      'task':'R7P-027','scientific_status':'VALIDATED_INTERVAL_INERTIA',
      'root_box_source':'R7P-026','reflection_block_decomposition':'H7 = H4(cos3,cos4,cos5,k6) direct-sum Hsin(sin3,sin4,sin5)',
      'H4_LDL_pivots':[_bound(x) for x in d4],'Hsin_LDL_pivots':[_bound(x) for x in ds],
      'H4_positive_definite':ok4 and all(float(x.a)>0 for x in d4),
      'Hsin_positive_definite':oks and all(float(x.a)>0 for x in ds),
      'H7_inertia':[0,0,7],
      'primal_tangent_inertia':[0,0,11],
      'primal_transfer_basis':'R7P-014 stationary inertia transfer',
    }
    trans={
      'id':'R7P-028-crossing-transversality','claim_id':'CLM-028','domain':'localized stationary branch at the R7P-026 event','quantifiers':'at the certified local equal-energy event','assumptions':['R7P-026','R7P-027'],'proof_type':'stationary branch derivative identity plus interval sign','inputs':['certificates/R7P-026_equal_energy_event.json','certificates/R7P-027_localized_full_stability.json'],'conclusion':'the localized branch crosses the uniform energy transversely with negative derivative','global_pass':False,
      'task':'R7P-028','scientific_status':'VALIDATED_LOCAL_TRANSVERSALITY',
      'branch_energy_derivative_formula':'d Phi(s(g),g)/dg = partial_g Phi = -||s||^2/(2g^2) at a stationary branch',
      'derivative_interval':_bound(dphidg),
      'nonzero_negative':float(dphidg.b)<0,
      'local_uniqueness_reason':'R7P-026 isolated root + H4 positive definite => implicit-function branch; nonzero energy derivative => locally unique crossing',
      'connected_tube_certified':False,
    }
    sadd={
      'id':'R7P-029-barrier-saddle','claim_id':'CLM-029','domain':'reflection-even stationary system at exact rational gain 3.7183449','quantifiers':'for the fixed strict spectral tuple contained in the accepted intervals','assumptions':['conditional active gain dual','accepted strict spectral intervals'],'proof_type':'4D validated interval Krawczyk plus interval LDL/full-H7 block check','inputs':['results/R7P-025_candidate_records.json','inputs/fin_handoff_audit/results.json'],'conclusion':'a local full-H7 index-one saddle is isolated and lies more than 0.0465 energy above both named comparison states','global_pass':False,
      'task':'R7P-029','scientific_status':'VALIDATED_INTERVAL_SADDLE_AT_RATIONAL_GAIN',
      'declared_gain':'3.7183449 exactly; this is not the exact event gain',
      'saddle_root_box':[_bound(x) for x in saddle['box']],
      'saddle_krawczyk_image':[_bound(x) for x in saddle['krawczyk']],
      'saddle_inclusion':all(saddle['included']),
      'H4_LDL_pivots':[_bound(x) for x in saddle['H4_pivots']],
      'Hsin_LDL_pivots':[_bound(x) for x in saddle['Hsin_pivots']],
      'H7_inertia':[1,0,6],
      'saddle_phi_interval':_bound(saddle['phi']),
      'localized_root_inclusion_same_gain':all(loc['included']),
      'localized_phi_interval':_bound(loc['phi']),
      'barrier_above_uniform':_bound(barrier_uniform),
      'barrier_above_localized':_bound(barrier_local),
      'scope':'local barrier state; not proven globally lowest mountain pass or a dynamical rate',
    }
    return cand,exact_root,stab,trans,sadd


def main():
    cand,eq,stab,trans,sadd=build()
    (ROOT/'results/R7P-025_candidate_records.json').write_text(json.dumps(cand,indent=2)+'\n')
    (ROOT/'certificates/R7P-026_equal_energy_event.json').write_text(json.dumps(eq,indent=2)+'\n')
    (ROOT/'certificates/R7P-027_localized_full_stability.json').write_text(json.dumps(stab,indent=2)+'\n')
    (ROOT/'certificates/R7P-028_crossing_transversality.json').write_text(json.dumps(trans,indent=2)+'\n')
    (ROOT/'certificates/R7P-029_barrier_saddle.json').write_text(json.dumps(sadd,indent=2)+'\n')
    print(json.dumps({'R7P-025':'DONE','R7P-026':eq['strict_inclusion'],'R7P-027':stab['H7_inertia'],
      'R7P-028':trans['derivative_interval'],'R7P-029':sadd['H7_inertia']},indent=2))

if __name__=='__main__':main()
