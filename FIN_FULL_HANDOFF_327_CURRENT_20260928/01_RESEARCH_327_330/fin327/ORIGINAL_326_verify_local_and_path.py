from pathlib import Path
import json, math
import numpy as np
import mpmath as mp
from scipy.optimize import root
OUT=Path('/mnt/data/fin326')
mp.mp.dps=90
mp.iv.dps=70
iv=mp.iv
Q=12; Gf=5.145228719489142
# computational model constants, reconstructed exactly as in 325 then frozen as decimal singletons
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),Q-abs(i-j))+0.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)],float)
L=np.diag(W.sum(1))-W; lam=np.fft.fft(L[0]).real[:7]; jj=np.arange(Q); cols=[]
for k in (3,4,5): cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
Xf=np.column_stack(cols); Af=Xf@Xf.T
# exact tangent basis used in 325 (normal in k4 plane is (sqrt3/2,-1/2))
Tc=np.zeros((7,6)); Tc[0,0]=1;Tc[1,1]=1;Tc[2,2]=.5;Tc[3,2]=math.sqrt(3)/2;Tc[4,3]=1;Tc[5,4]=1;Tc[6,5]=1
Bf=Xf@Tc
G=mp.mpf(str(Gf)); X=[[mp.mpf(str(x)) for x in row] for row in Xf]; A=[[mp.mpf(str(x)) for x in row] for row in Af]; B=[[mp.mpf(str(x)) for x in row] for row in Bf]
c=[mp.mpf(1 if j%3==0 else -1 if j%3==1 else 0) for j in range(Q)]
classes=[[j for j in range(Q) if j%3==a] for a in range(3)]

def soft_mp(v):
    m=max(v); e=[mp.exp(x-m) for x in v]; s=sum(e); return [x/s for x in e]
def p_y(y):
    h=[sum(B[j][k]*y[k] for k in range(6)) for j in range(Q)]
    S=[sum(mp.exp(h[j]) for j in classes[a]) for a in range(3)]
    z=mp.mpf('.5')*mp.log(S[1]/S[0])
    return soft_mp([h[j]+z*c[j] for j in range(Q)])
def F_y(*yy):
    y=list(yy);p=p_y(y);return tuple(y[k]/G-sum(B[j][k]*p[j] for j in range(Q)) for k in range(6))
def jac_point(y):
    p=p_y(y)
    C=[[ (p[i] if i==j else mp.mpf('0'))-p[i]*p[j] for j in range(Q)] for i in range(Q)]
    def bil(a,b): return sum(a[i]*C[i][j]*b[j] for i in range(Q) for j in range(Q))
    K=mp.matrix(6); u=[bil([B[i][a] for i in range(Q)],c) for a in range(6)]; cc=bil(c,c)
    for a in range(6):
        ba=[B[i][a] for i in range(Q)]
        for b in range(6):
            bb=[B[i][b] for i in range(Q)]
            K[a,b]=bil(ba,bb)-u[a]*u[b]/cc
    J=mp.eye(6)/G-K
    return J
# solve representative d4 constrained root
ystart=[mp.mpf(str(x)) for x in [2.57030486,0,1.29925167,.681792653,-1.18089952,1.97598638]]
yroot=list(mp.findroot(F_y, tuple(ystart), tol=mp.mpf('1e-70'), maxsteps=100))
J0=jac_point(yroot); R=J0**-1
# interval helpers
def I(x):
    if isinstance(x,(tuple,list)): return iv.mpf([mp.nstr(x[0],100),mp.nstr(x[1],100)])
    return iv.mpf(mp.nstr(x,100))
def bounds(q):
    s=mp.nstr(q,90)
    if s.startswith('['):
        a=s[1:s.index(',')]; b=s[s.index(',')+1:s.rindex(']')].strip(); return mp.mpf(a),mp.mpf(b)
    v=mp.mpf(s);return v,v
Biv=[[I(B[j][k]) for k in range(6)] for j in range(Q)]; civ=[I(x) for x in c]

def p_iv(Y):
    h=[sum(Biv[j][k]*Y[k] for k in range(6)) for j in range(Q)]
    S=[sum((iv.exp(h[j]) for j in classes[a]),I(0)) for a in range(3)]
    z=I(mp.mpf('.5'))*iv.log(S[1]/S[0])
    e=[iv.exp(h[j]+z*civ[j]) for j in range(Q)]; den=sum(e,I(0)); return [x/den for x in e]
def jac_iv(Y):
    p=p_iv(Y)
    C=[[ (p[i] if i==j else I(0))-p[i]*p[j] for j in range(Q)] for i in range(Q)]
    def bil(a,b): return sum((a[i]*C[i][j]*b[j] for i in range(Q) for j in range(Q)),I(0))
    u=[]
    for a in range(6):u.append(bil([Biv[i][a] for i in range(Q)],civ))
    cc=bil(civ,civ)
    J=[[I(0) for _ in range(6)] for __ in range(6)]
    for a in range(6):
        ba=[Biv[i][a] for i in range(Q)]
        for b in range(6):
            bb=[Biv[i][b] for i in range(Q)]
            K=bil(ba,bb)-u[a]*u[b]/cc
            J[a][b]=(I(1)/I(G) if a==b else I(0))-K
    return J
# Krawczyk at tiny radius
krad=mp.mpf('1e-18'); Y=[I((yroot[k]-krad,yroot[k]+krad)) for k in range(6)]
JX=jac_iv(Y); F0=[I(F_y(*yroot)[k]) for k in range(6)]
K0=[]
for i in range(6): K0.append(I(yroot[i])-sum((I(R[i,j])*F0[j] for j in range(6)),I(0)))
E=[[ (I(1) if i==j else I(0))-sum((I(R[i,k])*JX[k][j] for k in range(6)),I(0)) for j in range(6)] for i in range(6)]
dx=[I((-krad,krad)) for _ in range(6)]
KI=[K0[i]+sum((E[i][j]*dx[j] for j in range(6)),I(0)) for i in range(6)]
kraw=[];strict=True
for i in range(6):
    xl,xu=bounds(Y[i]);kl,ku=bounds(KI[i]); ok=(xl<kl and ku<xu);strict &= ok;kraw.append({'X':[str(xl),str(xu)],'K':[str(kl),str(ku)],'strict':bool(ok)})
# interval Hessian on local cube radius .02
hrad=mp.mpf('0.02'); Yh=[I((yroot[k]-hrad,yroot[k]+hrad)) for k in range(6)]; JH=jac_iv(Yh)
HL=np.zeros((6,6));HU=np.zeros((6,6))
for i in range(6):
    for j in range(6):
        a,b=bounds(JH[i][j]);HL[i,j]=float(a);HU[i,j]=float(b)
mid=(HL+HU)/2; rad=(HU-HL)/2
hess_lb=float(np.linalg.eigvalsh((mid+mid.T)/2)[0]-np.linalg.norm(rad,'fro')-1e-12)
# full stationary loc/sad in 7D theta for path certificate
def p_theta(th): return soft_mp([sum(X[j][a]*th[a] for a in range(7)) for j in range(Q)])
def Ft(*tt):
    th=list(tt);p=p_theta(th);return tuple(th[a]-G*sum(X[j][a]*p[j] for j in range(Q)) for a in range(7))
# starts from double model
loc0=np.array([.980692347,4.4082e-5,8.40916e-4,3.792984e-3,3.351279e-3,1.249181e-3,7.507681e-4,1.249181e-3,3.351279e-3,3.792984e-3,8.40916e-4,4.4082e-5])
sad0=np.array([.43331807,.0054905,.01074374,.0054905,.43331807,.00350416,.00985837,.0236854,.03754326,.0236854,.00985837,.00350416])
loc_th=[G*sum(X[j][a]*mp.mpf(str(loc0[j])) for j in range(Q)) for a in range(7)]
sad_th=[G*sum(X[j][a]*mp.mpf(str(sad0[j])) for j in range(Q)) for a in range(7)]
loc_th=list(mp.findroot(Ft,tuple(loc_th),tol=mp.mpf('1e-65'),maxsteps=100));sad_th=list(mp.findroot(Ft,tuple(sad_th),tol=mp.mpf('1e-65'),maxsteps=100))
ploc=p_theta(loc_th); psad=p_theta(sad_th); ploc4=ploc[-4:]+ploc[:-4]  # np.roll +4 equivalent new[j]=old[j-4]
# fix roll exact
ploc4=[ploc[(j-4)%12] for j in range(12)]
def V_mp(p): return sum(p[j]*mp.log(12*p[j]) for j in range(Q))-G*sum(p[i]*A[i][j]*p[j] for i in range(Q) for j in range(Q))/2
Vloc=V_mp(ploc);Vsad=V_mp(psad)
# interval path derivative certificates
def path_cert(a,b,sign=1,nmid=2000,edge=mp.mpf('0.002')):
    d=[b[i]-a[i] for i in range(Q)]
    dAd=sum(d[i]*A[i][j]*d[j] for i in range(Q) for j in range(Q))
    # second derivative on endpoint interval
    def d2_iv(tI):
        p=[I(a[i])+tI*I(d[i]) for i in range(Q)]
        return sum((I(d[i]*d[i])/p[i] for i in range(Q)),I(0))-I(G*dAd)
    left=d2_iv(I((0,edge))); right=d2_iv(I((1-edge,1)))
    ll,lu=bounds(left);rl,ru=bounds(right)
    # direct derivative middle intervals
    worst=mp.inf if sign>0 else -mp.inf; ok=True
    for r in range(nmid):
        aa=edge+(1-2*edge)*r/nmid;bb=edge+(1-2*edge)*(r+1)/nmid;tI=I((aa,bb))
        p=[I(a[i])+tI*I(d[i]) for i in range(Q)]
        grad=[]
        for i in range(Q):
            Ap=sum((I(A[i][j])*p[j] for j in range(Q)),I(0))
            grad.append(iv.log(I(12)*p[i])+I(1)-I(G)*Ap)
        dv=sum((I(d[i])*grad[i] for i in range(Q)),I(0));dl,du=bounds(dv)
        if sign>0: worst=min(worst,dl); ok &= dl>0
        else: worst=max(worst,du); ok &= du<0
    # endpoint logic: first segment sign + requires left d2>0 and right d2<0; reverse sign - has left d2<0 and right d2>0
    edgeok=(ll>0 and ru<0) if sign>0 else (lu<0 and rl>0)
    return {'middle_strict':bool(ok),'endpoint_curvature_strict':bool(edgeok),'middle_worst_bound':str(worst),'left_d2':[str(ll),str(lu)],'right_d2':[str(rl),str(ru)]}
path1=path_cert(ploc,psad,+1);path2=path_cert(psad,ploc4,-1)
# energies and output
out={'task':326,'root_y':[mp.nstr(x,50) for x in yroot],'krawczyk_radius':str(krad),'krawczyk_strict':bool(strict),'krawczyk':kraw,
     'local_hessian_cube_radius':'0.02','local_hessian_eigen_lower_bound':hess_lb,
     'localized_energy':mp.nstr(Vloc,50),'d4_energy':mp.nstr(Vsad,50),'barrier':mp.nstr(Vsad-Vloc,50),
     'path_loc_to_d4':path1,'path_d4_to_shift4loc':path2}
(OUT/'LOCAL_PATH_INTERVAL_CERT_326.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
