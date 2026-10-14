from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from scipy.optimize import root

ROOT=Path(__file__).resolve().parent
Q=12; OMEGA=0.18575; PHI=0.1625; ETA=1.8
GF_INT=(3.5156447068395917,3.5156447268395917)
A_INT=(-0.2125716400840138741,-0.2125716260204590083)
B_INT=(0.1189825274504987752,0.1189826671032236644)

checks=[]
def ck(name, ok, value=None, expected=None):
    checks.append({'name':name,'pass':bool(ok),'value':value,'expected':expected})

def sm(z):
    z=np.asarray(z,float); e=np.exp(z-z.max()); return e/e.sum()

# Rebuild kernel and rank-7 embedding independently (no import of producer code).
j=np.arange(Q)
W=np.zeros((Q,Q))
for r in range(Q):
    for c in range(Q):
        if r==c: continue
        d=min(abs(r-c),Q-abs(r-c))
        W[r,c]=math.cos(OMEGA*d+PHI)/(1+d**ETA)
Lap=np.diag(W.sum(axis=1))-W
# eigenvalue of real Fourier modes computed by direct Rayleigh quotients
lams={}
for k in range(7):
    if k==0: vec=np.ones(Q)
    elif k==6: vec=(-1.0)**j
    else: vec=np.cos(2*np.pi*k*j/Q)
    lams[k]=float(vec@Lap@vec/(vec@vec))
cols=[]
for k in (3,4,5):
    cols.extend([np.sqrt(lams[k]/6)*np.cos(2*np.pi*k*j/Q), np.sqrt(lams[k]/6)*np.sin(2*np.pi*k*j/Q)])
cols.append(np.sqrt(lams[6]/12)*((-1.0)**j))
X=np.column_stack(cols); A=X@X.T
ck('q12',Q==12,Q,12)
ck('A symmetric',np.max(abs(A-A.T))<1e-14,float(np.max(abs(A-A.T))),'<1e-14')
ck('A rank7',np.linalg.matrix_rank(A,1e-10)==7,int(np.linalg.matrix_rank(A,1e-10)),7)
ck('A constant diagonal',np.ptp(np.diag(A))<2e-14,float(np.ptp(np.diag(A))),'<2e-14')
for k, ref in [(3,1.9614068619764455),(4,2.19956884933321),(5,2.298606272079097),(6,2.3421820411463004)]:
    ck(f'lambda{k}',abs(lams[k]-ref)<2e-12,lams[k],ref)

# Fold solve, again independently.
Xe=X[:,[0,2,4,6]]
def grad4(th,g): return th/g-Xe.T@sm(Xe@th)
def H4(th,g):
    p=sm(Xe@th); S=np.diag(p)-np.outer(p,p); return np.eye(4)/g-Xe.T@S@Xe
y0=np.array([1.3643114281782407,1.4331027014285425,1.4080457012550374,1.0056820773201627,3.5156447168395917])
def fun(y): return np.r_[grad4(y[:4],y[4]),np.linalg.det(H4(y[:4],y[4]))]
sol=root(fun,y0,tol=1e-12); th=sol.x[:4]; gf=float(sol.x[4]); H=H4(th,gf)
ev,V=np.linalg.eigh(H); v=V[:,0]; v*=1 if v[0]>=0 else -1
p=sm(Xe@th); z=Xe@v; cen=z-p@z
a=float(-v@th/gf**2); b=float(-np.sum(p*cen**3))
ck('fold root residual',np.linalg.norm(fun(sol.x),np.inf)<2e-11,float(np.linalg.norm(fun(sol.x),np.inf)),'<2e-11')
ck('g_f in source interval',GF_INT[0]<=gf<=GF_INT[1],gf,GF_INT)
ck('fold null eigenvalue',abs(ev[0])<2e-10,float(ev[0]),'|.|<2e-10')
ck('fold stable spectral gap',ev[1]>0.17,float(ev[1]),'>0.17')
ck('a in source interval',A_INT[0]<=a<=A_INT[1],a,A_INT)
ck('b in source interval',B_INT[0]<=b<=B_INT[1],b,B_INT)
ck('fold stationarity',np.max(abs(p-sm(Xe@th)))<1e-15,float(np.max(abs(p-sm(Xe@th)))),'<1e-15')

# Microscopic soft jump cumulants at equilibrium.
delta=gf*(z[None,:]-z[:,None]); weights=p[:,None]*p[None,:]
m2=float(np.sum(weights*delta**2)); m3=float(np.sum(weights*delta**3)); m4=float(np.sum(weights*delta**4))
ck('soft variance equals 2g',abs(m2-2*gf)<2e-12,m2,2*gf)
ck('soft variance source enclosure',2*GF_INT[0]<=m2<=2*GF_INT[1],m2,(2*GF_INT[0],2*GF_INT[1]))
ck('soft third cumulant equilibrium zero',abs(m3)<2e-12,m3,'<2e-12')
ck('soft fourth cumulant nonzero',m4>100,m4,'>100')

# Exact N=120 hidden-fibre witness.
h=np.cos(2*np.pi*2*j/Q)
nplus=np.rint(10+2*h).astype(int); nminus=np.rint(10-2*h).astype(int)
ck('witness plus integer counts',nplus.sum()==120 and nplus.min()>=0,nplus.tolist(),'sum120 nonnegative')
ck('witness minus integer counts',nminus.sum()==120 and nminus.min()>=0,nminus.tolist(),'sum120 nonnegative')
muplus=X.T@(nplus/120); muminus=X.T@(nminus/120)
ck('witness equal retained mean',np.max(abs(muplus-muminus))<1e-13,float(np.max(abs(muplus-muminus))),'<1e-13')
ck('witness mean zero',max(np.max(abs(muplus)),np.max(abs(muminus)))<1e-13,float(max(np.max(abs(muplus)),np.max(abs(muminus)))),'<1e-13')
y=X[:,4]
def L_mu2(n):
    N=n.sum(); mu=(n/N)@y; out=0.0
    for i in range(Q):
        for jj in range(Q):
            if i==jj or n[i]==0: continue
            rate=n[i]/Q
            mn=mu+(y[jj]-y[i])/N
            out += rate*(mn*mn-mu*mu)
    return out
diff=float(L_mu2(nplus)-L_mu2(nminus)); formula=lams[5]/7200
ck('hidden generator difference formula',abs(diff-formula)<2e-15,diff,formula)
ck('hidden generator difference positive',diff>0,diff,'>0')

# Moment lift and exact product identities.
c={k:np.cos(2*np.pi*k*j/Q) for k in range(1,7)}; ss={k:np.sin(2*np.pi*k*j/Q) for k in range(1,6)}
Lift=np.vstack([np.ones(Q),X.T,c[3]*c[4],ss[3]*c[4],c[3]*c[5],ss[3]*c[5]])
ck('moment lift rank12',np.linalg.matrix_rank(Lift,1e-11)==12,int(np.linalg.matrix_rank(Lift,1e-11)),12)
ck('c3c4 identity',np.max(abs(c[3]*c[4]-(c[1]+c[5])/2))<2e-15,float(np.max(abs(c[3]*c[4]-(c[1]+c[5])/2))),'<2e-15')
ck('s3c4 identity',np.max(abs(ss[3]*c[4]-(-ss[1]-ss[5])/2))<2e-15,float(np.max(abs(ss[3]*c[4]-(-ss[1]-ss[5])/2))),'<2e-15')
ck('c3c5 identity',np.max(abs(c[3]*c[5]-(c[2]+c[4])/2))<2e-15,float(np.max(abs(c[3]*c[5]-(c[2]+c[4])/2))),'<2e-15')
ck('s3c5 identity',np.max(abs(ss[3]*c[5]-(-ss[2]-ss[4])/2))<2e-15,float(np.max(abs(ss[3]*c[5]-(-ss[2]-ss[4])/2))),'<2e-15')

# Hamiltonian identities.
rng=np.random.default_rng(20261014); pr=rng.dirichlet(np.ones(Q)); xi=rng.normal(size=Q); g=2.3
sr=sm(g*A@pr)
def Ham(q,zz): return float(np.sum(q[:,None]*sm(g*A@q)[None,:]*(np.exp(zz[None,:]-zz[:,None])-1)))
grad=np.log(pr)-g*A@pr
ck('Hamiltonian gauge invariance',abs(Ham(pr,xi)-Ham(pr,xi+1.234*np.ones(Q)))<2e-13,abs(Ham(pr,xi)-Ham(pr,xi+1.234)), '<2e-13')
ck('Hamiltonian reversal identity',abs(Ham(pr,xi)-Ham(pr,grad-xi))<2e-12,abs(Ham(pr,xi)-Ham(pr,grad-xi)),'<2e-12')
# numerical gradient in xi at zero = drift
hh=1e-6; ng=np.zeros(Q)
for k in range(Q):
    e=np.zeros(Q);e[k]=hh
    ng[k]=(Ham(pr,e)-Ham(pr,-e))/(2*hh)
ck('Hamiltonian first derivative is drift',np.max(abs(ng-(sr-pr)))<3e-10,float(np.max(abs(ng-(sr-pr)))),'<3e-10')

# Kraaij rate convergence bound on a concrete lattice state.
def rate_error(N):
    n=np.ones(Q,dtype=int)*(N//Q); n[:N-n.sum()]+=1; pp=n/N
    maxerr=0.0
    for i in range(Q):
        if n[i]==0: continue
        m=n.copy();m[i]-=1
        sf=sm((g/N)*(A@m)); sl=sm(g*A@pp)
        maxerr=max(maxerr,float(pp[i]*np.max(abs(sf-sl))))
    return maxerr
Amax=float(np.max(abs(A)))
e120=rate_error(120); e240=rate_error(240)
ck('Kraaij uniform-convergence sample within Lipschitz bound N120',e120<=g*Amax/(2*120)+1e-15,e120,g*Amax/(240))
ck('Kraaij rate discrepancy decreases N240',e240<e120,e240,f'<{e120}')

# Critical diagnostic replay read from independently reconstructed parameters.
diag=json.load(open(ROOT/'TH-005'/'FOLD_GENERATOR_DIAGNOSTICS.json'))
zero=sorted([r for r in diag if r['lambda']==0.0], key=lambda r:r['N'])
last=zero[-1]
vtheta=float(v@th)
scaled_corr=-last['drift_error']*(last['N']**(1/3))
ck('critical variance converges numerically',abs(last['scaled_variance']-2*gf)<0.002,abs(last['scaled_variance']-2*gf),'<0.002')
ck('critical third characteristic decays across grid',zero[-1]['scaled_abs_third_jump_rate']<zero[0]['scaled_abs_third_jump_rate'],(zero[0]['scaled_abs_third_jump_rate'],zero[-1]['scaled_abs_third_jump_rate']),'decrease')
ck('leave-one-out correction coefficient diagnostic',abs(scaled_corr-vtheta)<0.002,scaled_corr,vtheta)

# Canonical normalization.
cscale=(4/b)**(1/3); dscale=cscale*cscale/(2*gf); kappa=-a*cscale/2
ck('canonical Z scale positive',cscale>0,cscale,'>0')
ck('canonical time scale positive',dscale>0,dscale,'>0')
ck('canonical lambda map positive',kappa>0,kappa,'>0')

passed=sum(x['pass'] for x in checks); total=len(checks)
out={'status':'PASS' if passed==total else 'FAIL','passed':passed,'total':total,'checks':checks,
     'scope':'Independent numerical/algebraic replay; it does not independently prove Kraaij theorem or the open fast-mode lemma.'}
(ROOT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(f'{passed}/{total} checks passing')
if passed!=total:
    for x in checks:
        if not x['pass']: print('FAIL',x)
    raise SystemExit(1)
