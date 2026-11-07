"""Exact causal identities and finite-grid threshold profile controls."""
from pathlib import Path
import hashlib,json,math
import sympy as sp
import numpy as np
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,condition):
    checks.append({'name':name,'pass':bool(condition)})
    if not condition: raise AssertionError(name)
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
    check('input: '+row['path'],hashlib.sha256((HERE.parent/row['path']).read_bytes()).hexdigest()==row['sha256'])
z=sp.symbols('z',real=True)
moment=3/(2*sp.sqrt(2))*(3*sp.pi/16-sp.pi/8)/2
check('mixed profile threshold moment',sp.simplify(moment-3*sp.pi/(64*sp.sqrt(2)))==0)
check('threshold orthogonal to bound f',sp.integrate((3*z*z-1)/2,(z,-1,1))==0)
w=sp.symbols('w',positive=True)
lam=sp.sqrt(3);omega=sp.sqrt(w);m=sp.sqrt(4+w)
check('sum squared gap',sp.simplify((lam+omega)**2-m*m-(-1+2*sp.sqrt(3*w)))==0)
check('difference squared gap',sp.simplify(m*m-(lam-omega)**2-(1+2*sp.sqrt(3*w)))==0)
check('uniform rational gap certificates',sp.sqrt(15)>3 and sp.sqrt(3)+sp.sqrt(sp.Rational(7,4))+sp.sqrt(sp.Rational(23,4))<6)
# Abel-regularized Laplace transforms of sine and cosine.
a,b,nu=sp.symbols('a b nu',positive=True,real=True)
pole=a-sp.I*nu
position=1/(pole*pole+b*b);velocity=pole/(pole*pole+b*b)
check('causal susceptibility sign',sp.simplify(position-1/(b*b-(nu+sp.I*a)**2))==0)
check('causal momentum factor',sp.simplify(velocity+sp.I*(nu+sp.I*a)*position)==0)
t,r=sp.symbols('t r',nonnegative=True)
# Integrable example with nonzero A(0) and A' to test all endpoint signs.
A=lambda v: 1+v
P=sp.Rational(1,2);E=lambda v:sp.exp(-2*v)/2
lhs=sp.integrate(sp.exp(-2*r)*A(t-r),(r,0,t))
rhs=P*A(t)-E(t)*A(0)-sp.integrate(E(r),(r,0,t))
check('finite retarded identity all terms',sp.simplify(lhs-rhs)==0)
check('omitting switch-on violates initial data',(P*A(0))!=0 and sp.simplify(rhs.subs(t,0))==0)
check('tail L2 squared logarithm',sp.integrate(1/(1+r),(r,0,t))==sp.log(1+t))


def sample(n,L,w,t):
    x=np.arange(n)*2*L/n-L;k=2*np.pi*np.fft.fftfreq(n,d=2*L/n)
    K=np.tanh(x);sech=2*np.exp(-abs(x))/(1+np.exp(-2*abs(x)))
    f=np.sqrt(.75)*sech**2;g=np.sqrt(1.5)*sech*K
    B=K*f*g;B-=f*np.sum(f*B)*2*L/n
    bh=np.fft.fft(B);root=np.sqrt((k*k+1)*(k*k+4))
    tr=np.fft.fft(np.fft.ifft(-k*k*bh)+3*K*np.fft.ifft(1j*k*bh)+2*B)/root
    m=np.sqrt(4+w);freq=np.sqrt(k*k+m*m)
    def W(h):
        q=h/root
        return np.fft.ifft(-k*k*q)-3*K*np.fft.ifft(1j*k*q)+(3*K*K-1)*np.fft.ifft(q)
    pos=W(tr*np.sin(t*freq)/freq);vel=W(tr*np.cos(t*freq))
    r0=(3*K*K-1)/2;mu=3*np.pi/(64*np.sqrt(2));theta=m*t+np.pi/4
    leadpos=r0*mu*np.sin(theta)/np.sqrt(2*np.pi*m*t)
    leadvel=r0*mu*m*np.cos(theta)/np.sqrt(2*np.pi*m*t)
    weight=(1+x*x)**-3;dx=2*L/n
    def norm(q,p):
        q=weight*q;dq=np.fft.ifft(1j*k*np.fft.fft(q))
        return float(np.sqrt(dx*np.sum(abs(q)**2+abs(dq)**2))+np.sqrt(dx*np.sum(abs(weight*p)**2)))
    err=norm(pos-leadpos,vel-leadvel)
    return {'w':w,'t':t,'phase_norm':norm(pos,vel),'remainder_norm':err,
            'scaled_remainder':t**1.5*err,
            'threshold_moment_error':float(abs(np.sum(r0*B)*dx-mu)),
            'projection_error':float(abs(np.sum(f*B)*dx))}
rows=[sample(16384,512.,w,t) for w in [1.25,1.5,1.75] for t in [16.,32.,64.,128.,256.]]
check('threshold sample moments',max(v['threshold_moment_error'] for v in rows)<1e-13)
check('continuous source sample projection',max(v['projection_error'] for v in rows)<1e-13)
check('finite grid diagnostics finite',all(np.isfinite(v['scaled_remainder']) for v in rows))
base=next(v for v in rows if v['w']==1.5 and v['t']==256.)
fine=sample(32768,512.,1.5,256.);wide=sample(24576,768.,1.5,256.)
refine={kind:max(abs(v[k]-base[k]) for k in ['phase_norm','remainder_norm']) for kind,v in [('grid',fine),('domain',wide)]}
check('grid refinement consistency',refine['grid']<1e-8)
check('domain refinement consistency',refine['domain']<1e-8)
out={'status':'PASS','scope':'Exact identities and finite-grid linear spectral diagnostics; no nonlinear PDE validation or independent theorem certification.',
'check_count':len(checks),'checks':checks,'samples':rows,'refinement_norm_differences':refine,
'versions':{'numpy':np.__version__,'sympy':sp.__version__}}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks),'max_scaled_remainder':max(v['scaled_remainder'] for v in rows),'refinements':refine}))
