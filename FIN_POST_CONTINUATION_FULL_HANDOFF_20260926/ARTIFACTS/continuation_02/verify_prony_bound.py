import numpy as np

def moments(a,l):
    return np.array([sum(ai*li**n for ai,li in zip(a,l)) for n in range(4)],float)

def recover(m):
    A=np.array([[m[1],-m[0]],[m[2],-m[1]]],float)
    b=np.array([m[2],m[3]],float)
    s=np.linalg.solve(A,b)
    return np.sort(np.roots([1,-s[0],s[1]]).real),A,s

def certificate(m,eps):
    roots,A,s=recover(m)
    sig=np.linalg.svd(A,compute_uv=False)[-1]
    if sig<=2*eps:return dict(ok=False,sig=sig,reason='linear system uncertified')
    eta=(np.sqrt(2)*eps+2*eps*np.linalg.norm(s))/(sig-2*eps)
    delta=roots[1]-roots[0]
    r=2*abs(s[0])*eta+eta*eta+4*eta
    if r>=delta*delta:return dict(ok=False,sig=sig,eta_s=eta,delta=delta,r=r,reason='root gap uncertified')
    eroot=.5*(eta+r/delta)
    D=m[0]*m[2]-m[1]**2
    fro=np.linalg.norm(A,'fro')
    return dict(ok=True,sig=sig,eta_s=eta,delta=delta,r=r,root_bound=eroot,D=D,sig_lower_D=D/fro)

def trial(a,l,rel=1e-8):
    m=moments(a,l); eps=rel*np.max(abs(m)); c=certificate(m,eps)
    # fixed adversarial sign pattern only as reproducible witness, not worst-case proof
    mp=m+eps*np.array([1,-1,1,-1.])
    rr,_,_=recover(mp); actual=np.max(np.abs(rr-np.sort(l)))
    return m,eps,c,actual

for title,cases in [
 ('near_poles',[([1,1],[1,2]),([1,1],[1,1.2]),([1,1],[1,1.05]),([1,1],[1,1.01])]),
 ('weak_residue',[([1,1],[1,2]),([1,.2],[1,2]),([1,.05],[1,2]),([1,.01],[1,2]),([1,.001],[1,2])])]:
    print(title)
    for a,l in cases:
        m,eps,c,actual=trial(a,l)
        print('a',a,'lam',l,'eps',eps,'D2',m[0]*m[2]-m[1]**2,'actual',actual,'cert',c)
        if c.get('ok'): assert actual <= c['root_bound']*(1+1e-9)
