"""Fourier collocation replication with physical interior energy-flux balance."""
from pathlib import Path
import sys,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'fin_canonical_audit_pde_pilot_20261107'))
from canonical import pullback
HERE=Path(__file__).resolve().parent

def run(w,dx,dt,eps,preparation):
    L=130;T=100;n=round(2*L/dx);x=np.arange(n)*dx-L
    k=2*np.pi*np.fft.fftfreq(n,d=dx);K=np.tanh(x);s=1-K*K
    f=np.sqrt(.75)*s;g=np.sqrt(1.5)*np.sqrt(s)*K
    regions=[np.where(abs(x)<=R+1e-9)[0] for R in [10,20]]
    y=np.zeros((4,n));y[0]=eps*g;y[2]=(eps*f if preparation=="mixed" else 0);flux=np.zeros(2)
    def calc(y):
        e,p,h,q=y;eh=np.fft.fft(e);hh=np.fft.fft(h)
        ex=np.fft.ifft(1j*k*eh).real;hx=np.fft.ifft(1j*k*hh).real
        ee=np.fft.ifft(-k*k*eh).real;hh2=np.fft.ifft(-k*k*hh).real
        F=np.array([p,ee-(4-6*s)*e-6*K*(e*e+h*h)-6*e*h*h-2*e**3,
                    q,hh2-(4+w-6*s)*h-12*K*e*h-6*e*e*h-h**3-h**5])
        j=-p*(s+ex)-q*hx
        return F,np.array([j[idx[-1]]-j[idx[0]] for idx in regions]),ex,hx
    def obs(t):
        _,_,ex,hx=calc(y);e,p,h,q=y
        # Exact polynomial expansion of V(K+e,h)-V(K,0), retaining linear terms.
        density=.5*(p*p+q*q+ex*ex+hx*hx)+s*ex+2*K*(K*K-1)*e+(3*K*K-1)*e*e+2*K*e**3+.5*e**4
        density+=(w-2+6*K*K)*h*h/2+6*K*e*h*h+3*e*e*h*h+h**4/4+h**6/6
        local=[]
        for idx in regions:
            a=density[idx];assert len(a)%2==1
            local.append(float(dx*(a[0]+a[-1]+4*a[1:-1:2].sum()+2*a[2:-1:2].sum())/3))
        z,pp,c,r=[dx*np.dot(a,b) for a,b in [(h,f),(q,f),(e,g),(p,g)]]
        Z,P,C,R=pullback([z,pp,c,r],w,1.)
        return {'t':t,'I':float((P*P+w*Z*Z)/(2*np.sqrt(w))),'J':float((R*R+3*C*C)/(2*np.sqrt(3))),'global_energy':float(dx*density.sum()),'local_energy':local,'integrated_outward_flux':flux.tolist()}
    hist=[obs(0)];steps=round(T/dt)
    for j in range(1,steps+1):
        a,fa,_,_=calc(y);b,fb,_,_=calc(y+dt*a/2);c,fc,_,_=calc(y+dt*b/2);d,fd,_,_=calc(y+dt*c)
        y+=dt*(a+2*b+2*c+d)/6;flux+=dt*(fa+2*fb+2*fc+fd)/6
        if j%round(.5/dt)==0:hist.append(obs(j*dt))
    E0=hist[0]['global_energy'];tail=[r for r in hist if r['t']>=25]
    residual=[max(abs(r['local_energy'][i]-hist[0]['local_energy'][i]+r['integrated_outward_flux'][i]) for r in hist)/E0 for i in range(2)]
    return {'w':w,'epsilon':eps,'g4':1,'g6':1,'dx':dx,'dt':dt,'L':L,'T':T,'radii':[10,20],'preparation':preparation,
    'relative_global_energy_drift':max(abs(r['global_energy']-E0) for r in hist)/E0,
    'max_relative_local_balance_residual':residual,'samples':hist}

if __name__=='__main__':
    runs=[]
    for preparation in ['shape','mixed']:
        for eps,dx,dt in [(.15,.1,.02),(.1,.1,.02),(.1,.05,.01)]:
            r=run(1.5,dx,dt,eps,preparation);runs.append(r)
            (HERE/'results.json').write_text(json.dumps({'scope':'Auxiliary preparations outside historical AW theorem','runs':runs},indent=2)+'\n')
            print(preparation,eps,dx,r['relative_global_energy_drift'],flush=True)
