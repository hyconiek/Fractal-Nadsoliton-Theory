"""Full nonlinear RP evolution, fourth-order finite differences and RK4.

The solver does not use the averaged ODE as an evolution equation.
Fixed g4=g6=1 is a declared example, never a fit to the damping law.
"""
from pathlib import Path
import json,math
import numpy as np
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'fin_canonical_audit_pde_pilot_20261107'))
from canonical import pullback
HERE=Path(__file__).resolve().parent

def d2(u,dx):
    return (-np.roll(u,2)+16*np.roll(u,1)-30*u+16*np.roll(u,-1)-np.roll(u,-2))/(12*dx*dx)

def ka(w):
    k=np.sqrt(4*w-4);S2=np.pi*k/np.sinh(np.pi*k/2);S4=(4+k*k)*S2/6
    M=3j*k*(k*k+4)*S4/(32*(1+1j*k)*(2+1j*k))
    return float(36*abs(M)**2/(2*k*w))

def run(eps,dx,dt,extra=0,w=1.5):
    g4=g6=1.;T=1/eps**2;L=math.ceil(T+30+extra)
    n=2*math.ceil(L/dx);dx=2*L/n;x=np.arange(n)*dx-L
    K=np.tanh(x);sech=2*np.exp(-abs(x))/(1+np.exp(-2*abs(x)))
    s=sech**2;f=np.sqrt(.75)*s;g=np.sqrt(1.5)*sech*K
    omega=np.sqrt(w);lam=np.sqrt(3.)
    steps=math.ceil(T/dt);dt=T/steps
    y=np.zeros((4,n));y[2]=eps*f
    def rhs(y):
        e,ep,h,hp=y
        return np.array([ep,d2(e,dx)-(4-6*s)*e-6*K*(e*e+h*h)-6*e*h*h-2*e**3,
                         hp,d2(h,dx)-(4+w-6*s)*h-12*K*e*h-6*e*e*h-g4*h**3-g6*h**5])
    def observe(y,t):
        e,ep,h,hp=y
        z,p,c,r=[dx*np.dot(a,b) for a,b in [(h,f),(hp,f),(e,g),(ep,g)]]
        Z,P,C,R=pullback([z,p,c,r],w,g4)
        I=(P*P+w*Z*Z)/(2*omega);J=(R*R+3*C*C)/(2*lam)
        v=e-c*g;vm=ep-r*g;u=h-z*f;um=hp-p*f
        cont=.5*dx*np.sum(vm*vm+um*um-v*d2(v,dx)+(4-6*s)*v*v-u*d2(u,dx)+(4+w-6*s)*u*u)
        E=.5*dx*np.sum(ep*ep+hp*hp-e*d2(e,dx)+(4-6*s)*e*e-h*d2(h,dx)+(4+w-6*s)*h*h)
        E+=dx*np.sum(2*K*e**3+6*K*e*h*h+.5*e**4+3*e*e*h*h+g4*h**4/4+g6*h**6/6)
        reflection=(-np.arange(n))%n
        parity=max(float(np.max(abs(e+e[reflection]))),float(np.max(abs(h-h[reflection]))))
        return {'t':float(t),'I':float(I),'J':float(J),'total_energy':float(E),'continuum_quadratic_energy':float(cont),
                'interaction_and_chart_energy':float(E-omega*I-lam*J-cont),'parity_error':parity}
    history=[observe(y,0)];stride=max(1,steps//200)
    for step in range(1,steps+1):
        a=rhs(y);b=rhs(y+dt*a/2);c=rhs(y+dt*b/2);d=rhs(y+dt*c)
        y+=dt*(a+2*b+2*c+d)/6
        if step%stride==0 or step==steps:history.append(observe(y,step*dt))
    I0=history[0]['I'];E0=history[0]['total_energy'];k=ka(w)
    for row in history:row['leading_ode_I']=I0/(1+k*I0*row['t'])
    tail=[v for v in history if v['t']>=T/4]
    slope=float(np.polyfit([v['t'] for v in tail],[1/v['I'] for v in tail],1)[0])
    return {'epsilon':eps,'w':w,'g4':g4,'g6':g6,'dx':dx,'dt':dt,'L':L,'T':T,'points':n,
      'predicted_ka':k,'finite_window_inverse_action_slope':slope,'slope_over_ka':slope/k,
      'relative_energy_drift':max(abs(v['total_energy']-E0) for v in history)/E0,
      'max_parity_error':max(v['parity_error'] for v in history),
      'max_normalized_action_error':max(abs(v['I']-v['leading_ode_I']) for v in history)/eps**2,
      'history':history}


if __name__=='__main__':
    results=[]
    for w in [1.25,4/3-.01,4/3,4/3+.01,1.75]:
        for eps in [.2,.15,.1]:
            for dx,dt in [(.2,.04),(.1,.02)]:
                v=run(eps,dx,dt,w=w)
                h=v.pop('history')
                v['initial_I']=h[0]['I'];v['final_I']=h[-1]['I']
                v['continuum_energy_gain']=h[-1]['continuum_quadratic_energy']-h[0]['continuum_quadratic_energy']
                v['samples']=h[::20]+([h[-1]] if h[-1] not in h[::20] else [])
                results.append(v)
                (HERE/'results.json').write_text(json.dumps({'scope':'Finite-window parameter campaign, not empirical validation','runs':results},indent=2)+'\n')
                print(json.dumps({k:v[k] for k in ['w','epsilon','dx','slope_over_ka','relative_energy_drift']}),flush=True)
