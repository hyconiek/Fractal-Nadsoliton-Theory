from pathlib import Path
import json
import numpy as np
H=Path(__file__).resolve().parent
D=json.loads((H.parent/'fin_transverse_kink_stability_20261113/results.json').read_text())['critical_w']
def run(driven,dx,dt,w,phase):
 eps=.1;T=200;L=230;n=round(2*L/dx);x=np.arange(n)*dx-L;k=2*np.pi*np.fft.fftfreq(n,dx)
 K=np.tanh(x);s=1-K*K;f=np.sqrt(.75)*s;g=np.sqrt(1.5)*np.sqrt(s)*K;gap=w-D;V0=4-6*s+gap
 # e,p,h,q plus two transverse displacement/velocity pairs.
 y=np.zeros((8,n));y[0]=eps*g if driven else 0;y[2]=eps*np.cos(phase)*f if driven else 0;y[3]=-np.sqrt(w)*eps*np.sin(phase)*f if driven else 0;y[4]=f/np.sqrt(gap);y[7]=f
 def rhs(y):
  e,p,h,q,z,v,Z,V=y;lap=np.fft.ifft(-k*k*np.fft.fft(y[[0,2,4,6]],axis=1),axis=1).real
  pot=w-2-D+6*(K+e)**2+3*h*h+5*h*h*h*h
  return np.array([p,lap[0]-(4-6*s)*e-6*K*(e*e+h*h)-6*e*h*h-2*e*e*e,q,lap[1]-(4+w-6*s)*h-12*K*e*h-6*e*e*h-h*h*h-h*h*h*h*h,v,lap[2]-pot*z,V,lap[3]-pot*Z])
 def gram():
  z=y[[4,6]];v=y[[5,7]];zx=np.fft.ifft(1j*k*np.fft.fft(z,axis=1),axis=1).real
  return dx*(v@v.T+zx@zx.T+(z*V0)@z.T)/2
 G0=gram();a,U=np.linalg.eigh(G0);N=(U/np.sqrt(a))@U.T
 work=np.zeros((2,2));hist=[]
 def power(a):
  e,_,h,_,z,v,Z,V=a
  zv=np.array([z*v,Z*V])
  return dx*np.array([zv@(-12*K*e-6*e*e),zv@(-3*h*h-5*h*h*h*h)])
 def obs(t):
  G=gram();values=np.linalg.eigvalsh(N@G@N);hist.append([t,float(values[0]),float(values[1]),float(G[0,0]/G0[0,0]),float(G[1,1]/G0[1,1]),*(work/np.diag(G0)).ravel().tolist()])
 obs(0)
 for j in range(1,round(T/dt)+1):
  a=rhs(y);yb=y+dt*a/2;b=rhs(yb);yc=y+dt*b/2;c=rhs(yc);yd=y+dt*c;d=rhs(yd);work+=dt*(power(y)+2*power(yb)+2*power(yc)+power(yd))/6;y+=dt*(a+2*b+2*c+d)/6
  if j%round(.5/dt)==0:obs(j*dt)
 name=f'{"driven" if driven else "static"}_w{w}_phase{phase:.6f}_dx{dx}_dt{dt}'
 np.savetxt(H/(name+'.csv'),hist,delimiter=',',header='t,min_energy_gain,max_energy_gain,position_gain,velocity_gain,phi_work_position,phi_work_velocity,psi_work_position,psi_work_velocity',comments='')
 return dict(name=name,phase=phase,driven=driven,dx=dx,dt=dt,T=T,L=L,w=w,epsilon=eps if driven else 0,delta=D,gap=gap,max_energy_gain=max(r[2] for r in hist),final_max_energy_gain=hist[-1][2],static_max_deviation=max(abs(r[i]-1) for r in hist for i in [1,2]) if not driven else None)
if __name__=='__main__':
 rows=[]
 for case in [(False,.1,.02,1.25,0.),(True,.1,.02,1.25,np.pi/2),(True,.1,.01,1.25,np.pi/2)]:
  r=run(*case);rows.append(r);print(r,flush=True);(H/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
