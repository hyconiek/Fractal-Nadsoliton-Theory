"""Triangular packet/odd-response/even-feedback comparison, not full nonlinear RP."""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent

def run(R,dt,L,dx,m2=16/3,S=1.0):
 n=int(round(2*L/dx));dx=2*L/n;x=np.arange(n)*dx-L
 k=2*np.pi*np.fft.fftfreq(n,d=dx);om=np.sqrt(k*k+m2)
 dar=np.sqrt((1+k*k)*(4+k*k));K=np.tanh(x);ss=1/np.cosh(x)**2
 g=np.sqrt(1.5)/np.cosh(x)*K;g/=np.sqrt(dx*np.dot(g,g))
 f=np.sqrt(.75)*ss;f/=np.sqrt(dx*np.dot(f,f))
 reflect=(-np.arange(n))%n
 initial=R**(-.5)*np.pi**(-.25)*np.exp(-(x/R)**2/2)
 coeff=np.fft.fft(initial)/dar
 T=S*R*R;steps=int(round(T/dt));dt=T/steps
 if dt*np.sqrt(k.max()**2+max(4,m2))>=2:raise ValueError('Time step violates spectral wave stability bound')
 def packet(t):
  G=coeff*np.exp(-1j*om*t)
  out=(np.fft.ifft(-k*k*G)-3*K*np.fft.ifft(1j*k*G)+(3*K*K-1)*np.fft.ifft(G)).real
  return .5*(out+out[reflect])
 def oddforce(p):
  F=-6*K*p*p;F=.5*(F-F[reflect]);cg=dx*np.dot(g,F)
  return F-g*cg,cg
 def evenforce(p,y,c):
  F=np.array([-12*K*p*y,-12*K*p*c*g,-p**3])
  F=.5*(F+F[:,reflect]);return F-dx*np.outer(F@f,f)
 p0=packet(0);fy0,fc0=oddforce(p0)
 yo=np.zeros(n);yc=.5*dt*dt*fy0;co=0.;cc=.5*dt*dt*fc0
 zo=np.zeros((3,n));zc=.5*dt*dt*evenforce(p0,yo,co)
 test=np.zeros(n);inside=np.abs(x)<1+1e-12;test[inside]=(1+x[inside]**2)**(-.5)
 test[np.isclose(np.abs(x),1,atol=1e-12)]*=.5
 gram=np.zeros((3,3));work=np.zeros(3);maxf=np.zeros(3);maxg=0.;maxc=0.
 for j in range(1,steps+1):
  t=j*dt;p=packet(t);fy,fc=oddforce(p);fz=evenforce(p,yc,cc)
  yh=np.fft.fft(yc);yn=2*yc-yo+dt*dt*(np.fft.ifft(-k*k*yh).real-(4-6*ss)*yc+fy)
  cn=2*cc-co+dt*dt*(-3*cc+fc)
  zh=np.fft.fft(zc,axis=-1);zx=np.fft.ifft(1j*k*zh,axis=-1).real;zxx=np.fft.ifft(-k*k*zh,axis=-1).real
  zn=2*zc-zo+dt*dt*(zxx-(m2-6*ss)*zc+fz)
  zv=(zn-zo)/(2*dt)
  wd=np.fft.ifft(np.fft.fft(zxx+3*K*zx+2*zc,axis=-1)/dar,axis=-1).real
  q=dx*(wd@test);wt=.5 if j==steps else 1.
  gram+=wt*dt*np.outer(q,q);work+=wt*dt*dx*np.sum(fz*zv,axis=-1)
  maxf=np.maximum(maxf,np.abs(dx*(zc@f)));maxg=max(maxg,abs(dx*np.dot(g,yc)));maxc=max(maxc,abs(cc))
  if j==steps:
   energies=.5*dx*np.sum(zv*zv+zx*zx+(m2-6*ss)*zc*zc,axis=-1)
   edge=np.abs(x)>.85*L
   edge_norms=dx*np.sum(zc[:,edge]**2+zv[:,edge]**2,axis=-1)
  yo,yc=yc,yn;co,cc=cc,cn;zo,zc=zc,zn
 return {'R':R,'dt':dt,'dx':dx,'half_domain':L,'T':T,'q_gram':gram.tolist(),
         'q_total_g4_one_squared_integral':float(np.sum(gram)),
         'final_even_energies':energies.tolist(),'integrated_even_source_work':work.tolist(),
         'relative_energy_work_defects':(np.abs(energies-work)/np.maximum(np.maximum(np.abs(energies),np.abs(work)),1e-15)).tolist(),
         'max_even_bound_projections':maxf.tolist(),'max_odd_shape_projection':float(maxg),
         'max_shape_amplitude':float(maxc),'final_outer_band_phase_squared':edge_norms.tolist()}

if __name__=='__main__':
 results=[]
 for cfg in [(4,.02,64,.0625),(4,.01,64,.0625),(4,.01,64,.03125),(4,.01,96,.0625),(8,.02,128,.0625),(8,.01,128,.0625)]:
  r=run(*cfg);results.append(r)
  print(json.dumps({k:r[k] for k in ['R','dt','dx','half_domain','q_total_g4_one_squared_integral','relative_energy_work_defects']}),flush=True)
 (HERE/'feedback_preview.json').write_text(json.dumps({'scope':'Triangular linear comparison with real packet, odd continuum and shape, then three even responses; not full nonlinear RP','numpy_version':np.__version__,'mass_squared':'16/3','S':1,'channel_order':['odd_continuum_feedback','shape_feedback','direct_cubic_g4_one'],'results':results},indent=2)+'\n')
