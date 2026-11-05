"""Numerical comparison: exact linear packet drives an odd linear continuum.

This is not a simulation of the prepared nonlinear RP trajectory.
Periodic Fourier space and centered time differences are convergence-tested.
"""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent

def run(R,dt,half_length,dx,m2=16/3,slow_time=1.0):
    n=int(round(2*half_length/dx));dx=2*half_length/n
    x=np.arange(n)*dx-half_length
    k=2*np.pi*np.fft.fftfreq(n,d=dx)
    omega=np.sqrt(k*k+m2)
    dar=np.sqrt((1+k*k)*(4+k*k))
    K=np.tanh(x);s=1/np.cosh(x)**2
    g=np.sqrt(1.5)/np.cosh(x)*K
    g/=np.sqrt(dx*np.dot(g,g))
    reflect=(-np.arange(n))%n
    initial=R**(-0.5)*np.pi**(-0.25)*np.exp(-(x/R)**2/2)
    initial_hat=np.fft.fft(initial)
    coeff=initial_hat/dar
    weight=(1+x*x)**(-3)
    weight_prime=-6*x*(1+x*x)**(-4)
    T=slow_time*R*R;steps=int(round(T/dt));dt=T/steps
    def packet(t):
        G=coeff*np.exp(-1j*omega*t)
        u=(np.fft.ifft(-k*k*G)-3*K*np.fft.ifft(1j*k*G)+(3*K*K-1)*np.fft.ifft(G)).real
        return 0.5*(u+u[reflect])
    def force(t):
        u=packet(t);f=-6*K*u*u
        f=0.5*(f-f[reflect])
        return f-g*(dx*np.dot(g,f))
    packet0=packet(0)
    old=np.zeros(n);cur=0.5*dt*dt*force(0)
    local_integral=0.0;work=0.0;max_bound=0.0
    initial_norm=dx*np.dot(packet0,packet0)
    checkpoints=[]
    for step in range(1,steps+1):
        t=step*dt
        yh=np.fft.fft(cur)
        yxx=np.fft.ifft(-k*k*yh).real
        yx=np.fft.ifft(1j*k*yh).real
        f=force(t)
        nxt=2*cur-old+dt*dt*(yxx-(4-6*s)*cur+f)
        vel=(nxt-old)/(2*dt)
        energy_local=dx*np.sum((weight*cur)**2+(weight*yx+weight_prime*cur)**2+(weight*vel)**2)
        quadrature_weight=0.5 if step==steps else 1.0
        local_integral+=dt*quadrature_weight*energy_local
        work+=dt*quadrature_weight*dx*np.dot(f,vel)
        max_bound=max(max_bound,abs(dx*np.dot(g,cur)))
        if step==steps or step%(max(1,steps//4))==0:
            checkpoints.append({'time':t,'local_phase_norm_squared':float(energy_local)})
        if step==steps:
            energy_global=0.5*dx*np.sum(vel*vel+yx*yx+(4-6*s)*cur*cur)
            edge_mask=np.abs(x)>0.85*half_length
            edge_phase=dx*np.sum((cur[edge_mask])**2+(vel[edge_mask])**2)
        old,cur=cur,nxt
    return {'R':R,'dt':dt,'dx':dx,'half_length':half_length,'points':n,'T':T,
            'integrated_local_phase_squared':float(local_integral),
            'final_global_energy':float(energy_global),'integrated_source_work':float(work),
            'energy_work_relative_defect':float(abs(energy_global-work)/max(abs(energy_global),abs(work),1e-15)),
            'initial_packet_L2_squared':float(initial_norm),'max_shape_projection':float(max_bound),
            'final_outer_band_phase_squared':float(edge_phase),'checkpoints':checkpoints}

if __name__=='__main__':
    configurations=[(4,.02,64,.0625),(4,.01,64,.0625),(4,.01,64,.03125),
                    (4,.01,96,.0625),(8,.02,128,.0625),(8,.01,128,.0625)]
    results=[]
    for cfg in configurations:
        r=run(*cfg);results.append(r)
        print(json.dumps({k:r[k] for k in ['R','dt','dx','half_length','integrated_local_phase_squared','energy_work_relative_defect']}),flush=True)
    out={'scope':'Finite periodic linear packet-driven comparison, not nonlinear RP evolution',
         'numpy_version':np.__version__,'mass_squared':'16/3','slow_time':1.0,'results':results}
    (HERE/'cascade_preview.json').write_text(json.dumps(out,indent=2)+'\n')
