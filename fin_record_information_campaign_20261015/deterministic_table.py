"""Deterministic finite-difference / finite-volume diagnostic of the fixed RI experiment.
Not a certified PDE solver: do not use its printed decimal values as interval bounds.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.ndimage import gaussian_filter1d, map_coordinates
import json,sys


def run(hz=0.10, hy=0.10, dt=0.025, Lz=4.5, Ly=3.8):
    nz=int(round(2*Lz/hz))+1; ny=int(round(2*Ly/hy))+1
    z=np.linspace(-Lz,Lz,nz); y=np.linspace(-Ly,Ly,ny)
    hz=z[1]-z[0]; hy=y[1]-y[0]
    logrho=-z**4/12
    rp=np.exp(np.diff(logrho)/2)/hz**2
    rm=np.exp(-np.diff(logrho)/2)/hz**2
    Q=sparse.diags([rm, -(np.r_[rp,0]+np.r_[0,rm]),rp],[-1,0,1],shape=(nz,nz),format='csr')
    # test conservativity
    assert np.max(abs(Q@np.ones(nz)))<1e-12
    mass=np.exp(logrho-logrho.max());mass/=mass.sum()
    # matrix with row z and column y: evolve Z using forward chain, then Y gaussian drift.
    yi0=np.argmin(abs(y))
    yindex=np.arange(ny)[None,:]
    zshift=z[:,None]*dt/hy
    # new output at y_j comes from previous value y_j-z*dt
    coords_y=np.broadcast_to(yindex-zshift, (nz,ny))
    coords_z=np.broadcast_to(np.arange(nz)[:,None],(nz,ny))
    sigma=np.sqrt(2*dt)/hy
    # scipy finite volume with mass splitting (not positivity conserving precisely when interpolate)
    def step(M):
        out=expm_multiply(Q.T*(dt/2), M)
        out=gaussian_filter1d(out, sigma=sigma, axis=1, mode='constant',truncate=5)
        out=map_coordinates(out, [coords_z,coords_y], order=1,mode='constant',cval=0.,prefilter=False)
        out=expm_multiply(Q.T*(dt/2),out)
        return out
    def integrate_mass(M, duration):
        n=int(round(duration/dt));assert abs(n*dt-duration)<1e-9
        for _ in range(n): M=step(M)
        return M
    m=np.zeros((nz,ny));m[:,yi0]=mass
    m=integrate_mass(m,0.25)
    classes=[y<=-0.5, (y>-0.5)&(y<=0.5),y>0.5]
    stage1=[m[:,sel].sum(axis=1) for sel in classes]
    s1=float(sum(x.sum() for x in stage1))
    after=[]
    for p1 in stage1:
        mm=np.zeros((nz,ny));mm[:,yi0]=p1
        mm=integrate_mass(mm,0.25)
        after.extend([mm[:,sel].sum(axis=1) for sel in classes])
    s2=float(sum(x.sum() for x in after))
    # probability of positive future state at time 1 given position at 0.5
    plus=np.where(z>0,1.,np.where(z==0,.5,0.))
    pred=expm_multiply(Q*0.5,plus)
    table=np.zeros((3,3,3)) # index record q1,q2, target minus/zero/plus
    for i,v in enumerate(after):
        a,b=divmod(i,3)
        table[a,b,0]=v@(1-pred)
        table[a,b,2]=v@pred
    total=table.sum()
    tab=table/total
    risk=1-sum(np.max(tab[a,b,:]) for a in range(3) for b in range(3))
    predictor=np.argmax(tab,axis=2).tolist()
    return dict(params={"hz":hz,"hy":hy,"dt":dt,"Lz":Lz,"Ly":Ly,"nz":nz,"ny":ny}, table=tab.tolist(), risk=float(risk),mass_loss_stage1=1-s1,mass_loss_stage2=1-s2,unnormalized_total=total, predictor=predictor,baseline_risk=0.5)

if __name__=='__main__':
    vals=[float(x) for x in sys.argv[1:]]
    r=run(*vals) if vals else run()
    print(json.dumps(r,indent=2))
