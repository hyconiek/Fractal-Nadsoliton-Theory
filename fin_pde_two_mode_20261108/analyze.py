from pathlib import Path
import json,hashlib
import numpy as np
H=Path(__file__).resolve().parent

def coefficients(w):
    def density(k,channel):
        D=(1+1j*k)*(2+1j*k);S2=np.pi*k/np.sinh(np.pi*k/2);S3=np.pi*(1+k*k)/(2*np.cosh(np.pi*k/2));S4=(4+k*k)*S2/6
        M={'a':3j*k*(k*k+4)*S4/(32*D),'b':1j*k*(k*k+4)*(2-k*k)*S2/(32*D),'x':(3/(2*np.sqrt(2)))*(1+k*k)*(3+k*k)*S3/(24*D)}[channel]
        return abs(M)**2/(2*k)
    return np.array([36*density(np.sqrt(4*w-4),'a')/w,12*density(np.sqrt(8),'b'),72*density(np.sqrt(2*np.sqrt(3*w)-1),'x')/np.sqrt(3*w)])

def estimates(r,start):
    h=[v for v in r['samples'] if v['t']>=start];t=np.array([v['t'] for v in h]);I=np.array([v['I'] for v in h]);J=np.array([v['J'] for v in h]);ka,kb,kx=coefficients(r['w'])
    if r['preparation']=='shape':return {'kb':float(np.polyfit(t,1/J,1)[0]),'max_I':float(max(I))}
    def integ(v):return np.r_[0,np.cumsum(np.diff(t)*(v[1:]+v[:-1])/2)]
    x=integ(I*J)
    yi=I[0]-I-ka*integ(I*I);yj=J[0]-J-kb*integ(J*J)
    return {'kx_from_I':float(np.polyfit(x,yi,1)[0]),'kx_from_J':float(np.polyfit(x,yj,1)[0]),'mixed_work_integral':float(x[-1])}
if __name__=='__main__':
    runs=json.loads((H/'results.json').read_text())['runs'];assert len(runs)==6
    rows=[]
    for r in runs:
        rows.append({'preparation':r['preparation'],'epsilon':r['epsilon'],'dx':r['dx'],'dt':r['dt'],'estimates_25_100':estimates(r,25),'estimates_50_100':estimates(r,50),'energy_drift':r['relative_global_energy_drift'],'flux_residual':max(r['max_relative_local_balance_residual'])})
    out={'theory':dict(zip(['ka','kb','kx'],coefficients(1.5))),'runs':rows,'scope':'Fixed self coefficients in mixed regression; not a simultaneous unbiased identification of all coefficients.'}
    (H/'comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
