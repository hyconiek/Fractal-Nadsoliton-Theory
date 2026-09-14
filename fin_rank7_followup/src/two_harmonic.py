"""R7P-034/035 invariant k3-cos + k6 alternating stationary family."""
from __future__ import annotations
import math,json
from pathlib import Path
import numpy as np
from scipy.optimize import root
from .model import feature_spaces,dual7
ROOT=Path(__file__).resolve().parents[1]
W,A,L,X,C,A7=feature_spaces(); a=L[3]/6; b=L[6]/12

def xy(J,K):
    den=math.cosh(J)+math.exp(-2*K)
    return math.sinh(J)/den,(math.cosh(J)-math.exp(-2*K))/den

def residual(z,g):
    J,K=z;x,y=xy(J,K)
    return np.array([J-g*a*x,K-g*b*y])

def theta_from_JK(J,K):
    th=np.zeros(7);th[0]=J/math.sqrt(a);th[6]=K/math.sqrt(b);return th

def candidate(J,K,g):
    th=theta_from_JK(J,K);phi,grad,H,p=dual7(th,g,X)
    vals=np.linalg.eigvalsh(H)
    return {'g':g,'J':J,'K':K,'theta7':th.tolist(),'stationarity_residual_2d':float(np.linalg.norm(residual((J,K),g))),
            'stationarity_residual_7d':float(np.linalg.norm(grad)),'phi':phi,'pmax':float(p.max()),
            'H7_eigenvalues':vals.tolist(),'H7_index':int(np.sum(vals<-1e-9)),'H7_nullity':int(np.sum(abs(vals)<=1e-9))}

def roots_at(g,grid=15):
    Js=np.linspace(-g*a,g*a,grid);Ks=np.linspace(-g*b,g*b,grid)
    found=[]
    for z0 in ((J,K) for J in Js for K in Ks):
        sol=root(lambda z:residual(z,g),z0,method='hybr')
        if not sol.success or np.linalg.norm(residual(sol.x,g))>1e-10:continue
        J,K=sol.x
        if abs(J)>g*a+1e-7 or abs(K)>g*b+1e-7:continue
        if not any(np.linalg.norm(sol.x-np.array([q['J'],q['K']]))<1e-7 for q in found):
            found.append(candidate(float(J),float(K),g))
    found.sort(key=lambda q:(q['J'],q['K']))
    return found

def build():
    gains=[3.7,3.71834489812038,4.0,5.0,5.2,6.0]
    allr={str(g):roots_at(g,19) for g in gains}
    out={'exact_reduction':{
        'x':'sinh(J)/(cosh(J)+exp(-2K))','y':'(cosh(J)-exp(-2K))/(cosh(J)+exp(-2K))',
        'stationarity':['J=g*(lambda3/6)*x','K=g*(lambda6/12)*y'],
        'omitted_modes':'zero exactly because the probability is period-4 and reflection-even'},
        'gains':allr}
    (ROOT/'results/R7P-034_035_two_harmonic_stationary_atlas.json').write_text(json.dumps(out,indent=2)+'\n')
    return out
if __name__=='__main__':
    d=build()
    for g,rs in d['gains'].items():
        print('g',g,'roots',len(rs))
        for q in rs: print('  J,K',q['J'],q['K'],'index',q['H7_index'],'eig2',q['H7_eigenvalues'][:3],'res',q['stationarity_residual_7d'])
