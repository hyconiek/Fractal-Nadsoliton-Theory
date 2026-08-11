#!/usr/bin/env python3
"""Dependency-minimal independent replay of the ST108 Krawczyk certificate.

Third-party dependencies are limited to NumPy and mpmath.  No FIN research
module, SciPy solver, or prior interval helper is imported.  The numerical
center is read from the published ST108 certificate; every interval kernel,
equation, Jacobian, and Krawczyk enclosure is recomputed here.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "FIN_ST108_Transcendental_Fold_Certificate.json"
OUTPUT = ROOT / "FIN_ST120_Minimal_Transcendental_Replay_Certificate.json"
N = 12


def iv(value):
    if isinstance(value, tuple):
        return mp.iv.mpf([str(value[0]), str(value[1])])
    if isinstance(value, (float, np.floating)):
        return mp.iv.mpf([str(float(np.nextafter(value,-np.inf))),str(float(np.nextafter(value,np.inf)))])
    return mp.iv.mpf(str(value))


def bounds(x):
    return float(np.nextafter(float(x.a),-np.inf)),float(np.nextafter(float(x.b),np.inf))


def prod(a0,a1,b0,b1):
    vals=[a0*b0,a0*b1,a1*b0,a1*b1]
    return float(np.nextafter(min(vals),-np.inf)),float(np.nextafter(max(vals),np.inf))


def matvec(alo,ahi,xlo,xhi):
    lo=np.zeros(alo.shape[0]); hi=np.zeros(alo.shape[0])
    for i in range(alo.shape[0]):
        l=h=0.0
        for j in range(alo.shape[1]):
            q0,q1=prod(alo[i,j],ahi[i,j],xlo[j],xhi[j])
            l=np.nextafter(l+q0,-np.inf); h=np.nextafter(h+q1,np.inf)
        lo[i],hi[i]=l,h
    return lo,hi


def left_product(r,alo,ahi):
    lo=np.zeros((r.shape[0],alo.shape[1])); hi=np.zeros_like(lo)
    for i in range(r.shape[0]):
        for j in range(alo.shape[1]):
            l=h=0.0
            for k in range(r.shape[1]):
                q0,q1=prod(r[i,k],r[i,k],alo[k,j],ahi[k,j])
                l=np.nextafter(l+q0,-np.inf); h=np.nextafter(h+q1,np.inf)
            lo[i,j],hi[i,j]=l,h
    return lo,hi


def strict_interval_matrix():
    omega=iv("0.18575"); phi=iv("0.16250"); eta=iv(9)/iv(5)
    w={d:mp.iv.cos(omega*d+phi)/(1+iv(d)**eta) for d in range(1,7)}
    s=2*sum((w[d] for d in range(1,6)),iv(0))+w[6]
    a=[[iv(0) for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            a[i][j]=s if i==j else -w[min((i-j)%N,(j-i)%N)]
    return a,w,s


def interval_system(center,radius,a):
    q=[iv((center[i]-radius,center[i]+radius)) for i in range(7)]
    kappa=iv((center[7]-radius,center[7]+radius))
    v=[iv((center[8+i]-radius,center[8+i]+radius)) for i in range(7)]
    u=[]; w=[]
    for site in range(N):
        index=site if site<=6 else N-site
        u.append(q[index]); w.append(v[index])
    au=[sum((a[i][j]*u[j] for j in range(N)),iv(0)) for i in range(N)]
    aw=[sum((a[i][j]*w[j] for j in range(N)),iv(0)) for i in range(N)]
    h=[]; rdiag=[]; drdu=[]
    for x in u:
        rho=x**2; den=1+rho/2; qfun=rho/den
        qp=den**-2; qpp=-(den**-3); qppp=iv("1.5")*den**-4
        hh=-qfun*qp+iv("0.075"); hp=-(qp**2+qfun*qpp); hpp=-(3*qp*qpp+qfun*qppp)
        h.append(hh); rdiag.append(2*hh+4*rho*hp); drdu.append(2*x*(6*hp+4*rho*hpp))
    g=[kappa*au[i]+2*u[i]*h[i] for i in range(7)]
    g += [kappa*aw[i]+rdiag[i]*w[i] for i in range(7)]
    g += [iv("0.5")*(sum((x*x for x in v),iv(0))-1)]
    jac=[[iv(0) for _ in range(15)] for _ in range(15)]
    for i in range(7):
        for col in range(7):
            total=iv(0)
            for site in range(N):
                if (site if site<=6 else N-site)==col: total += kappa*a[i][site]
            if i==col: total += rdiag[i]
            jac[i][col]=total; jac[7+i][8+col]=total
        jac[i][7]=au[i]; jac[7+i][i]=drdu[i]*w[i]; jac[7+i][7]=aw[i]; jac[14][8+i]=v[i]
    glo=np.array([bounds(x)[0] for x in g]); ghi=np.array([bounds(x)[1] for x in g])
    jlo=np.array([[bounds(x)[0] for x in row] for row in jac]); jhi=np.array([[bounds(x)[1] for x in row] for row in jac])
    return glo,ghi,jlo,jhi


def run_replay(write: bool=True):
    mp.iv.dps=70
    source=json.loads(INPUT.read_text(encoding="utf-8")); center=np.array(source["center"],dtype=float)
    a,w,s=strict_interval_matrix()
    g0lo,g0hi,j0lo,j0hi=interval_system(center,0.0,a)
    c=np.linalg.inv(0.5*(j0lo+j0hi))
    attempts=[]; accepted=None
    for radius in [1e-8,3e-9,1e-9,3e-10]:
        _,_,jlo,jhi=interval_system(center,radius,a)
        cglo,cghi=matvec(c,c,g0lo,g0hi)
        ylo=np.nextafter(center-cghi,-np.inf); yhi=np.nextafter(center-cglo,np.inf)
        cjlo,cjhi=left_product(c,jlo,jhi); mlo,mhi=-cjhi,-cjlo
        for i in range(15):
            mlo[i,i]=np.nextafter(mlo[i,i]+1,-np.inf); mhi[i,i]=np.nextafter(mhi[i,i]+1,np.inf)
        mdlo,mdhi=matvec(mlo,mhi,np.full(15,-radius),np.full(15,radius))
        klo=np.nextafter(ylo+mdlo,-np.inf); khi=np.nextafter(yhi+mdhi,np.inf)
        margin=float(min(np.min(klo-(center-radius)),np.min((center+radius)-khi)))
        row={"radius":radius,"minimum_strict_inclusion_margin":margin,"included":margin>0,
             "maximum_Krawczyk_half_width":float(np.max((khi-klo)/2))}
        attempts.append(row)
        if accepted is None and margin>0: accepted=row
    result={
        "dependencies":["Python standard library","NumPy","mpmath"],
        "forbidden_dependencies_not_imported":["SciPy","FIN project research modules"],
        "input_center_file":INPUT.name,
        "declared_parameters":{"omega":"0.18575","phi":"0.16250","eta":"9/5"},
        "row_sum_interval":list(bounds(s)),
        "kernel_widths":{str(d):bounds(w[d])[1]-bounds(w[d])[0] for d in w},
        "attempts":attempts,"accepted":accepted,
        "independent_replay_matches_ST108_acceptance":bool(accepted and source["accepted"]["included"]),
        "scope":"Independent recomputation of intervals, equations, Jacobian and Krawczyk image; the published numerical center remains an input.",
    }
    if write: OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True),encoding="utf-8")
    return result


if __name__=="__main__":
    print(json.dumps(run_replay(True),indent=2,sort_keys=True))
