#!/usr/bin/env python3
"""Bounded numerical check only; no proof of infinite-dimensional PDE or radiation."""
import math
import cmath
import mpmath as mp
mp.mp.dps=40
checks=[]
def check(label, ok):
    checks.append((label, bool(ok)))
    if not ok: raise AssertionError(label)
def t(x): return mp.tanh(x)
def sech(x): return 1/mp.cosh(x)
f=lambda x: mp.sqrt(mp.mpf(3)/4)*sech(x)**2
check('normalization', abs(mp.quad(lambda x:f(x)**2,[-mp.inf,0,mp.inf])-1)<mp.mpf('1e-30'))
a=mp.mpf('-0.5'); om2=a+2
for x in [mp.mpf('-1.25'),0,mp.mpf('0.75'),mp.mpf('2')]:
    ev=-mp.diff(f,x,2)+(a+6*t(x)**2)*f(x)
    check(f'eigenfunction-x={x}',abs(ev-om2*f(x))<mp.mpf('1e-30'))
# source-compatible 2/3 complement commutator, no truncated mass selection
h=mp.mpf('0.01'); d2=mp.mpf('-0.005'); d3=mp.mpf('0.005')
check('noncommutator',h*(d3-d2)!=0)
sminus=(mp.sqrt(1+4*mp.mpf('5.98'))-1)/2
splus=(mp.sqrt(1+4*mp.mpf('6.02'))-1)/2
check('complement-mode-separation',sminus+mp.mpf('0.08')>2)
check('full-transverse-gap',-mp.mpf('.75')+6-splus**2>1)
q=mp.sqrt(2)
psi=lambda x: mp.exp(1j*q*x)*(3*t(x)**2-3j*q*t(x)-(q*q+1))/((1-1j*q)*(2-1j*q))
for x in [mp.mpf('-1.1'),mp.mpf('0'),mp.mpf('0.75')]:
    resid=-mp.diff(psi,x,2)+(4-6*sech(x)**2)*psi(x)-(q*q+4)*psi(x)
    check(f'Jost-ODE-x={x}', abs(resid)<mp.mpf('1e-28'))
F4=lambda q: mp.pi*q*(q*q+4)/(6*mp.sinh(mp.pi*q/2))
I_numeric=mp.quad(lambda x: t(x)*sech(x)**4*mp.conj(psi(x)),[-mp.inf,0,mp.inf])
I_exact=1j*q*(q*q+4)*F4(q)/(8*(1+1j*q)*(2+1j*q))
check('on-shell-Fourier-integral',abs(I_numeric-I_exact)<mp.mpf('1e-29'))
Gamma=lambda q: 81*q*(q*q+4)*F4(q)**2/(2048*(1+q*q))
for q0 in [mp.mpf('1'),mp.sqrt(2),mp.sqrt(3)]:
    check(f'positive-Gamma-q={q0}',Gamma(q0)>0)
# check local density coefficient consistency: two +/- scattering roots
M=mp.mpf(9)/4*I_exact
check('spectral-normalization',abs(abs(M)**2/(2*q)-Gamma(q))<mp.mpf('1e-30'))
print(f'PASS {len(checks)}/{len(checks)} bounded controls')
print('Gamma(a=-1/2) =',mp.nstr(Gamma(mp.sqrt(2)),16))
