#!/usr/bin/env python3
import numpy as np, math
from scipy.optimize import root_scalar

l4=2.19956884933321
l3=1.96140686197644
l5=2.2986062720790903
l6=2.3421820411463

def R(t):
    e=np.exp(t)
    return (e-1)/(e+2)

def Rp(t):
    e=np.exp(t)
    return 3*e/(e+2)**2

sol=root_scalar(lambda t:R(t)-t*Rp(t),bracket=(0.1,5),
                xtol=1e-14)
t=sol.root
c=t/R(t)
gfold=4*c/l4
A4=math.sqrt(l4/6)
sfold=2*t/(3*A4)

assert abs(gfold-4.993057757778402)<2e-13
assert abs(sfold-1.140334204309495)<1e-11

thresholds={
    6:12/l6,
    5:12/l5,
    4:12/l4,
    3:12/l3,
}
assert thresholds[6] < thresholds[5] < thresholds[4] < thresholds[3]

# Exact moment/cumulant signs on discrete uniform carrier.
j=np.arange(12)
def cumulant4(k,scale):
    x=scale*np.cos(2*np.pi*k*j/12) if k!=6 else scale*((-1.)**j)
    x=x-x.mean()
    return np.mean(x**4)-3*np.mean(x**2)**2

A5=math.sqrt(l5/6)
A6=math.sqrt(l6/12)
kap5=cumulant4(5,A5)
kap6=cumulant4(6,A6)
assert kap5<0 and kap6<0

print("PASS")
print("pure_k4_fold_t",t)
print("pure_k4_fold_g",gfold)
print("pure_k4_fold_s",sfold)
print("uniform_thresholds",thresholds)
print("k5_D4Phi",-kap5)
print("k6_D4Phi",-kap6)
