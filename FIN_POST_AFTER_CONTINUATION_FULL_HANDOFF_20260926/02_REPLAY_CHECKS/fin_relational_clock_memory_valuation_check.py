#!/usr/bin/env python3
import numpy as np

def fit_gamma(S,D,n=30):
    x=np.log(S[:n]); y=np.log(D[:n])
    return float(np.polyfit(x,y,1)[0])

# Generic base p=t^m(a+t^q b), wave-like m=4,q=2.
a=np.array([1.,2.,3.])
b=np.array([.2,-.1,.05])
m,q=4,2

# t=1e-4...1e-1 remains asymptotic for this polynomial fixture while avoiding
# catastrophic cancellation in normalized-profile differences.
ts=np.logspace(-4,-1,300)

def profile(P, leading):
    S=P.sum(axis=1)
    R=P/S[:,None]
    r0=leading/leading.sum()
    D=np.linalg.norm(R-r0,axis=1)
    return S,D

P=ts[:,None]**m*(a+ts[:,None]**q*b)
S,D=profile(P,a)
g0=fit_gamma(S,D)
assert abs(g0-.5)<2e-4

# Memoryless invertible map: order q/m survives.
M=np.array([[1.,.1,0.],[.05,1.,.1],[0.,.05,1.]])
Y=P@M.T
Sy,Dy=profile(Y,M@a)
gm=fit_gamma(Sy,Dy)
assert abs(gm-.5)<2e-4

# Direct feedthrough + generic constant memory kernel K0.
K0=np.array([[0.,1.,0.],[0.,0.,.5],[.2,0.,0.]])
Y2=P.copy()
Y2 += ts[:,None]**(m+1)/(m+1)*(K0@a)
Y2 += ts[:,None]**(m+q+1)/(m+q+1)*(K0@b)
S2,D2=profile(Y2,a)
gd=fit_gamma(S2,D2)
assert abs(gd-.25)<2e-3

# Shape-preserving K0=cI: the +1 term is parallel to the leading ray.
Kp=.7*np.eye(3)
Y3=P.copy()
Y3 += ts[:,None]**(m+1)/(m+1)*(Kp@a)
Y3 += ts[:,None]**(m+q+1)/(m+q+1)*(Kp@b)
S3,D3=profile(Y3,a)
gp=fit_gamma(S3,D3)
assert abs(gp-.5)<2e-3

# Pure-memory filter K(u)=K0+K1 u. Generic K1 changes shape at +1.
K1=np.array([[.1,0.,0.],[0.,.3,0.],[0.,0.,.6]])
Y4=ts[:,None]**(m+1)/(m+1)*(K0@a)
Y4 += ts[:,None]**(m+2)/((m+1)*(m+2))*(K1@a)
S4,D4=profile(Y4,K0@a)
gpm=fit_gamma(S4,D4)
assert abs(gpm-.2)<2e-3

print("PASS")
print("base_gamma",g0)
print("memoryless_gamma",gm)
print("direct_generic_memory_gamma",gd)
print("direct_shape_preserving_memory_gamma",gp)
print("pure_memory_gamma",gpm)
