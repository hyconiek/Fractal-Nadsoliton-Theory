#!/usr/bin/env python3
import math, numpy as np
from scipy.special import softmax

N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
L=np.fft.fft(A[0]).real[:7]
j=np.arange(N)

cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/N),
        np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/N)
    ]
cols += [np.sqrt(L[6]/12)*(-1.)**j]
X=np.column_stack(cols)

hidden=[]
for k in (1,2):
    hidden += [
        np.sqrt(2/N)*np.cos(2*np.pi*k*j/N),
        np.sqrt(2/N)*np.sin(2*np.pi*k*j/N)
    ]
Y=np.column_stack(hidden)

s=np.array([
  1.8199035812800826819822480178121342,
  1.9139895546687240691680856665430389,
  1.9145691325468480852004194017926659,
  1.3672032801955039010722866317201310
])
theta=np.array([s[0],0,s[1],0,s[2],0,s[3]])
p=softmax(X@theta)
S=np.diag(p)-np.outer(p,p)

F=X.T@S@X
H=Y.T@S@X
G=Y.T@S@Y
Sigma=G-H@np.linalg.solve(F,H.T)

# Uniform exact structure numerically to roundoff.
pu=np.ones(N)/N
Su=np.diag(pu)-np.outer(pu,pu)
Fu=X.T@Su@X
Hu=Y.T@Su@X
Gu=Y.T@Su@Y
assert np.linalg.norm(Hu)<1e-14
assert np.linalg.norm(Gu-np.eye(4)/12)<1e-14

# Positive Schur complement and covariance block diagonalization.
assert np.linalg.eigvalsh(F).min()>0
assert np.linalg.eigvalsh(Sigma).min()>0
K=H@np.linalg.inv(F)

# Joint covariance transformed by z=y-Kx.
C=np.block([[F,H.T],[H,G]])
T=np.block([
    [np.eye(7),np.zeros((7,4))],
    [-K,np.eye(4)]
])
Ct=T@C@T.T
assert np.linalg.norm(Ct[:7,7:])<1e-12
assert np.linalg.norm(Ct[7:,7:]-Sigma)<1e-12

# Schur loss is PSD.
loss=G-Sigma
assert np.linalg.eigvalsh(loss).min()>-1e-12

# Whitening.
evals,evecs=np.linalg.eigh(Sigma)
Sinvhalf=evecs@np.diag(1/np.sqrt(evals))@evecs.T
white=Sinvhalf@Sigma@Sinvhalf.T
assert np.linalg.norm(white-np.eye(4))<1e-11

# Phase tangent Fisher matrix used by report 44.
t=np.array([0,3*s[0],0,4*s[1],0,5*s[2],0])
vtheta=X@t
Z=np.column_stack([vtheta,Y])
GF=Z.T@S@Z
hidden_eigs=np.linalg.eigvalsh(G)
phase_cross=np.linalg.norm(GF[0,1:])

assert abs(GF[0,0]-2.644071152609)<2e-12
assert abs(np.linalg.cond(G)-8.951269366)<2e-8

print("PASS")
print("localized_probability_max",float(p.max()))
print("phase_fisher",float(GF[0,0]))
print("raw_hidden_fisher_eigs",hidden_eigs.tolist())
print("phase_hidden_cross_norm",float(phase_cross))
print("raw_hidden_condition",float(np.linalg.cond(G)))
print("visible_hidden_H_opnorm",float(np.linalg.norm(H,2)))
print("conditional_hidden_eigs",np.linalg.eigvalsh(Sigma).tolist())
print("conditional_hidden_condition",float(np.linalg.cond(Sigma)))
print("uniform_hidden_error",float(np.linalg.norm(Gu-np.eye(4)/12)))
print("whitening_error",float(np.linalg.norm(white-np.eye(4))))
