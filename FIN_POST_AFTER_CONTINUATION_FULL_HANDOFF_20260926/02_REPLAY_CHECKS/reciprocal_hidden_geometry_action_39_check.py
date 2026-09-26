#!/usr/bin/env python3
import numpy as np

rng=np.random.default_rng(20260926)
n=9

# Connected cycle with positive nonuniform weights.
w=rng.uniform(.4,1.5,size=n)
W=np.zeros((n,n))
for i in range(n):
    j=(i+1)%n
    W[i,j]=W[j,i]=w[i]
A=np.diag(W.sum(axis=1))-W

u=rng.normal(size=n)
psi=rng.normal(size=n)
beta=.37
eps=.08

D=np.diag(u)
dA=beta/2*(D@A+A@D-np.diag(A@u))

# Edge-form interaction.
Eint=0.0
for i in range(n):
    j=(i+1)%n
    Eint += eps*beta/4*W[i,j]*(u[i]+u[j])*(psi[i]-psi[j])**2

# Matrix identity: Eint = eps/2 psi^T dA psi.
Emat=eps*.5*psi@dA@psi
assert abs(Eint-Emat)<1e-12

# Gradient wrt psi equals eps dA psi.
def edge_energy(p,uu):
    val=0.0
    for i in range(n):
        j=(i+1)%n
        val += .5*W[i,j]*(1+eps*beta*(uu[i]+uu[j])/2)*(p[i]-p[j])**2
    return val

h=1e-6
gnum=np.zeros(n)
for k in range(n):
    e=np.zeros(n);e[k]=h
    gnum[k]=(edge_energy(psi+e,u)-edge_energy(psi-e,u))/(2*h)
gan=(A+eps*dA)@psi
assert np.max(np.abs(gnum-gan))<2e-9

# Gradient wrt u of interaction.
gunum=np.zeros(n)
for k in range(n):
    e=np.zeros(n);e[k]=h
    gunum[k]=(edge_energy(psi,u+e)-edge_energy(psi,u-e))/(2*h)

guan=np.zeros(n)
for i in range(n):
    for j in range(n):
        if W[i,j]>0:
            guan[i]+=eps*beta/4*W[i,j]*(psi[i]-psi[j])**2
assert np.max(np.abs(gunum-guan))<2e-9

# Positivity lower bound.
rho=abs(eps*beta)*np.max(np.abs(u))
assert rho<1
E0=.5*psi@A@psi
Eu=edge_energy(psi,u)
assert Eu >= (1-rho)*E0-1e-12

print("PASS")
print("interaction_matrix_identity_error",float(abs(Eint-Emat)))
print("psi_gradient_error",float(np.max(np.abs(gnum-gan))))
print("u_gradient_error",float(np.max(np.abs(gunum-guan))))
print("small_coupling_rho",float(rho))
print("energy_lower_bound_margin",float(Eu-(1-rho)*E0))
