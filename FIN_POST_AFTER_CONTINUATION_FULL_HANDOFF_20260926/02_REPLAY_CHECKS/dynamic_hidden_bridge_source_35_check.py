#!/usr/bin/env python3
import math
import numpy as np

N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
j=np.arange(N)

hidden=[]
for k in (1,2):
    hidden += [
        np.sqrt(2/N)*np.cos(2*np.pi*k*j/N),
        np.sqrt(2/N)*np.sin(2*np.pi*k*j/N)
    ]
Y=np.column_stack(hidden)

# Minimal grammar structural theorem on a nonconstant hidden mode.
u=Y[:,0]
D=np.diag(u)

alpha=.371
J=alpha*(D@A+A@D-np.diag(A@u))
assert np.linalg.norm(J-J.T)<1e-12
assert np.linalg.norm(J@np.ones(N))<1e-12

J_bad_sym=.3*D@A+.4*A@D-.3*np.diag(A@u)
assert np.linalg.norm(J_bad_sym-J_bad_sym.T)>1e-6

J_bad_cons=.3*(D@A+A@D)-.2*np.diag(A@u)
assert np.linalg.norm(J_bad_cons@np.ones(N))>1e-6

# Six shell-resolved bridges considered as LINEAR MAPS H4 -> matrices.
shell_map_columns=[]
shell_maps=[]

for d in range(1,7):
    Wd=np.zeros_like(W)
    for i in range(N):
        for k in range(N):
            dist=min(abs(i-k),N-abs(i-k))
            if i!=k and dist==d:
                Wd[i,k]=W[i,k]

    responses=[]
    basis_outputs=[]
    for a in range(4):
        ua=Y[:,a]
        dW=.5*Wd*(ua[:,None]+ua[None,:])
        Ld=np.diag(dW.sum(axis=1))-dW
        assert np.linalg.norm(Ld-Ld.T)<1e-12
        assert np.linalg.norm(Ld@np.ones(N))<1e-12
        responses.append(Ld.reshape(-1))
        basis_outputs.append(Ld)
    shell_map_columns.append(np.concatenate(responses))
    shell_maps.append(basis_outputs)

Map=np.column_stack(shell_map_columns)
rank=np.linalg.matrix_rank(Map,tol=1e-11)
sv=np.linalg.svd(Map,compute_uv=False)
assert rank==6

# Fractional bridge equals the sum of shell bridges, basis direction by basis direction.
max_sum_err=0.0
for a in range(4):
    ua=Y[:,a]
    Da=np.diag(ua)
    Jhalf=.5*(Da@A+A@Da-np.diag(A@ua))
    shell_sum=sum(shell_maps[d][a] for d in range(6))
    max_sum_err=max(max_sum_err,np.linalg.norm(Jhalf-shell_sum))
assert max_sum_err<1e-12

print("PASS")
print("minimal_bridge_symmetry_error",float(np.linalg.norm(J-J.T)))
print("minimal_bridge_zero_mode_error",float(np.linalg.norm(J@np.ones(N))))
print("shell_local_map_family_rank",rank)
print("shell_map_singular_values",sv.tolist())
print("fractional_bridge_shell_sum_error",max_sum_err)
