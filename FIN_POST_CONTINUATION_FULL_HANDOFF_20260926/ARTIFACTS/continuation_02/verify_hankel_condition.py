import numpy as np

def moments(a,l):
    return np.array([sum(ai*li**n for ai,li in zip(a,l)) for n in range(4)],float)

def recover(m):
    A=np.array([[m[1],-m[0]],[m[2],-m[1]]],float)
    b=np.array([m[2],m[3]],float)
    s1,s2=np.linalg.solve(A,b)
    return np.sort(np.roots([1,-s1,s2]).real)

print('near-pole campaign')
for sep in [1.0,0.2,0.05,0.01]:
    a=[1.,1.]; l=[1.,1.+sep]; m=moments(a,l)
    D=m[0]*m[2]-m[1]**2
    eps=1e-8*np.max(abs(m)); mp=m+eps*np.array([1,-1,1,-1])
    err=np.max(abs(recover(mp)-l))
    print(sep,D,err)
    assert abs(D-a[0]*a[1]*sep**2)<1e-12

print('weak-residue campaign')
for a2 in [1.,.2,.05,.01,.001]:
    a=[1.,a2]; l=[1.,2.]; m=moments(a,l)
    D=m[0]*m[2]-m[1]**2
    eps=1e-8*np.max(abs(m)); mp=m+eps*np.array([1,-1,1,-1])
    err=np.max(abs(recover(mp)-l))
    print(a2,D,err)
    assert abs(D-a2)<1e-12
