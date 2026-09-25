import numpy as np
# Irregular rooted tree: leaves 0,1 under u; leaf2 and node v under root; leaves3,4 under v.
# masses and kappa strictly increasing down branches.
m=np.array([0.7,1.1,0.9,1.4,0.8])
# nodes: leaves 0..4, u=5, v=6, root=7
children={5:[0,1],6:[3,4],7:[5,2,6]}
kappa={7:0.8,5:1.7,6:2.3}
M={i:m[i] for i in range(5)}
for v in (5,6,7): M[v]=sum(M[c] for c in children[v])
# edge conductances according to theorem
edges=[]
for v,cs in children.items():
    for c in cs:
        if c<5:g=kappa[v]*M[c]
        else:g=M[c]/(1/kappa[v]-1/kappa[c])
        edges.append((v,c,g))
# graph Laplacian
n=8; L=np.zeros((n,n))
for a,b,g in edges:
    L[a,a]+=g;L[b,b]+=g;L[a,b]-=g;L[b,a]-=g
B=list(range(5)); I=[5,6,7]
S=L[np.ix_(B,B)]-L[np.ix_(B,I)]@np.linalg.inv(L[np.ix_(I,I)])@L[np.ix_(I,B)]
# build target hierarchy quadratic by basis probing
def target_energy(f):
    mean={i:f[i] for i in range(5)}
    for v in (5,6,7): mean[v]=sum(M[c]*mean[c] for c in children[v])/M[v]
    E=0.0
    for v,cs in children.items():
        E+=0.5*kappa[v]*sum(M[c]*(mean[c]-mean[v])**2 for c in cs)
    return E
T=np.zeros((5,5))
for i in range(5):
    e=np.zeros(5);e[i]=1
    T[i,i]=2*target_energy(e)
for i in range(5):
  for q in range(i+1,5):
    e=np.zeros(5);e[i]=e[q]=1
    T[i,q]=T[q,i]=target_energy(e)-0.5*T[i,i]-0.5*T[q,q]
print('irregular_max_error',np.max(np.abs(S-T)))
# regular pole/residue ratios
b=2.; r=np.sqrt(2.); A0=r*r/(b*(r*r-1))
for d in range(1,5):
    lam=r**(2*d)*(1+r*r)/(r*r-1)
    mult=(b-1)*b**(d-1)
    Wtot=(b-1)*A0*A0*r**(4*d)
    print(d,lam,mult,Wtot,Wtot/(lam*lam))
