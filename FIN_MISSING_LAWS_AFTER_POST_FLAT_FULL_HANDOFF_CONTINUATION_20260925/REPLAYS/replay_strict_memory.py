import math, numpy as np
from scipy.special import softmax
N=12
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/(1+min(abs(i-j),N-abs(i-j))**1.8) for j in range(N)] for i in range(N)])
A=np.diag(W.sum(1))-W
L=np.fft.fft(A[0]).real[:7]
j=np.arange(N)
cols=[]
for k in (3,4,5):
    cols += [np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/N),np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/N)]
cols += [np.sqrt(L[6]/12)*(-1.0)**j]
X=np.column_stack(cols)
Y=[]
for k in (1,2):
    Y += [np.sqrt(2/N)*np.cos(2*np.pi*k*j/N),np.sqrt(2/N)*np.sin(2*np.pi*k*j/N)]
Y=np.column_stack(Y)
def ps(s):
    th=np.zeros(7); th[[0,2,4,6]]=s
    return softmax(X@th)
states={
'uniform':np.ones(N)/N,
'saddle':ps([0.9409570673214491,1.0014394023775437,0.9621088536282347,0.6864149839950513]),
'localized':ps([1.8199035812800828,1.913989554668724,1.914569132546848,1.367203280195504])}
Ts=np.array([X.T@np.diag(Y[:,a])@X for a in range(4)])
GT=np.einsum('aij,bij->ab',Ts,Ts)
for name,p in states.items():
    S=np.diag(p)-np.outer(p,p)
    F=X.T@S@X; H=Y.T@S@X; G=Y.T@S@Y
    Sig=G-H@np.linalg.solve(F,H.T)
    M=float(np.sum(Sig*GT))
    print(name,'Sigma_eigs',np.linalg.eigvalsh(Sig),'M',M)
