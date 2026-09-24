import numpy as np, math
n=24; h=2*np.pi/n
M=h*np.eye(n); K=np.zeros((n,n))
for i in range(n):
 j=(i+1)%n; c=1/h; K[i,i]+=c; K[j,j]+=c; K[i,j]-=c; K[j,i]-=c
lam=np.linalg.eigvalsh(np.linalg.solve(M,K)); l=float(lam[1])
# same spatial M,K, two internal kinetic tensors; relative frequencies differ
H=np.diag([1.,2.,3.,4.])
for G in [np.eye(4),np.diag([1.,2.,4.,8.])]:
    vals=np.linalg.eigvalsh(np.linalg.solve(G,H))
    print('internal omega ratios',np.sqrt(vals/vals[0]))
# energy-conserving gyroscopic family additionally shows nonuniqueness if T reversal is not imposed
for g in [0,.5,1,2]:
    d=math.sqrt(g*g+4*l); wp=(d+g)/2; wm=(d-g)/2
    print('g',g,'split',wp/wm)
print('NL-04 PASS replay: inequivalent kinetic choices exist.')
