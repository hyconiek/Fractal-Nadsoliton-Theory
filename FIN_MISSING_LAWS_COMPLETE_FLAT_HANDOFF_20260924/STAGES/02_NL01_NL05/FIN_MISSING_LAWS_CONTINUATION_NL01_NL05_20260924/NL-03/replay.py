import itertools, math, numpy as np
R=np.array([[2,-1,0,0,0],[1,1,-1,0,0],[1,0,1,-1,0],[1,0,0,1,-1]],int)
k=np.arange(1,6)
assert np.all(R@k==0) and np.linalg.matrix_rank(R)==4
mins=[]
for cols in itertools.combinations(range(5),4): mins.append(round(np.linalg.det(R[:,cols])))
g=0
for x in mins: g=math.gcd(g,abs(x))
assert g==1
beta=np.array([0,0,1,-2,1]); gamma=np.array([0,0,4,-3,0])
assert np.all(np.array([0,0,1,-1])@R==beta)
assert np.all(np.array([-1,-1,3,0])@R==gamma)
ev=np.linalg.eigvalsh(R.T@R)
assert abs(ev[0])<1e-10 and ev[1]>.999999
print('NL-03 PASS replay',ev)
