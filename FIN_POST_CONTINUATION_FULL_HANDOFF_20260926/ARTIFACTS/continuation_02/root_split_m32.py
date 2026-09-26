import math,numpy as np
from scipy.special import logsumexp
alpha=.745
# same V
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]])
G=V@V.T
def dec(h):
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def norm(c):a=np.array(c,float);return float(a@G@a)
L={}
beta16=alpha/2
with open('/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02/level16_a0745.txt') as f:
 for line in f:
  h,z1,z2,z4=line.split();c=dec(h);L[c]=float(z1)-beta16*norm(c)/32
root=(4,)*8
rows=[]
for a,lya in L.items():
 if sum(a)!=16 or any(x>4 for x in a):continue
 b=tuple(4-x for x in a)
 lyb=L[b]
 rows.append((a,lya+lyb,norm(a)/8.0,sum((x-2)**2 for x in a))) # m32 root delta norm/8
lw=np.array([r[1] for r in rows]);ln=logsumexp(lw);p=np.exp(lw-ln)
H=-float(np.sum(p*np.log(p)));ipr=float(np.sum(p*p))
print('ordered_states',len(rows),'entropy',H,'Neff',math.exp(H),'invIPR',1/ipr)
for col,name in [(2,'root_delta'),(3,'imbalance')]:
 x=np.array([r[col] for r in rows],float);mean=float(np.sum(p*x));sd=float(np.sqrt(np.sum(p*(x-mean)**2)));o=np.argsort(x);cs=np.cumsum(p[o]);qs=[float(x[o[np.searchsorted(cs,q)]]) for q in [.05,.25,.5,.75,.95]];print(name,mean,sd,qs)
for ii in np.argsort(p)[-12:][::-1]: print('top',p[ii],rows[ii][0],'delta',rows[ii][2],'imb',rows[ii][3])
