import numpy as np
for q in [12,18,24,30]:
    for L in range(0,min(5,(q-2)//2)+1):
        th=2*np.pi*np.arange(q)/q
        rows=[np.ones(q)]
        for k in range(1,L+1): rows += [np.cos(k*th),np.sin(k*th)]
        C=np.vstack(rows)
        rank=np.linalg.matrix_rank(C,tol=1e-10)
        assert rank==2*L+1
        assert q-rank==q-(2*L+1)
print('NL-01 PASS replay: hierarchy verified; no L=2 selector encoded.')
