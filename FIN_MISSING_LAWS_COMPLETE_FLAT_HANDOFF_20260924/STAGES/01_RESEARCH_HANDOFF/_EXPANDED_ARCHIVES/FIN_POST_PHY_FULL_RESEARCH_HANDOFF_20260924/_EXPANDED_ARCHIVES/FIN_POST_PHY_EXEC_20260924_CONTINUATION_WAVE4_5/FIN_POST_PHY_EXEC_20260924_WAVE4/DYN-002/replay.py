#!/usr/bin/env python3
import numpy as np,json
pi=np.array([.4,.6]);a=.3;b=a*pi[0]/pi[1]
Q=np.array([[-a,a],[b,-b]])
assert np.allclose(pi@Q,0)
rows=[]
for c in [.2,1,7]:
 Qc=c*Q; assert np.allclose(pi@Qc,0); vals=np.linalg.eigvals(Qc); gap=-min(vals.real)
 rows.append({'c':c,'nonzero_decay_rate':float(gap)})
print(json.dumps(rows,indent=2))
