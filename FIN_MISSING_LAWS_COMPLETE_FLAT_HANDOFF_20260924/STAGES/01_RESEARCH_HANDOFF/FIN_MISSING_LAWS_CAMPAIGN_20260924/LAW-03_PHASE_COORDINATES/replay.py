#!/usr/bin/env python3
import numpy as np, json
M=np.array([[-1,1,0],[1,-2,1],[4,-3,0]],dtype=int);t=np.array([3,4,5])
out={'det':int(round(np.linalg.det(M))),'translation_response':(M@t).tolist(),'inverse':np.linalg.inv(M).tolist(),'locked_test':(M@np.array([3,4,5])).tolist(),'identity_energy_statement':'for Delta alpha=ell>0 and c=kappa/ell, edge energy=0.5*kappa*ell'}
print(json.dumps(out,indent=2,sort_keys=True))
