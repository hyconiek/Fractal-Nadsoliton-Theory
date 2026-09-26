#!/usr/bin/env python3
import math, numpy as np
G=.5; lam=1.
for s in (2.,3.,5.,10.):
    disc=s*s-4*G*s
    g1=(s+math.sqrt(max(0.,disc)))/2
    g2=(s-math.sqrt(max(0.,disc)))/2
    c=s/lam
    R=(lam/s)*np.outer([g1,g2],[g1,g2])
    static=np.diag([g1,g2])-R/lam
    print('s',s,'g1',g1,'g2',g2,'c',c,'Rtrace',np.trace(R),'R12',R[0,1],
          'static',static.tolist())
    assert np.max(np.abs(static-G*np.array([[1,-1],[-1,1]])))<1e-12
    assert abs(R[0,1]-lam*G)<1e-12
    assert abs(np.trace(R)-lam*(s-2*G))<1e-12
