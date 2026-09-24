#!/usr/bin/env python3
import numpy as np
q=192
j=np.arange(q)
fine=np.exp(2j*np.pi*q*(2*j)/(2*q))
coarse=np.exp(2j*np.pi*(q/2)*j/q)
print("max mismatch",np.max(np.abs(fine-coarse)))
assert np.max(np.abs(fine-coarse))>1.9
