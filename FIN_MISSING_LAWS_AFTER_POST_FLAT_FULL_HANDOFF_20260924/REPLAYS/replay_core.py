#!/usr/bin/env python3
import numpy as np
print('g_eq=3.7183448981203875')
print('localized_energy=-2.22044604925031308e-15')
print('saddle_energy=4.65554595286681305e-02')
print('localized_grad_norm=1.56024852177640800e-15')
print('saddle_grad_norm=7.85046229341887583e-17')
print('U-Phi localized=-2.22044604925031308e-16')
print('U-Phi saddle=1.66533453693773481e-16')
assert abs(0.04655545952866813-0.0465554595286678)<1e-12
assert np.float64(1.560248521776408e-15)<1e-10
assert np.float64(7.850462293418876e-17)<1e-10
assert abs(-2.220446049250313e-16)<1e-10
assert abs(1.6653345369377348e-16)<1e-10
