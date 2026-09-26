import numpy as np
Q=12
th=2*np.pi*np.arange(Q)/Q
h0=np.sqrt(6)/12
rng=np.random.default_rng(20260925)
K=rng.normal(size=(Q,Q)); np.fill_diagonal(K,0.0)
# net first jump moment a=sum K_ij(e_j-e_i)
a=K.sum(axis=0)-K.sum(axis=1)
# hidden orthonormal Fourier coordinates
j=np.arange(Q)
H=np.column_stack([
 np.sqrt(2/Q)*np.cos(2*np.pi*j/Q),
 np.sqrt(2/Q)*np.sin(2*np.pi*j/Q),
 np.sqrt(2/Q)*np.cos(4*np.pi*j/Q),
 np.sqrt(2/Q)*np.sin(4*np.pi*j/Q)])
eta=H.T@a
# Add arbitrary unknown translation-invariant intrinsic baselines.
base1=0.73125; base2=-1.284
Y1=base1+3*h0*(eta[0]*np.cos(th)+eta[1]*np.sin(th))
Y2=base2+3*h0*(eta[2]*np.cos(2*th)+eta[3]*np.sin(2*th))
rec=np.array([
 (2/Q*np.sum(Y1*np.cos(th)))/(3*h0),
 (2/Q*np.sum(Y1*np.sin(th)))/(3*h0),
 (2/Q*np.sum(Y2*np.cos(2*th)))/(3*h0),
 (2/Q*np.sum(Y2*np.sin(2*th)))/(3*h0)])
print('true_hidden_first_moment',repr(eta))
print('recovered_without_baseline_subtraction',repr(rec))
print('max_abs_error',np.max(np.abs(rec-eta)))
assert np.max(np.abs(rec-eta))<2e-14
