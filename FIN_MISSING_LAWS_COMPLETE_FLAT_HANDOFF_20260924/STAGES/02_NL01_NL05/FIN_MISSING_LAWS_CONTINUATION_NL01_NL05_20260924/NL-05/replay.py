import numpy as np
N=512; L=2*np.pi; dx=L/N; x=np.arange(N)*dx
k=2*np.pi*np.fft.fftfreq(N,d=dx); lap=4*np.sin(.5*k*dx)**2/dx**2
omega=np.sqrt(1+lap)
d=((x-np.pi+np.pi)%(2*np.pi))-np.pi
q0=np.exp(-.5*(d/.18)**2)*np.cos(12*d); qh=np.fft.fft(q0)
ipr=[]
for t in [0,.5,1,2,4,8]:
 q=np.fft.ifft(qh*np.cos(omega*t)).real; w=q*q; w/=w.sum(); ipr.append(float((w*w).sum()))
 print(t,ipr[-1],float(np.max(np.abs(q))))
assert ipr[-1] < .5*ipr[0]
print('NL-05 PASS negative-baseline replay')
