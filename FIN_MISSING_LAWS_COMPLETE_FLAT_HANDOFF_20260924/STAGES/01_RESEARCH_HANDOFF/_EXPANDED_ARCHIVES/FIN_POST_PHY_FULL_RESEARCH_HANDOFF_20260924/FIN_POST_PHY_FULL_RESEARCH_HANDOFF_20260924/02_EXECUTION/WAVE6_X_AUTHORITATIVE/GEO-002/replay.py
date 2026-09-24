import math
for a in [1,2,3]:
 prev=None
 for m in [20,40,80,160,320]:
  nx=a*m; ny=m; N=nx*ny
  lx=N*2*(1-math.cos(2*math.pi/nx)); ly=N*2*(1-math.cos(2*math.pi/ny))
  ex=abs(lx-4*math.pi**2/a); ey=abs(ly-4*math.pi**2*a)
  if prev is not None: assert ex+ey < prev
  prev=ex+ey
 assert abs((4*math.pi**2*a)/(4*math.pi**2/a)-a*a)<1e-12
print('PASS GEO-002 rectangular stationary-mixture counterexample')
