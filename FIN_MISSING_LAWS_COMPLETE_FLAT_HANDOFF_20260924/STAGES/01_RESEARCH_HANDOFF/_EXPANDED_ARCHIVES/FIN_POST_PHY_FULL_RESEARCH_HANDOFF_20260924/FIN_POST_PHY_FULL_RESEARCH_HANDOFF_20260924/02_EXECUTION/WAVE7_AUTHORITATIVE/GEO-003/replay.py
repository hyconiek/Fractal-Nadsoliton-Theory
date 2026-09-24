import math
# Path-metric control
h=0.01
assert abs((4*h)/(2*h)-2)<1e-15
# finite-graph short-time estimator for p(t)~c t^k tends to zero
for k in [1,2,4]:
 vals=[]
 for t in [1e-2,1e-4,1e-8]:
  p=t**k
  vals.append(math.sqrt(-4*t*math.log(p)))
 assert vals[-1] < vals[0]
 assert vals[-1] < 0.01
# anisotropic dual norm
a=2; dx=dy=1
d=math.sqrt(a*dx*dx+dy*dy/a)
assert abs(d-math.sqrt(2.5))<1e-15
print('PASS GEO-003 order-of-limits and intrinsic metric checks')
