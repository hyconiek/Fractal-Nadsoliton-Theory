import math, json
from pathlib import Path
w=[0.4699856726450201, 0.1920435516901028, 0.09142861427792495, 0.0470291687456504, 0.02413122336363006, 0.011070817321442113]
k2=sum((i+1)**2*x for i,x in enumerate(w)); k4=sum((i+1)**4*x for i,x in enumerate(w))
assert abs(k2-3.8153141154998296)<1e-14
prev=None
for M in [24,48,96,192,384]:
 h=1/M; p=2*math.pi
 lam=2/h**2*sum(x*(1-math.cos(p*(i+1)*h)) for i,x in enumerate(w))
 lim=k2*p*p; err=abs(lam-lim); bound=h*h*p**4*k4/12
 assert err <= bound*(1+1e-12)
 if prev is not None: assert err < prev
 prev=err
theta=math.pi/2
ratio=(2*sum(x*(1-math.cos((i+1)*theta)) for i,x in enumerate(w)))/(k2*theta**2)
assert abs(ratio-0.2107039949061925)<1e-14 and abs(ratio-1)>0.5
print('PASS REF-006',k2,ratio)
