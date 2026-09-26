import sympy as sp
z,lam,d,a,eps,R,lp=sp.symbols('z lam d a eps R lp', positive=True)
x=z+lam
fnear=sp.factor(a/(x-d)+a/(x+d)-2*a/x)
print('near_collision_difference',fnear)
assert sp.factor(fnear-2*a*d**2/(x*(x**2-d**2)))==0
fweak=sp.factor(eps*R/(z+lam)-eps*R/(z+lp))
print('weak_residue_difference',fweak)
assert sp.factor(fweak-eps*R*(lp-lam)/((z+lam)*(z+lp)))==0
# four-leaf additive tree: split 12|34 has four-point margin 2w.
w,x1,x2,x3,x4=sp.symbols('w x1 x2 x3 x4', positive=True)
d12=x1+x2; d34=x3+x4
d13=x1+w+x3;d24=x2+w+x4;d14=x1+w+x4;d23=x2+w+x3
same=sp.expand(d12+d34);cross1=sp.expand(d13+d24);cross2=sp.expand(d14+d23)
print('four_point',same,cross1,cross2,'margin',sp.simplify(cross1-same))
assert sp.simplify(cross1-same)==2*w and sp.simplify(cross2-same)==2*w
