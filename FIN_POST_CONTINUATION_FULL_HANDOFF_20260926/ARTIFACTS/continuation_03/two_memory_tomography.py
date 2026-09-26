import sympy as sp
Q=12
A,B,t=sp.symbols('A B t', positive=True, real=True)
th=[2*sp.pi*i/Q for i in range(Q)]
def c(k,i):return sp.cos(k*th[i])
def s(k,i):return sp.sin(k*th[i])
def S(h,f):return sp.factor(sp.simplify(sum(h[i]*sp.Rational(1,Q)*(f[j]-f[i])**4 for i in range(Q) for j in range(Q))))
f=[A*c(3,i)+t*B*c(4,i) for i in range(Q)]
K={}
for l in (1,2):
    kc=S([c(l,i) for i in range(Q)],f); ks=S([s(l,i) for i in range(Q)],f)
    K[l]=kc
    print('l',l,'cos',kc,'sin',ks)
assert sp.factor(K[1]-3*A*B*t*(10*A*A+9*B*B*t*t))==0
assert sp.factor(K[2]-9*A*A*B*B*t*t)==0
# parity separation
assert sp.simplify((K[1].subs(t,-t)+K[1])/2)==0
assert sp.simplify((K[2].subs(t,-t)-K[2])/2)==0
# exact 12-shift character orthogonality / tomography Gram
for l in (1,2):
    M=sp.Matrix([[sp.cos(l*2*sp.pi*a/Q),sp.sin(l*2*sp.pi*a/Q)] for a in range(Q)])
    G=sp.simplify(M.T*M)
    print('l',l,'translation_gram',G,'rank',M.rank())
    assert G==sp.eye(2)*6 and M.rank()==2
# Every nontrivial hidden Fourier linear term is killed by uniform shift average.
for l in (1,2):
    zc=sp.simplify(sum(sp.cos(l*2*sp.pi*a/Q) for a in range(Q)))
    zs=sp.simplify(sum(sp.sin(l*2*sp.pi*a/Q) for a in range(Q)))
    print('character_sum',l,zc,zs); assert zc==0 and zs==0
