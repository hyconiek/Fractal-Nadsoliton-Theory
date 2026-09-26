import sympy as sp
Q=12
L=sp.symbols('L', positive=True, real=True)

def c(k,j): return sp.cos(2*sp.pi*k*j/Q)
def s(k,j): return sp.sin(2*sp.pi*k*j/Q)

def coupling(k,Lk,hidden='c',probe='c',exceptional=False):
    h=[c(2,i) if hidden=='c' else s(2,i) for i in range(Q)]
    if exceptional:
        f=[sp.sqrt(Lk/12)*(-1)**i for i in range(Q)]
    else:
        base=c if probe=='c' else s
        f=[sp.sqrt(Lk/6)*base(k,i) for i in range(Q)]
    return sp.simplify(sum(h[i]*sp.Rational(1,Q)*(f[j]-f[i])**4
                           for i in range(Q) for j in range(Q)))

for k in (3,4,5):
    print('cos-hidden cos-probe',k,coupling(k,L))
print('cos-hidden k6',coupling(6,L,exceptional=True))
assert coupling(3,L)==0 and coupling(4,L)==0 and coupling(6,L,exceptional=True)==0
assert sp.simplify(coupling(5,L)-L**2/3)==0

# Exact phase law is checked by expanding rotated basis symbolically in sampled
# exact Z12 sums for independent symbols C=cos(varphi), S=sin(varphi).
C,S=sp.symbols('C S', real=True)
f=[sp.sqrt(L/6)*(C*c(5,i)-S*s(5,i)) for i in range(Q)]
hc=[c(2,i) for i in range(Q)]
hs=[s(2,i) for i in range(Q)]
Sc=sp.factor(sum(hc[i]*sp.Rational(1,Q)*(f[j]-f[i])**4 for i in range(Q) for j in range(Q)))
Ss=sp.factor(sum(hs[i]*sp.Rational(1,Q)*(f[j]-f[i])**4 for i in range(Q) for j in range(Q)))
print('hidden_cos polynomial',Sc)
print('hidden_sin polynomial',Ss)
# On C=cos(varphi),S=sin(varphi): Sc=L^2/3*(C^2-S^2), Ss=2L^2/3*C*S.
assert sp.simplify(Sc-L**2*(C**2-S**2)*(C**2+S**2)/3)==0
assert sp.simplify(Ss-2*L**2*C*S*(C**2+S**2)/3)==0
print('unit-circle law: Sc=L^2/3 cos(2varphi), Ss=L^2/3 sin(2varphi)')
