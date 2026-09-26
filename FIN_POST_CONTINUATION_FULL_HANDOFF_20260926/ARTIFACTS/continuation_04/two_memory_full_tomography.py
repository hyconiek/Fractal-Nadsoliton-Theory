import sympy as sp
Q=12
theta=[2*sp.pi*j/Q for j in range(Q)]
t=sp.Rational(1,3)
# arbitrary exact hidden coefficients
a1,b1,a2,b2=map(sp.Rational,[2,3,5,7])
def c(k,j): return sp.cos(k*theta[j])
def s(k,j): return sp.sin(k*theta[j])
h=[a1*c(1,j)+b1*s(1,j)+a2*c(2,j)+b2*s(2,j) for j in range(Q)]

def response(shift,tt):
    # translated mixed retained probe; A=B=1 is enough for rank/tomography proof
    f=[c(3,(j-shift)%Q)+tt*c(4,(j-shift)%Q) for j in range(Q)]
    return sp.simplify(sum(h[i]*sp.Rational(1,Q)*(f[j]-f[i])**4 for i in range(Q) for j in range(Q)))
Rp=[response(r,t) for r in range(Q)]
Rm=[response(r,-t) for r in range(Q)]
Odd=[sp.simplify((x-y)/2) for x,y in zip(Rp,Rm)]
Even=[sp.simplify((x+y)/2) for x,y in zip(Rp,Rm)]
# DFT projections. Gram sum cos^2=sin^2=6.
def proj(v,k,kind):
    basis=[c(k,r) if kind=='c' else s(k,r) for r in range(Q)]
    return sp.simplify(sum(vr*br for vr,br in zip(v,basis))/6)
print('odd_k1_cos',proj(Odd,1,'c'))
print('odd_k1_sin',proj(Odd,1,'s'))
print('odd_k2_norm',sp.simplify(proj(Odd,2,'c')**2+proj(Odd,2,'s')**2))
print('even_k2_cos',proj(Even,2,'c'))
print('even_k2_sin',proj(Even,2,'s'))
print('even_k1_norm',sp.simplify(proj(Even,1,'c')**2+proj(Even,1,'s')**2))
print('uniform_odd',sp.simplify(sum(Odd)/Q),'uniform_even_hidden',sp.simplify(sum(Even)/Q))
# expected K1=3 t(10+9t^2), K2=9t^2 for A=B=1, signs depend on shift convention.
K1=sp.factor(3*t*(10+9*t*t));K2=sp.factor(9*t*t)
print('K1',K1,'K2',K2)
# magnitudes must reconstruct hidden pairs exactly up to the known shift-sign convention.
assert sp.simplify(proj(Odd,1,'c')**2+proj(Odd,1,'s')**2-K1**2*(a1*a1+b1*b1))==0
assert sp.simplify(proj(Even,2,'c')**2+proj(Even,2,'s')**2-K2**2*(a2*a2+b2*b2))==0
assert sp.simplify(sum(Odd))==0 and sp.simplify(sum(Even))==0
