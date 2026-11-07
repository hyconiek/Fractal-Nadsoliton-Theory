"""Independent reconstruction of the finite-rank canonical generators."""
from functools import lru_cache
import sympy as s
import numpy as np
z,p,c,r,w,g4,d=s.symbols('z p c r w g4 d',real=True)
coords=(z,p,c,r)
def bracket(a,b):
    return s.diff(a,z)*s.diff(b,p)-s.diff(a,p)*s.diff(b,z)+s.diff(a,c)*s.diff(b,r)-s.diff(a,r)*s.diff(b,c)
H2=(p*p+w*z*z+r*r+3*c*c)/2
H3=d*(6*z*z*c+4*c**3)
H4=9*g4*z**4/70+18*z*z*c*c/35+9*c**4/70
Dlt=3-4*w
A=2*(2*w-3)/Dlt
G3=d*(A*z*z*r+12/Dlt*z*p*c+4/Dlt*p*p*r-s.Rational(4,3)*c*c*r-s.Rational(8,27)*r**3)
basis=[z**a*p**b*c**h*r**k for a in range(5) for b in range(5-a) for h in range(5-a-b) for k in [4-a-b-h] if (a+b)%2==0]
def average(poly):
    ans=0
    for powers,coef in s.Poly(s.expand(poly),*coords).terms():
        if any(v%2 for v in powers):continue
        a,b,h,k=[v//2 for v in powers]
        def factor(i,j):return s.factorial(2*i)*s.factorial(2*j)/(4**(i+j)*s.factorial(i)*s.factorial(j)*s.factorial(i+j))
        ans+=coef*factor(a,b)*factor(h,k)*w**(-a)*s.Rational(3)**(-h)*(p*p+w*z*z)**(a+b)*(r*r+3*c*c)**(h+k)
    return s.expand(ans)
def vector(poly):
    poly=s.Poly(s.expand(poly),*coords)
    return s.Matrix([poly.coeff_monomial(v) for v in basis])
@lru_cache(None)
def construct():
    D=s.Matrix.hstack(*(vector(bracket(b,H2)) for b in basis))
    Pi=s.Matrix.hstack(*(vector(average(b)) for b in basis))
    Q=s.expand(H4+bracket(H3,G3)/2)
    Z=average(Q)
    # The three oscillator-degree blocks are independent.
    rhs=vector(Q-Z);answer=s.zeros(19,1)
    blocks=[]
    for degree in [0,2,4]:
        inds=[i for i,b in enumerate(basis) if sum(s.Poly(b,*coords).monoms()[0][:2])==degree]
        M=(D+Pi).extract(inds,inds);blocks.append(M)
        solution=M.inv()*rhs.extract(inds,[0])
        for i,val in zip(inds,solution):answer[i]=s.factor(val)
    G4=s.expand(sum(v*b for v,b in zip(answer,basis)))
    return D,Pi,Q,Z,G4,blocks
@lru_cache(None)
def vector_fields(wval,g4val):
    G4=construct()[4];dv=3*np.pi*np.sqrt(1.5)/64
    out=[]
    for G in [G3,G4]:
        grad=[s.diff(G,p),-s.diff(G,z),s.diff(G,r),-s.diff(G,c)]
        out.append(s.lambdify(coords,[v.subs({w:wval,g4:g4val,d:dv}) for v in grad],'numpy'))
    return out

def pullback(state,wval,g4val,steps=16):
    """OLD=Phi3(Phi4(NEW)): first undo Phi3, then undo Phi4."""
    y=np.array(state,dtype=float);h=-1/steps
    for fun in vector_fields(wval,g4val):
        def F(q):return np.asarray(fun(*q),dtype=float)
        for _ in range(steps):
            k1=F(y);k2=F(y+h*k1/2);k3=F(y+h*k2/2);k4=F(y+h*k3)
            y+=h*(k1+2*k2+2*k3+k4)/6
    return y
