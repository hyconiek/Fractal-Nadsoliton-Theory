import sympy as sp
Q=12; u=sp.Rational(1,Q)
th=[2*sp.pi*j/Q for j in range(Q)]
def c(k,j): return sp.cos(k*th[j])
def s(k,j): return sp.sin(k*th[j])
def p12(x): return sp.simplify(x)

def soft_taylor(h,scale):
    hs=[scale*x for x in h]
    m2=sp.simplify(sum(x*x for x in hs)/Q)
    return [u*x for x in hs],[u*sp.Rational(1,2)*(x*x-m2) for x in hs]

def run(h):
    A1,A2=soft_taylor(h,sp.Integer(1))
    B1,B2=soft_taylor(h,sp.Rational(1,2))
    firstA=[];secondA=[];firstB=[];secondB=[]
    for a in range(Q):
        f=[c(5,(j-a)%Q) for j in range(Q)]
        r1a=r2a=r1b=r2b=sp.Integer(0)
        for i in range(Q):
            for j in range(Q):
                K=(f[j]-f[i])**4
                # A: p=softmax(eps h), q=u
                r1a += A1[i]*u*K
                r2a += A2[i]*u*K
                # B: p=q=softmax(eps h/2)
                r1b += (B1[i]*u+u*B1[j])*K
                r2b += (B2[i]*u+u*B2[j]+B1[i]*B1[j])*K
        firstA.append(sp.expand_trig(r1a)); secondA.append(sp.expand_trig(r2a))
        firstB.append(sp.expand_trig(r1b)); secondB.append(sp.expand_trig(r2b))
    diffs=[sp.simplify(sp.trigsimp(firstA[a]-firstB[a])) for a in range(Q)]
    meanA=sp.trigsimp(sum(secondA)/Q)
    meanB=sp.trigsimp(sum(secondB)/Q)
    return diffs,sp.simplify(meanA),sp.simplify(meanB)
for hk in (1,2):
  for name,fun in [('c',c),('s',s)]:
    h=[fun(hk,j) for j in range(Q)]
    d,a,b=run(h)
    ok=all(x==0 for x in d)
    print(name,hk,'first_equal',ok,'second_shift_mean_A',a,'second_shift_mean_B',b,'difference',sp.simplify(b-a))
