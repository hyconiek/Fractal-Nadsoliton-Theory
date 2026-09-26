import sympy as sp
z=sp.symbols('z', positive=True)

def series(g,Y): return sp.factor(g*Y/(g+Y))
def build(gs,cs):
    Y=z*cs[-1]+gs[-1]
    for k in range(len(cs)-2,-1,-1):
        Y=z*cs[k]+series(gs[k+1],Y)
    return sp.factor(series(gs[0],Y))

def reconstruct(Y,m):
    gs=[]; cs=[]; cur=sp.factor(Y)
    # g0 is visible from the high-frequency input admittance
    g0=sp.simplify(sp.limit(cur,z,sp.oo)); gs.append(g0)
    cur=sp.factor(g0*cur/(g0-cur))
    for k in range(m):
        c=sp.simplify(sp.limit(cur/z,z,sp.oo)); cs.append(c)
        rem=sp.factor(cur-z*c)
        g=sp.simplify(sp.limit(rem,z,sp.oo)); gs.append(g)
        if k<m-1: cur=sp.factor(g*rem/(g-rem))
    return gs,cs

gs=list(map(sp.Rational,[2,3,5,7])); cs=list(map(sp.Rational,[11,13,17]))
Y=build(gs,cs); rgs,rcs=reconstruct(Y,len(cs))
print('Y',Y); print('recovered_g',rgs); print('recovered_c',rcs)
assert rgs==gs and rcs==cs
# Zero-storage node: adjacent edges collapse to series conductance.
ga,gb=sp.Rational(3),sp.Rational(5); geq=sp.factor(ga*gb/(ga+gb))
Yzero=series(sp.Rational(2),series(ga,series(gb,z*sp.Rational(7)+sp.Rational(11))))
Ycol=series(sp.Rational(2),series(geq,z*sp.Rational(7)+sp.Rational(11)))
print('zero_storage_geq',geq,'difference',sp.factor(Yzero-Ycol))
assert sp.factor(Yzero-Ycol)==0
