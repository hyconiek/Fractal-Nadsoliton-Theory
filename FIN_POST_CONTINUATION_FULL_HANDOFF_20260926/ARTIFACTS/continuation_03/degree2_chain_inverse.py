import sympy as sp
z=sp.symbols('z', positive=True)

def series(g,Y): return sp.factor(g*Y/(g+Y))
def build(gs,cs):
    # gs length m+1; cs length m
    Y=z*cs[-1]+gs[-1]
    for k in range(len(cs)-2,-1,-1):
        Y=z*cs[k]+series(gs[k+1],Y)
    return sp.factor(series(gs[0],Y))

def reconstruct(Y,g0,m):
    gs=[sp.Rational(g0)]; cs=[]; cur=sp.factor(g0*Y/(g0-Y))
    for k in range(m):
        c=sp.simplify(sp.limit(cur/z,z,sp.oo)); cs.append(c)
        rem=sp.factor(cur-z*c)
        g=sp.simplify(sp.limit(rem,z,sp.oo)); gs.append(g)
        if k<m-1: cur=sp.factor(g*rem/(g-rem))
    return gs,cs

gs=list(map(sp.Rational,[2,3,5,7]));cs=list(map(sp.Rational,[11,13,17]))
Y=build(gs,cs)
rgs,rcs=reconstruct(Y,gs[0],len(cs))
print('Y=',Y)
print('gs=',rgs,'cs=',rcs)
assert rgs==gs and rcs==cs
# zero-storage degree-2 vertex collapses: g1,g2 in series.
ga,gb=sp.Rational(3),sp.Rational(5)
geq=sp.factor(ga*gb/(ga+gb))
Yzero=series(sp.Rational(2), series(ga,series(gb,z*sp.Rational(7)+sp.Rational(11))))
Ycollapsed=series(sp.Rational(2), series(geq,z*sp.Rational(7)+sp.Rational(11)))
print('zero_storage_series_geq',geq,'difference',sp.factor(Yzero-Ycollapsed))
assert sp.factor(Yzero-Ycollapsed)==0
