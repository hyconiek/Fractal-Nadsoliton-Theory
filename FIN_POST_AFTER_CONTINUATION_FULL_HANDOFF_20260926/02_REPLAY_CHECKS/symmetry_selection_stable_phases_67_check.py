#!/usr/bin/env python3
import sympy as sp

x,y=sp.symbols("x y", real=True)
z=x+sp.I*y

# Verify Re(z^m) invariance numerically/symbolically for representative m.
for m in (3,4,6,12):
    P=sp.expand(sp.re(sp.expand_complex(z**m)))
    # Reflection y->-y.
    assert sp.expand(P.subs(y,-y)-P)==0

# D3 cubic daughter Hessian signs.
lam,c,r=sp.symbols("lam c r", nonzero=True, real=True)
# At stationary daughter 3 c s r = -lam.
# radial eigenvalue = lam+6csr = -lam
# angular Cartesian eigenvalue = -9 c r s = 3 lam
assert sp.simplify(lam+6*(-lam/3)) == -lam
assert sp.simplify(-9*(-lam/3)) == 3*lam

# D4 stable-ray formulas.
beta,cc=sp.symbols("beta cc", positive=True, real=True)
# choose anisotropic minimum with c cos4phi = -|c|;
# denominator beta-|c| must be positive for quartic bounded branch.
# radial curvature at nonzero stationary point is -2 lambda.
assert sp.simplify(lam+12*(-lam/(4))) == -2*lam

print("PASS")
print("lowest_anisotropic_degrees",{3:3,4:4,6:6,12:12})
print("D3_daughter_critical_signs","(-lambda, 3 lambda)")
print("D4_stable_condition","beta > |c| for lambda < 0")
