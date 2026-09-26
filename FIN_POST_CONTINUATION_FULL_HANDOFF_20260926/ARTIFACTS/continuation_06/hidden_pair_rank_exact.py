import sympy as sp
Q=12
j=range(Q)
E=[]
for k in (3,4,5):
    E += [sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.cos(2*sp.pi*k*n/Q) for n in j]),
          sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.sin(2*sp.pi*k*n/Q) for n in j])]
E += [sp.Matrix([(-1)**n/sp.sqrt(Q) for n in j])]
H=[]
for k in (1,2):
    H += [sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.cos(2*sp.pi*k*n/Q) for n in j]),
          sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.sin(2*sp.pi*k*n/Q) for n in j])]
cols=[];labels=[]
for a in range(7):
    prod=sp.matrix_multiply_elementwise(E[a],E[a])
    cols.append(sp.Matrix([sp.simplify(h.dot(prod)) for h in H]));labels.append((a,a))
for a in range(7):
  for b in range(a+1,7):
    prod=sp.sqrt(2)*sp.matrix_multiply_elementwise(E[a],E[b])
    cols.append(sp.Matrix([sp.simplify(h.dot(prod)) for h in H]));labels.append((a,b))
B=sp.Matrix.hstack(*cols)
BB=sp.simplify(B*B.T)
print('rank',B.rank())
print('BBstar')
sp.print_latex(BB)
print(BB)
print('eigenvals',BB.eigenvals())
# four explicit witness products
for pair in [(0,2),(0,3),(2,6),(3,6)]:
 a,b=pair; prod=sp.matrix_multiply_elementwise(E[a],E[b]); v=[sp.simplify(h.dot(prod)) for h in H]
 print('pair',pair,'hidden_coords',v)
