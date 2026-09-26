import sympy as sp
Q=12; jj=range(Q)
l3,l4,l5,l6=sp.symbols('l3 l4 l5 l6', positive=True)
ls=[l3,l3,l4,l4,l5,l5,l6]
E=[]
for k in (3,4,5):
 E += [sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.cos(2*sp.pi*k*j/Q) for j in jj]),
       sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.sin(2*sp.pi*k*j/Q) for j in jj])]
E += [sp.Matrix([(-1)**j/sp.sqrt(Q) for j in jj])]
Y=[]
for k in (1,2):
 Y += [sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.cos(2*sp.pi*k*j/Q) for j in jj]),
       sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.sin(2*sp.pi*k*j/Q) for j in jj])]
# X_r=sqrt(lambda_r) E_r
T=[]
for ya in Y:
 M=sp.zeros(7)
 for r in range(7):
  for s in range(7):
   M[r,s]=sp.simplify(sp.sqrt(ls[r]*ls[s])*sum(ya[j]*E[r][j]*E[s][j] for j in jj))
 T.append(M)
Gram=sp.Matrix(4,4,lambda a,b:sp.simplify(sp.trace(T[a]*T[b])))
print('Gram=');print(Gram)
print('diag factored',[sp.factor(Gram[i,i]) for i in range(4)])
Mstrength=sp.factor(sum(Gram[i,i] for i in range(4))/Q)
print('M_uniform=',Mstrength)
# nonzero entries per channel
for a,M in enumerate(T):
 print('T',a)
 for r in range(7):
  for s in range(r,7):
   if sp.simplify(M[r,s])!=0: print(r,s,sp.simplify(M[r,s]))
