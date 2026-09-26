import sympy as sp
Q=12; u=sp.Rational(1,Q); t=sp.Rational(1,3)
th=[2*sp.pi*j/Q for j in range(Q)]
def c(k,j):return sp.cos(k*th[j])
def s(k,j):return sp.sin(k*th[j])
P=tuple(map(sp.Rational,[2,3,5,7])); T=tuple(map(sp.Rational,[-1,4,6,-2]))
hp=[P[0]*c(1,j)+P[1]*s(1,j)+P[2]*c(2,j)+P[3]*s(2,j) for j in range(Q)]
hq=[T[0]*c(1,j)+T[1]*s(1,j)+T[2]*c(2,j)+T[3]*s(2,j) for j in range(Q)]
def response(power,shift,tt):
    f=[c(3,(j-shift)%Q)+tt*c(4,(j-shift)%Q) for j in range(Q)]
    return sp.simplify(sum(u*u*(hp[i]+hq[j])*(f[j]-f[i])**power for i in range(Q) for j in range(Q)))
def parity(power,odd):
    rp=[response(power,r,t) for r in range(Q)];rm=[response(power,r,-t) for r in range(Q)]
    return [sp.simplify((a-b)/2 if odd else (a+b)/2) for a,b in zip(rp,rm)]
def proj(v,k,kind):
    b=[c(k,r) if kind=='c' else s(k,r) for r in range(Q)]
    return sp.simplify(sum(x*y for x,y in zip(v,b))/6)
Q1=parity(4,True); Q2=parity(4,False); C1=parity(3,False); C2=parity(3,True)
KQ1=sp.factor(t*(9*t*t+10)/4);KQ2=sp.factor(3*t*t/4);KC1=sp.factor(-3*t*t/8);KC2=sp.factor(-3*t/4)
S1=(sp.simplify(proj(Q1,1,'c')/KQ1),sp.simplify(proj(Q1,1,'s')/KQ1))
S2=(sp.simplify(proj(Q2,2,'c')/KQ2),sp.simplify(proj(Q2,2,'s')/KQ2))
D1=(sp.simplify(proj(C1,1,'c')/KC1),sp.simplify(proj(C1,1,'s')/KC1))
D2=(sp.simplify(proj(C2,2,'c')/KC2),sp.simplify(proj(C2,2,'s')/KC2))
Prec=tuple(sp.simplify((x+y)/2) for pair in zip(S1+S2,D1+D2) for x,y in [pair])
Trec=tuple(sp.simplify((x-y)/2) for pair in zip(S1+S2,D1+D2) for x,y in [pair])
print('sum_k1',S1,'sum_k2',S2)
print('diff_k1',D1,'diff_k2',D2)
print('recovered_departure',Prec)
print('recovered_target',Trec)
assert Prec==P and Trec==T
# Shift-averaged O(eps^2) cross term for a pure c5 quartic.
def k5(sh):return [c(5,(j-sh)%Q) for j in range(Q)]
cross=sp.simplify(sum(sum((u*hp[i])*(u*hq[j])*(k5(r)[j]-k5(r)[i])**4 for i in range(Q) for j in range(Q)) for r in range(Q))/Q)
pred=sp.Rational(3,16)*(P[2]*T[2]+P[3]*T[3])
print('second_order_shift_mean_cross',cross,'predicted',pred)
assert sp.simplify(cross-pred)==0
