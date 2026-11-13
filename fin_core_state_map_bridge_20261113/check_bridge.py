"""Exact symbolic checks for the conditional A11 -> RP state-map bridge."""
from pathlib import Path
import json
import sympy as s
H=Path(__file__).resolve().parent
n=12
q,phi,w,kappa=s.symbols('q phi w kappa',real=True)
m,mp,lam,c4,c6,y=s.symbols('m mp lam c4 c6 y',real=True)
k=s.symbols('k1:12',real=True)
K=s.Matrix(n,n,lambda i,j:0 if i==j else k[(j-i)%n-1])
M=(K+K.T)/2
u=s.ones(n,1)/s.sqrt(n)
v=s.Matrix([s.cos(s.pi*i/6)/s.sqrt(6) for i in range(n)])
P=lambda z:s.eye(n)-z*z.T
simp=lambda a:a.applyfunc(s.simplify) if isinstance(a,s.MatrixBase) else s.simplify(a)
checks={}
def check(name,value):
    checks[name]=bool(value)
    assert checks[name],name
check('uniform_unit_norm',simp((u.T*u)[0])==1)
check('cosine_unit_norm',simp((v.T*v)[0])==1)
check('circulant_uniform_eigenvector',simp(M*u-sum(k)*u)==s.zeros(n,1))
check('circulant_cosine_eigenvector',simp(P(v)*M*v)==s.zeros(n,1))
def power(z,p):return z.applyfunc(lambda a:a**p)
for p in [3,5]:check(f'uniform_power_{p}_closure',simp(P(u)*power(u,p))==s.zeros(n,1))
# Pull back all local terms and the full (possibly directed) mixing potential.
fields=s.symbols('z0:12',real=True);z=s.Matrix(fields)
V=mp*phi**2/2+lam*phi**4/4+sum(m*a*a/2+c4*a**4/4+c6*a**6/6+y*phi**2*a*a for a in fields)+(z.T*K*z)[0]/2
sub={fields[i]:u[i]*q for i in range(n)}
Vred=simp(V.subs(sub,simultaneous=True))
expected=mp*phi**2/2+lam*phi**4/4+(m+sum(k))*q*q/2+c4*q**4/48+c6*q**6/864+y*phi**2*q*q
check('restricted_potential',simp(Vred-expected)==0)
F=s.Matrix([s.diff(V,a).subs(sub,simultaneous=True) for a in fields])
check('all_twelve_full_gradient_components',simp(F-u*s.diff(expected,q))==s.zeros(n,1))
check('phi_full_gradient',simp(s.diff(V,phi).subs(sub,simultaneous=True)-s.diff(expected,phi))==0)
Vr=simp(expected.subs({mp:-2,lam:2,y:3,c4:12,c6:144,m:w-2-sum(k)})+s.Rational(1,2))
Vtarget=(phi**2-1)**2/2+(w-2+6*phi**2)*q*q/2+q**4/4+q**6/6
check('RP_exact_potential_match_with_added_coefficients',simp(Vr-Vtarget)==0)
check('RP_psi_force',simp(s.diff(Vr,q)-((w-2+6*phi**2)*q+q**3+q**5))==0)
check('RP_phi_force',simp(s.diff(Vr,phi)-(2*phi*(phi**2-1)+6*phi*q*q))==0)
check('fourth_response',s.diff(Vr,q,4).subs(q,0)==6)
check('sixth_response',s.diff(Vr,q,6)==120)
check('mixed_response',s.diff(Vr,phi,2,q,2)==12)
r3=simp(P(v)*power(v,3));norm3=simp((r3.T*r3)[0])
check('cosine_cubic_leakage_squared',norm3==s.Rational(1,576))
check('cosine_quartic_projection',simp(sum(a**4 for a in v))==s.Rational(1,8))
check('cosine_sextic_projection',simp(sum(a**6 for a in v))==s.Rational(5,288))
delta=s.Symbol('delta',real=True);shift=s.zeros(n,1);shift[0]=delta*u[0]**3
rh=simp(P(u)*shift)
check('inhomogeneous_quartic_breaks_uniform_closure',simp((rh.T*rh)[0])==11*delta**2/20736)
# Change an independently supplied local family, not the kernel.
DV=q**4/4
check('different_quartic_same_linearization',all(s.diff(DV,q,j).subs(q,0)==0 for j in range(3)))
check('different_quartic_different_force',s.diff(DV,q)==q**3)
check('kernel_has_no_direct_fourth_response',s.diff((z.T*K*z)[0]/2,fields[0],4)==0)
out={'status':'CONDITIONAL_EXACT_EMBEDDING_NOT_STRICT_COEFFICIENT_SOURCE','checks':checks,'check_count':len(checks),'uniform_effective_coefficients':{'mass_squared':'m + sum(k1,...,k11)','g4':'c4/12','g6':'c6/144','Phi2_psi2_coefficient':'y'},'cosine_cubic_transverse_norm_squared':str(norm3),'cosine_projected_g4_factor':'1/8','cosine_projected_g6_factor':'5/288','nonuniform_quartic_transverse_norm_squared':'11*delta**2/20736'}
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
