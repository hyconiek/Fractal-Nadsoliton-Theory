"""Exact virial fluxes, constants and critical-window accounting controls."""
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,condition):
    checks.append({'name':name,'pass':bool(condition)})
    if not condition: raise AssertionError(name)
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
    check('input: '+row['path'],hashlib.sha256((HERE.parent/row['path']).read_bytes()).hexdigest()==row['sha256'])
x=s.symbols('x',real=True);m=s.symbols('m',positive=True)
w=s.Function('w')(x);p=s.Function('p')(x);u=s.Function('u')(x)
phi=s.atan(x);zeta=1/s.sqrt(1+x*x);a=(1+x*x)**-3
v=u/zeta
raw=s.diff(phi,x)*s.diff(v,x)**2-s.diff(phi,x,3)*v*v/4
reduced=s.diff(u,x)**2+(x*x-1)*u*u/(2*(1+x*x)**2)
check('coercivity density differs by exact flux',s.simplify(raw-reduced-s.diff(x*u*u/(1+x*x),x))==0)
r,y=s.symbols('r y',real=True)
coef=s.integrate(s.integrate(y-r,(y,x,1)),(r,-1,x))/2
check('variance Fubini coefficient',s.simplify(coef-(1-x*x)/2)==0)
check('position coercivity residual',s.Rational(3,4)-s.Rational(122,7)/24==s.Rational(1,42))
check('trace coercivity residual',s.Rational(1,4)-s.Rational(33,14)/24==s.Rational(17,112))
B=s.diff(phi,x)*s.diff(w,x)**2-s.diff(phi,x,3)*w*w/4
Idot=(phi*s.diff(p,x)+s.diff(phi,x)*p/2)*p+(phi*s.diff(w,x)+s.diff(phi,x)*w/2)*(s.diff(w,x,2)-m*m*w)
flux=phi*p*p/2+phi*s.diff(w,x)**2/2+s.diff(phi,x)*w*s.diff(w,x)/2-s.diff(phi,x,2)*w*w/4-m*m*phi*w*w/2
check('I virial identity with spatial flux',s.simplify(Idot+B-s.diff(flux,x))==0)
Jdot=a*p*p+a*w*(s.diff(w,x,2)-m*m*w)
Jbulk=a*p*p-a*s.diff(w,x)**2+s.diff(a,x,2)*w*w/2-m*m*a*w*w
check('J identity with spatial flux',s.simplify(Jdot-Jbulk-s.diff(a*w*s.diff(w,x)-s.diff(a,x)*w*w/2,x))==0)
ratio=s.simplify(s.diff(a,x,2)/a)
check('weight curvature ratio',s.simplify(ratio-(-6+42*x*x)/(1+x*x)**2)==0)
z=s.symbols('z',nonnegative=True);R=(-6+42*z)/(1+z)**2
check('curvature maximum critical point',s.simplify(s.diff(R,z)-6*(9-7*z)/(z+1)**3)==0)
check('curvature extrema within bound',R.subs(z,0)==-6 and R.subs(z,s.Rational(9,7))==s.Rational(147,16) and s.Rational(147,16)<12)
check('J coefficient includes mass',s.Rational(23,4)+6+s.Rational(9,2)==s.Rational(65,4)<17)
check('delta leaves positive position budget',s.Rational(1,24)-s.Rational(17,816)==s.Rational(1,48))
check('trace coefficient after absorption',s.Rational(1,2)/(s.Rational(1,816)/2)==816)
t,T,beta,eps,mu,S=s.symbols('t T beta eps mu S',positive=True)
check('weighted constant endpoint',s.simplify(s.exp(-beta*T)+beta*s.integrate(s.exp(-beta*t),(t,0,T)))==1)
linear=s.exp(-beta*T)*T+beta*s.integrate(t*s.exp(-beta*t),(t,0,T))
check('weighted linear budget',s.simplify(linear-(1-s.exp(-beta*T))/beta)==0)
F=t+s.sin(t)/2;G=t+s.Rational(1,2)
check('differentiating cumulative bound fails',s.diff(F,t).subs(t,0)>s.diff(G,t).subs(t,0))
check('critical forcing budget',s.simplify(eps**4*S/eps**2)==S*eps**2)
check('unweighting factor',s.simplify(2*mu*eps**2*S/eps**2)==2*mu*S)
check('quartic trace exponent',s.simplify(eps**8*(S/eps**2)**s.Rational(5,2))==S**s.Rational(5,2)*eps**3)
result={'status':'PASS','scope':'Exact symbolic identities, constants, exponent accounting and input integrity; not independent proof certification.','check_count':len(checks),'checks':checks,'sympy':s.__version__}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks)}))
