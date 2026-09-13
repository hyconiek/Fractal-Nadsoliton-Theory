"""Independent intake audit of the September 13 rank-seven handoff.

No supplied CSV is executed or treated as a proof. Exact certificates use
Fraction intervals; numerical reconstruction is explicitly a separate layer.
"""
from fractions import Fraction as F
from pathlib import Path
import csv
import hashlib
import itertools
import math
import re
import sys

import numpy as np
import sympy as sp
from scipy.optimize import root, minimize_scalar
from scipy.special import logsumexp, softmax

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fin_projected_learning.research import strict, certify_strict_spectrum

BUNDLE = ROOT / 'FIN_research_artifacts_pre_and_post_Discord'


class Interval:
    """Small exact rational interval implementation, no floating decisions."""
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError('Reversed interval')

    @staticmethod
    def cast(x):
        return x if isinstance(x, Interval) else Interval(x)

    def __add__(self, other):
        o = self.cast(other)
        return Interval(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        o = self.cast(other)
        v = [x*y for x in (self.lo, self.hi) for y in (o.lo, o.hi)]
        return Interval(min(v), max(v))
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self.cast(other)
        if o.lo <= 0 <= o.hi:
            raise ValueError('Division by interval containing zero')
        return self * Interval(1/o.hi, 1/o.lo)

    def __rtruediv__(self, other):
        return self.cast(other) / self

    def __pow__(self, n):
        if n < 0:
            return 1 / self**(-n)
        result = Interval(1)
        for _ in range(n):
            result = result * self
        return result

    def strings(self):
        return [str(self.lo), str(self.hi)]


def evaluate(expr, values):
    expr = sp.sympify(expr)
    if expr.is_Rational:
        return Interval(F(int(expr.p), int(expr.q)))
    if expr.is_Symbol:
        return values[expr]
    if expr.is_Add:
        return sum((evaluate(x, values) for x in expr.args), Interval(0))
    if expr.is_Mul:
        out = Interval(1)
        for x in expr.args:
            out *= evaluate(x, values)
        return out
    if expr.is_Pow and expr.exp.is_Integer:
        return evaluate(expr.base, values)**int(expr.exp)
    raise ValueError(f'Unsupported exact expression: {expr}')


def bernstein(poly, variable, left, right):
    """Symbolic power-to-Bernstein conversion including variable endpoints."""
    z = sp.Symbol('z')
    p = sp.Poly(sp.expand(poly.subs(variable, left+(right-left)*z)), z)
    n = p.degree()
    return [sp.factor(sum(p.nth(k)*sp.binomial(i,k)/sp.binomial(n,k)
                         for k in range(i+1))) for i in range(n+1)]


def spectrum_intervals():
    w = [Interval(*x) for x in certify_strict_spectrum()['eigenvalue_intervals']]
    return [Interval(0)] + [w[0]-x for x in w[1:]]


def exact_certificates():
    L = spectrum_intervals()
    sigma = (2*L[3]*(L[4]+L[5])-L[4]*L[5])/(24*L[3])
    aI = L[3]/6
    cI = L[6]/(3*sigma)
    t2I = 1-sigma/aI
    assert t2I.lo > 0 and t2I.hi < 1
    assert (1-2*t2I).lo > 0  # First local parity-mixing shift is negative.
    assert (sigma-L[4]/12).lo > 0 and (sigma-L[5]/12).lo > 0
    assert (1-cI/4).lo > 0  # Schur denominator for every q in [0,1].
    l3,l4,l5=sp.symbols('l3 l4 l5',positive=True)
    sigma_symbol=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    t2_symbol=1-6*sigma_symbol/l3
    # Exact double-root: the 45 block has sigma as an eigenvalue;
    # its other eigenvalue is smaller by the interval inequalities above.
    assert sp.factor((sigma_symbol-l4/12)*(sigma_symbol-l5/12)
                     -l4*l5*t2_symbol/144)==0
    a,s,c,r,q = sp.symbols('a s c r q')
    values = {a:aI, s:sigma, c:cI}

    # r = sech(sqrt(lambda3/6)*s3), NOT exp(-sqrt(lambda3/6)*s3).
    d0 = s*(1+r)-a*r**2
    N = sp.expand((1+r)**2*d0-a*r*(1-r)*(1+r)**2-c*r*d0)
    assert sp.degree(N,r) == 3
    face1 = []
    for lo,hi in [(0,sp.Rational(1,2)),(sp.Rational(1,2),1)]:
        row = [evaluate(x,values) for x in bernstein(N,r,lo,hi)]
        assert all(x.lo > 0 for x in row)
        face1.append([x.strings() for x in row])
    assert (2*sigma-aI).lo > 0

    eta = 1-c*q*(1-q)
    # (sigma - scalar Schur channel)*eta, extremized over physical x^2.
    def R(x2):
        return sp.expand(s*eta-a*q*eta-a*x2*(c*(1-q)-1))
    qmin, qcrit = 1-s/(2*a), 1-1/c
    assert evaluate(qmin-sp.Rational(1,2),values).lo > 0
    assert evaluate(1-qcrit,values).lo > 0
    assert evaluate(qcrit-qmin,values).lo > 0
    low = [evaluate(x,values) for x in bernstein(R(2*q-1),q,qmin,qcrit)]
    assert all(x.lo > 0 for x in low)
    # Cancel the exact endpoint zero BEFORE interval evaluation.
    upper_poly = R(1-s/a)
    assert sp.factor(upper_poly.subs(q,1)) == 0
    quotient = sp.cancel(upper_poly/(1-q))
    assert sp.denom(quotient) == 1
    upper = [evaluate(x,values) for x in bernstein(quotient,q,qcrit,1)]
    assert all(x.lo > 0 for x in upper)

    # Rigorous log12 from log3+2log2, atanh series and geometric tails.
    def log_int(n):
        z = F(n-1,n+1); terms=80
        lo=2*sum(z**(2*k+1)/(2*k+1) for k in range(terms))
        tail=2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
        return Interval(lo,lo+tail)
    threshold = 6*(log_int(3)+2*log_int(2))
    rank6 = 2*(L[3]+L[4]+L[5]); rank7=rank6+L[6]
    assert (threshold-rank6).lo > 0 and (rank7-threshold).lo > 0
    assert (1/Interval(F('3.718345'))-sigma).lo > 0
    # Full-seven-coordinate counterexample h_j=2 cos(pi*j/2).
    # x=tanh(1)>3/4 since sum_{k=0}^5 2^k/k! > 7.
    assert sum(F(2)**k/math.factorial(k) for k in range(6)) > 7
    assert L[4].lo > F('2.19') and L[5].lo > F('2.29')
    assert F('2.23')**2 < F('2.19')*F('2.29')
    rayleigh_lower=(F('2.19')+F('2.29')+2*F('2.23')*F(3,4))/24
    assert rayleigh_lower > F(10,37)
    return {
        'laplacian_intervals':[x.strings() for x in L],
        'sigma_interval':sigma.strings(), 't_star_squared_interval':t2I.strings(),
        'schur_eta_lower':str((1-cI/4).lo),
        'resolvent_numerator_bernstein':face1,
        'extreme_face_lower_piece_bernstein':[x.strings() for x in low],
        'extreme_face_upper_quotient_bernstein':[x.strings() for x in upper],
        'rank6_budget_deficit_lower':str((threshold-rank6).lo),
        'rank7_budget_surplus_lower':str((rank7-threshold).lo),
        'full7_two_direction_covariance_lower':str(rayleigh_lower),
        'full7_negative_hessian_margin_at_g_ge_3_7':str(rayleigh_lower-F(10,37)),
        'scope':'Exact strict-spectrum face bounds, NOT off-face/global rank-seven closure.'}


def audit_manifest():
    records=[]
    for line in (BUNDLE/'MANIFEST.txt').read_text().splitlines():
        match=re.fullmatch(r'(.+)\t(\d+) bytes\tsha256=([0-9a-f]{64})',line)
        if not match:
            continue
        name,size,digest=match.groups(); data=(BUNDLE/'artifacts'/name).read_bytes()
        assert len(data)==int(size) and hashlib.sha256(data).hexdigest()==digest
        if name.endswith('.csv'):
            rows=list(csv.reader(data.decode().splitlines()))
            assert rows and all(len(x)==len(rows[0]) for x in rows)
        records.append({'file':name,'sha256':digest,'bytes':int(size)})
    assert len(records)==41
    assert {x['file'] for x in records}=={p.name for p in (BUNDLE/'artifacts').iterdir()}
    return records


def model():
    W=strict(); A=np.diag(W.sum(axis=1))-W
    L=np.fft.fft(A[0]).real[:7]
    j=np.arange(12)
    C=np.column_stack([np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/12) for k in [3,4,5]]
                      +[np.sqrt(L[6]/12)*(-1.)**j])
    X=np.column_stack([f for k in [3,4,5] for f in
                       [np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/12),
                        np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/12)]]
                      +[C[:,3]])
    return W,A,L,C,X


def moments(features, p):
    mu=p@features
    centered=features-mu
    return mu, centered.T@(p[:,None]*centered)


def dual(s,g,C):
    p=softmax(C@s); mu,cov=moments(C,p)
    return s@s/(2*g)-logsumexp(C@s)+math.log(len(C)), s/g-mu, np.eye(len(s))/g-cov,p


def parity(C,s):
    p=softmax(C@s); q=p[::2].sum()
    mp,cp=moments(C[::2],p[::2]/q)
    mm,cm=moments(C[1::2],p[1::2]/(1-q))
    W=q*cp+(1-q)*cm
    b=np.sqrt(q*(1-q))*(mp-mm)
    return q,W,b


def numerical_reconstruction():
    W,A,L,C,X=model()
    # Root plus equal-energy equation, no global optimizer/exhaustion assertion.
    def equations(z):
        value,grad,_,_=dual(z[:4],z[4],C)
        return np.r_[grad,value]
    sol=root(equations,[1.8199,1.914,1.9146,1.3672,3.71834],tol=1e-11)
    assert sol.success and np.max(np.abs(equations(sol.x)))<1e-10
    s,g=sol.x[:4],sol.x[4]
    value,grad,H,p=dual(s,g,C)
    tangent=np.linalg.qr(np.column_stack([np.ones(12),np.eye(12)[:,1:]]))[0][:,1:]
    primal_H=tangent.T@(np.diag(1/p)-g*X@X.T)@tangent
    A7=X@X.T
    inactive=p-1/12-X@np.linalg.solve(X.T@X,X.T@(p-1/12))
    sad=root(lambda z:dual(z,g,C)[1],[.941,1.001,.962,.686],tol=1e-11)
    assert np.max(np.abs(dual(sad.x,g,C)[1]))<1e-10
    sigma=(2*L[3]*(L[4]+L[5])-L[4]*L[5])/(24*L[3])
    t2=1-sigma/(L[3]/6)
    def scalarS(r):
        return 1-(L[3]/6)*r*(1-r)/(sigma*(1+r)-(L[3]/6)*r*r)-L[6]/(3*sigma)*r/(1+r)**2
    opt=minimize_scalar(scalarS,bounds=(.01,.99),method='bounded',options={'xatol':1e-14})
    s3=math.acosh(1/opt.x)/math.sqrt(L[3]/6)
    q,Wp,b=parity(C,np.array([s3,0,0,0]))
    direct=1-b@np.linalg.solve(sigma*np.eye(4)-Wp,b)
    assert abs(direct-opt.fun)<1e-12

    edges=list(itertools.combinations(range(12),2)); D=np.zeros((66,12))
    for z,(i,j) in enumerate(edges): D[z,i]=1; D[z,j]=-1
    G=np.diag([W[i,j] for i,j in edges])
    tree=G@D@np.linalg.pinv(A)@D.T@G; cyc=G-tree
    E=np.arange(0,12,2); O=np.arange(1,12,2)
    hidden=A[np.ix_(O,O)]; coupling=A[np.ix_(E,O)]
    memory=coupling@np.linalg.solve(hidden,coupling.T)
    schur=A[np.ix_(E,E)]-memory
    full7_p=softmax(2*np.cos(np.pi*np.arange(12)/2))
    full7_cov=moments(X,full7_p)[1]
    return {
        'coexistence_candidate_g':float(g),'candidate_s':s.tolist(),
        'candidate_p':p.tolist(),'stationarity_energy_residual':float(np.max(np.abs(equations(sol.x)))),
        'primal_tangent_min_eigenvalue':float(np.linalg.eigvalsh(primal_H)[0]),
        'inactive_power_fraction':float(inactive@inactive/np.sum((p-1/12)**2)),
        'saddle_s':sad.x.tolist(),'saddle_dual_hessian_eigenvalues':np.linalg.eigvalsh(dual(sad.x,g,C)[2]).tolist(),
        'sigma_candidate':float(sigma),'t_star_squared':float(t2),
        'resolvent_face_r_sech':float(opt.x),'resolvent_face_s3':s3,'resolvent_face_min':float(opt.fun),
        'resolvent_direct_minus_formula':float(direct-opt.fun),
        'wrong_exp_parameter_resolvent':float(scalarS(math.exp(-math.sqrt(L[3]/6)*s3))),
        'cycle_rank':int(np.linalg.matrix_rank(cyc,tol=1e-10)),
        'tree_rank':int(np.linalg.matrix_rank(tree,tol=1e-10)),
        'cycle_divergence_residual':float(np.linalg.norm(D.T@cyc)),
        'tree_contraction_residual':float(np.linalg.norm(D.T@tree@D-A)),
        'memory_trace':float(np.trace(memory)), 'memory_rank':int(np.linalg.matrix_rank(memory,tol=1e-10)),
        'effective_spectrum':np.linalg.eigvalsh(schur).tolist(),
        'rank7_row':A7[0].tolist(),
        'full7_counterexample_hessian_eigenvalues':np.linalg.eigvalsh(np.eye(7)/g-full7_cov).tolist(),
        'scope':'Floating reconstructions; local roots and samples are NOT interval existence or globality proofs.'}


def run():
    return {'manifest':audit_manifest(), 'exact':exact_certificates(), 'numerical':numerical_reconstruction()}


if __name__=='__main__':
    import json
    print(json.dumps(run(),indent=2))
