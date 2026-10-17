#!/usr/bin/env python3
"""Independent symbolic/bounded audit. Never execute historical source writers."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ARCHIVE = ROOT / 'material_dowodowy/korpus_qw_pozostaly/raporty_json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_math(text, names, functions):
    """Whitelist historical arithmetic; no eval/sympify of source strings."""
    def visit(n):
        if isinstance(n, ast.Constant) and type(n.value) is int:
            return s.Integer(n.value)
        if isinstance(n, ast.Name) and n.id in names:
            return names[n.id]
        if isinstance(n, ast.Tuple):
            return tuple(visit(v) for v in n.elts)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.UAdd, ast.USub)):
            return visit(n.operand) * (-1 if isinstance(n.op, ast.USub) else 1)
        if isinstance(n, ast.BinOp):
            a, b = visit(n.left), visit(n.right)
            for cls, op in ((ast.Add, lambda: a+b), (ast.Sub, lambda: a-b),
                            (ast.Mult, lambda: a*b), (ast.Div, lambda: a/b),
                            (ast.Pow, lambda: a**b)):
                if isinstance(n.op, cls):
                    return op()
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and not n.keywords:
            if n.func.id == 'Derivative':
                return s.diff(*(visit(v) for v in n.args))
            if n.func.id in functions:
                return functions[n.func.id](*(visit(v) for v in n.args))
        raise ValueError('Unlicensed source expression: ' + ast.dump(n)[:160])
    return visit(ast.parse(text, mode='eval').body)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    checks = []

    def ck(name, condition, **scope):
        assert bool(condition), name
        checks.append(dict(name=name, passed=True, **scope))

    paths = [ROOT/'AGENTS.md', ROOT/'TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex']
    paths += [ROOT/f for f in (
        'fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md',
        'fundamental_action_reconstruction/K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md',
        'fundamental_action_reconstruction/F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md',
        'fundamental_action_reconstruction/F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md',
        'fundamental_action_reconstruction/S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md',
        'SUMMARY_GROK.md',
        'fundamental_action_reconstruction/STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md',
        'fundamental_action_reconstruction/generated/p2685_s1635_strict_lagrangian_eom_reverse_closure_obstruction_matrix.md',
        'fundamental_action_reconstruction/generated/p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.md',
        'fundamental_action_reconstruction/generated/p3075_s2025_variational_source_obstruction_table.md')]
    reports = {}
    for n in (2163, 2165, 2166):
        paths += list(ROOT.glob(f'QW_{n}_*.py'))
        rp = next(ARCHIVE.glob(f'report_qw{n}_*.json'))
        paths.append(rp)
        paths += list((ROOT/'material_dowodowy/lean_fin_dowody/lean').glob(f'*QW{n}.lean'))
        reports[n] = json.loads(rp.read_text())
    before = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    expected = json.loads((HERE/'MODEL_AND_SOURCE_CONTRACT.json').read_text())['immutable_AGENTS_sha256']
    ck('AGENTS exact baseline remains preserved', sha(ROOT/'AGENTS.md') == expected)
    ck('historical full-source completion flags remain false', all(
        any(('complete_fin_action_completed' in k) and v is False for k, v in r['flags'].items())
        for r in reports.values()))

    q = list(s.symbols('psi0:12', real=True)) + [s.Symbol('phi', real=True)]
    p = list(s.symbols('p0:13', real=True))
    m = list(s.symbols('m2_psi0:12', real=True)) + [s.Symbol('m2_phi', real=True)]
    g4 = list(s.symbols('g4_psi0:12', real=True))
    g6 = list(s.symbols('g6_psi0:12', real=True))
    y = list(s.symbols('gY0:12', real=True))
    lam = s.Symbol('lam_phi', real=True)
    K = s.Matrix(12, 12, lambda i,j: 0 if i == j else s.Symbol(f'K_{i}_{j}', real=True))
    Ks = (K+K.T)/2
    selfV = [m[i]*q[i]**2/2+g4[i]*q[i]**4/4+g6[i]*q[i]**6/6 for i in range(12)]
    selfV += [m[12]*q[12]**2/2+lam*q[12]**4/4]
    pair = {(i,j): Ks[i,j]*q[i]*q[j] for i in range(12) for j in range(i+1,12)}
    pair.update({(i,12): y[i]*q[i]**2*q[12]**2 for i in range(12)})
    V = sum(selfV)+sum(pair.values())
    grad = [s.diff(V,a) for a in q]
    H = s.hessian(V,q)
    x = s.Symbol('x', real=True)
    functions = {str(a):s.Function(str(a)) for a in q}
    funcs = [functions[str(a)](x) for a in q]
    names = {str(a):a for a in V.free_symbols | {x}}
    dynamic = dict(zip(q,funcs))
    L = sum(s.diff(f,x)**2 for f in funcs)/2 - V.subs(dynamic)
    for n in (2163,2165):
        ck(f'QW{n} archived full Lagrangian agrees independently',
           s.expand(parse_math(reports[n]['model']['lagrangian_density'],names,functions)-L) == 0)
    for i,a in enumerate(q):
        key = 'eom_phi' if i==12 else f'eom_psi{i}'
        eom = s.diff(funcs[i],x,2)+grad[i].subs(dynamic)
        ck('full source equation '+key, s.expand(parse_math(reports[2165]['model'][key],names,functions)-eom)==0)
    ck('QW2166 archived potential agrees independently',
       s.expand(parse_math(reports[2166]['model']['potential_total'],names,{})-V)==0)
    manual = s.zeros(13)
    for i in range(12):
        manual[i,i] = m[i]+3*g4[i]*q[i]**2+5*g6[i]*q[i]**4+2*y[i]*q[12]**2
        manual[i,12] = manual[12,i] = 4*y[i]*q[i]*q[12]
        for j in range(i+1,12): manual[i,j]=manual[j,i]=Ks[i,j]
    manual[12,12] = m[12]+3*lam*q[12]**2+2*sum(y[i]*q[i]**2 for i in range(12))
    ck('all169 Hessian entries match declared formula',H==manual)
    vacuum=list(s.symbols('vpsi0:12',real=True))+[s.Symbol('vphi',real=True)]
    eta=[s.Function(f'eta{i}')(x) for i in range(12)]+[s.Function('eta_phi')(x)]
    alt_names=names | {str(a):a for a in vacuum}
    alt_functions=functions | {str(a.func):a.func for a in eta}
    Hv=H.subs(dict(zip(q,vacuum)))
    for i,label in ((0,'eta0'),(6,'eta6'),(12,'eta_phi')):
        expr=s.diff(eta[i],x,2)+sum(Hv[i,j]*eta[j] for j in range(13))
        ck('archived Hessian action '+label,s.expand(parse_math(reports[2166]['model']['sample_eom_'+label],alt_names,alt_functions)-expr)==0)
    ck('arbitrary Hessian point need not be stationary',s.diff(q[12]**2/2+q[12]**4/4,q[12]).subs(q[12],1)==2)
    Edot = sum((s.diff(V,q[i])*p[i] - p[i]*grad[i]) for i in range(13))
    ck('Hamiltonian derivative vanishes for complete26dim phase state',s.expand(Edot)==0)
    def along(expr):
        return s.expand(sum(s.diff(expr,q[i])*p[i]-s.diff(expr,p[i])*grad[i] for i in range(13)))
    errors=[]
    for a in range(13):
        ea=p[a]**2/2+selfV[a]; currents=0
        for (i,j),U in pair.items():
            if a in (i,j):
                b=j if a==i else i
                ea+=U/2
                currents+=(p[b]*s.diff(U,q[b])-p[a]*s.diff(U,q[a]))/2
        errors.append(along(ea)-currents)
    ck('all13 component balances and antisymmetric pair currents',all(s.expand(e)==0 for e in errors))
    hidden=list(s.symbols('v0:12', real=True)); zero=dict.fromkeys(hidden,0)
    radial={q[i]**2:q[i]**2+hidden[i]**2 for i in range(12)}
    complexV=sum(m[i]*(q[i]**2+hidden[i]**2)/2+g4[i]*(q[i]**2+hidden[i]**2)**2/4+g6[i]*(q[i]**2+hidden[i]**2)**3/6+y[i]*q[12]**2*(q[i]**2+hidden[i]**2) for i in range(12))+selfV[12]
    complexV+=sum(Ks[i,j]*(q[i]*q[j]+hidden[i]*hidden[j]) for i in range(12) for j in range(i+1,12))
    ck('realification restricts to template and imaginary forces vanish',s.expand(complexV.subs(zero)-V)==0 and all(s.diff(complexV,h).subs(zero)==0 for h in hidden))
    ck('imaginary Hermitian mixing can break the real slice',s.diff(q[0]*hidden[1]-hidden[0]*q[1],hidden[1]).subs(zero)==q[0])
    a,b,eps=s.symbols('a b eps',real=True)
    V1=a*a/2+a**4/4; V2=a*a/2+a**4/2
    ck('same positive Hessian different nonlinear forces',s.diff(V1,a,2).subs(a,0)==s.diff(V2,a,2).subs(a,0)==1 and s.diff(V2-V1,a)==a**3)
    U=eps*a*a*b*b
    ck('distinct composition completions preserve isolated axes and zero Hessian',U.subs(a,0)==U.subs(b,0)==0 and s.hessian(U,(a,b)).subs({a:0,b:0})==s.zeros(2) and s.diff(U,a)!=0)
    ck('nonlinear samefield composition residual nonzero',s.expand((a+b)**3-a**3-b**3-3*a*b*(a+b))==0)
    ck('unrestricted negative quartic permits finite travel time',s.factor(a**4/2-a*a-4)==(a-2)*(a+2)*(a*a+2)/2,
       scope='analytic integral endpoint and infinity bounds are in the proof; no simulated blowup')
    z=s.symbols('z',real=True); mu,ll=s.symbols('mu ll',positive=True)
    v=mu/s.sqrt(ll); k=v*s.tanh(z); D=lambda f:mu/s.sqrt(2)*s.diff(f,z)
    simp=lambda f:s.simplify(s.expand_trig(f).rewrite(s.exp))
    ck('existing legacy kink solves spatial not homogeneous-time equation',simp(-D(D(k))-mu**2*k+ll*k**3)==0)
    f0=1/s.cosh(z)**2; f1=s.tanh(z)/s.cosh(z)
    op=lambda f:-s.diff(f,z,2)+(4-6/s.cosh(z)**2)*f
    ck('translation zero mode and shape eigenvalue3',simp(op(f0))==0 and simp(op(f1)-3*f1)==0)
    h=s.Function('h')(z)
    fact=-s.diff(s.diff(h,z)+2*s.tanh(z)*h,z)+2*s.tanh(z)*(s.diff(h,z)+2*s.tanh(z)*h)
    ck('Darboux positive factorization exact',simp(fact-op(h))==0)
    ell=s.symbols('ell',positive=True)
    wave=1/s.cosh(z)**ell
    ck('kernel-sensitive transverse seed groundstate and factorization',
       simp(-s.diff(wave,z,2)-ell*(ell+1)*wave/s.cosh(z)**2+ell**2*wave)==0,
       scope='commuting conditional control, not the actual noncommuting full field theorem')
    t=s.symbols('t',positive=True)
    phiR=v*(s.tanh(z+t)-s.tanh(z-t)-1)
    residual=(-D(D(phiR))-mu**2*phiR+ll*phiR**3).subs(z,0)
    ck('twowall samefield Cauchy datum has nonzero interaction residual',simp(residual-6*mu**2*v*s.tanh(t)*(1-s.tanh(t))**2)==0)
    yy=s.symbols('yy',real=True)
    norm=s.integrate(yy**2,(yy,-1,1))
    i1=s.integrate(yy**6*(1-yy**2),(yy,-1,1))
    # y=sin(theta): exact even beta integral, no floating quadrature.
    i2=s.beta(s.Rational(5,2),s.Rational(3,2)).expand(func=True)
    ck('normalized shape and continuum source integrals exact',norm==s.Rational(2,3) and i1==s.Rational(4,63) and s.simplify(i2-s.pi/16)==0)
    pi_upper=s.integrate(yy**4*(1-yy)**4/(1+yy**2),(yy,0,1))
    ck('positive integral supplies pi upper and rational continuum margin',s.simplify(pi_upper-(s.Rational(22,7)-s.pi))==0 and s.Rational(22,7)**2<10 and s.Rational(4,63)-s.Rational(30,512)==s.Rational(79,16128))
    ck('kink BPS energy exact',s.simplify(s.sqrt(ll/2)*s.integrate(v*v-yy**2,(yy,-v,v))-2*s.sqrt(2)*mu**3/(3*ll))==0)
    # Conservative energy cannot be inferred from reversible equilibrium.
    Q=s.Matrix([[-2,2,0],[1,-4,3],[0,2,-2]])
    weights=s.Matrix([2,4,6]); ee=s.Matrix([1,2,4])
    direct=(ee.T*s.diag(*weights)*(-Q)*ee)[0]
    form=sum(weights[i]*Q[i,j]*(ee[j]-ee[i])**2 for i in range(3) for j in range(3) if i!=j)/2
    ck('finite reversible energy obstruction exact',Q.rank()==2 and (weights.T*Q)==s.zeros(1,3) and direct==form>0,
       scope='independent rational connected chain; general maximum-principle proof is in the text')
    ck('same retained Phi different hidden state different acceleration',s.diff(-2*y[0]*q[12]*q[0]**2,q[0])!=0)
    omega=s.Rational(743,4000);phase=s.Rational(13,80)
    d0=(s.pi/2-phase)/omega
    ck('strict continuous extension firstzero defeats dyadic exact scaling',s.simplify(omega*d0+phase)==s.pi/2 and s.simplify(omega*d0/2+phase)==s.pi/4+phase/2 and 0<phase<s.pi/2,
       scope='allpositive-distance extension and global homogeneity only, not every possible scaling limit')
    for name in ('__import__("os")','psi0.__class__','open("x")','1.0'):
        try: parse_math(name,names,functions)
        except ValueError: pass
        else: raise AssertionError(name)
    ck('archive expression parser rejects execution and floating literals',True)
    after={str(p.relative_to(ROOT)):sha(p) for p in paths}
    ck('all original scientific sources remain bytepreserved',before==after)
    result={'date':'2026-10-17','status':'PASS_BOUNDED_SOURCE_AND_EVOLUTION_AUDIT',
            'checks':checks,'input_hashes':before,'source_hashes_after':after,
            'symbolic_scope':'13 archived EOM, all169 Hessian entries,13 energy balances and exact limited counterexamples; not PDE formalization or physical validation',
            'strict_joint_origin':'NOT_ESTABLISHED','conditional_sector':'legacy Phi4 on supplied1+1 geometry and signs',
            'empirical_evidence':False,'sympy_version':s.__version__}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'bounded_checks':len(checks)},indent=2))


if __name__=='__main__':
    main()
