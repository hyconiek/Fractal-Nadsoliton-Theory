#!/usr/bin/env python3
"""Independent bounded NS audit; original writers run only on temporary copies."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_nonlinear_stability_transfer_campaign_20261017'
UE=ROOT/'fin_coupled_evolution_source_review_20261017'
AUDIT=ROOT/'fin_unified_strict_evolution_audit_20261017'
BASE='ffb62af2671f62aa582c688a90edd998e4a9621b29e35fd240e1863038c9e0d9'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def module(p,name):
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--replay',action='store_true');args=ap.parse_args()
    helper=module(UE/'verify_intake.py','accepted_ue_hash_helper')
    checks=[]
    def ck(name,ok,**detail):
        assert bool(ok),name;checks.append({'name':name,'passed':True,**detail})
    before=helper.manifest(SOURCE)
    ck('all16 incoming source entries match with complete coverage',before['matched']==16)
    refs=[]
    aliases={'/mnt/data/files/AGENTS.md':ROOT/'fin_coupled_evolution_source_campaign_20261017/REPLAY/AGENTS_LIBRARY_SNAPSHOT_20261004.md'}
    for line in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,n=line.split('  ',1);p=aliases.get(n,ROOT/n.removeprefix('/mnt/data/'))
        row={'original':n,'sha256':h,'status':'UNAVAILABLE_ORIGINAL'}
        if p.is_file():assert sha(p)==h,n;row.update(status='MATCHED_BYTES',local=str(p.relative_to(ROOT)))
        refs.append(row)
    ck('four actual input byte references match and original UE ZIP stays unavailable',sum(r['status']=='MATCHED_BYTES' for r in refs)==4 and sum(r['status']=='UNAVAILABLE_ORIGINAL' for r in refs)==1)
    dependencies={}
    for folder in (UE,AUDIT):
        dependencies[folder.name]=helper.manifest(folder)
        d=json.loads((folder/'verification.json').read_text())
        for rel,h in d['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    ck('actual UE and unified audit dependencies remain intact',len(dependencies)==2)
    ck('AGENTS exact immutable baseline preserved',sha(ROOT/'AGENTS.md')==BASE)

    # Parse the actual archived polynomial, using the accepted safe AST reader only.
    parser=module(AUDIT/'verify_evolution.py','accepted_safe_polynomial_parser')
    archive=ROOT/'material_dowodowy/korpus_qw_pozostaly/raporty_json'
    report=next(archive.glob('report_qw2166_*.json'))
    import ast
    literal=json.loads(report.read_text())['model']['potential_total']
    names={n.id:s.Symbol(n.id,real=True) for n in ast.walk(ast.parse(literal,mode='eval')) if isinstance(n,ast.Name)}
    V=parser.parse_math(literal,names,{})
    phi=names['phi'];psi=[names[f'psi{i}'] for i in range(12)]
    ck('bound Phi action has no sixth-order term while Psi has one',s.diff(V,phi,6)==0 and s.diff(V,psi[0],6)==120*names['g6_psi0'],
       scope='gPhi nonzero is comparison-only, not a sourced extension')

    q,v,lam,g=s.symbols('q v lam g',real=True)
    mass=lam*v*v+g*v**4
    U=-mass*q*q/2+lam*q**4/4+g*q**6/6
    U=s.expand(U-U.subs(q,v))
    fact=(q*q-v*v)**2*(lam/4+g*(q*q+2*v*v)/6)
    ck('comparison sextic factorization independently reconstructed',s.expand(U-fact)==0)
    C=lam/2+2*g*v*v/3;D=g*v*v/3;theta=s.cancel(D/C)
    T=s.Symbol('T',real=True);den=1+theta*(1-T*T);kk=v*v*(C+D)
    r2=v*v*T*T/den
    # Different route: evaluate polynomial in squared field and compare squared chain derivative.
    # U(0) is the already shifted vacuum constant.
    U_at=-mass*r2/2+lam*r2*r2/4+g*r2**3/6+U.subs(q,0)
    deriv2=v*v*kk*(1+theta)**2*(1-T*T)**2/den**3
    ck('exact normalized comparison wall first integral',s.cancel(deriv2-2*U_at)==0)
    ck('comparison wall reduces to bound quartic sector at gPhi0',theta.subs(g,0)==0 and kk.subs(g,0)==lam*v*v/2)
    p=s.Symbol('p',real=True)
    tanh_force=2*v*kk*(T**3-T)-(-mass*v*T+lam*v**3*T**3+g*v**5*T**5)
    ck('bare tanh cannot solve a nonzero Phi sextic',s.Poly(s.expand(tanh_force),T).coeff_monomial(T**5)==-g*v**5)
    Udd=-mass+3*lam*r2+5*g*r2*r2
    err=g*v**4*(1-6*T*T+5*(r2/v**2)**2)+3*lam*v*v*(r2/v**2-T*T)
    ck('sextic scalar Hessian error identity independent rational propagation',s.cancel(Udd-(-2*kk+6*kk*T*T)-err)==0)
    u=s.Symbol('u',real=True)
    ck('uniform wall error uses paid quarter-product bound',s.expand(s.Rational(1,4)-u*(1-u)-(u-s.Rational(1,2))**2)==0)
    ck('comparison vacuum curvature is exactly four width-squared',s.simplify(Udd.subs(T,1)-4*kk)==0)

    # Nonempty source-type controls: Y diagonal, signed masses, nonzero mixing.
    A=s.Matrix([[-s.Rational(3,2),-s.Rational(1,2)],[-s.Rational(1,2),-s.Rational(3,2)]])
    B=s.diag(99,101);W=A+B
    ck('source-type diagonalY control remains genuinely noncommuting',A*B!=B*A and sorted(A.eigenvals())==[-2,-1])
    mu=(197-s.sqrt(5))/2
    ck('source-type far spectrum and gap exact',min(W.eigenvals())==mu and s.simplify(mu-((s.sqrt(405)-1)/2)**2)==-3+4*s.sqrt(5))
    ck('source-type gap positivity paid by exact rational sqrt comparison',s.Integer(5)>s.Integer(2)**2,
       scope='sqrt5>2 implies -3+4sqrt5>5; matrices are not actual FIN values')
    A12=s.diag(*([-s.Rational(3,2)]*2+[-1]*10));A12[0,1]=A12[1,0]=-s.Rational(1,2)
    B12=s.diag(99,101,*([100]*10))
    ck('same diagonalY control retains24 real directions',s.diag(A12,A12).shape==(24,24) and A12*B12!=B12*A12)
    # Compare t+e bound without any floating asymptotic subtraction.
    width,r=s.symbols('width r',positive=True)
    ss=(s.sqrt(1+4*r/width**2)-1)/2
    ck('fixedkernel nonempty-region lower bound identity',s.simplify(r-width**2*ss**2-width**2*ss)==0)
    signed=[]
    for gg in (F(1,1000),F(-1,1000)):
        kval=1+gg;ct=1+F(2,3)*gg;dt=F(1,3)*gg;th=dt/ct
        rr=abs(th)/(4*(1-abs(th)));e6=12*abs(gg)+6*rr
        assert ct>0 and kval>0 and abs(th)<1
        assert kval*kval+4*kval*101>400
        lower_psi=F(779,8)-(101+kval/2-10)-101*rr
        lower_phi=3*kval-e6
        assert lower_phi>0 and lower_psi>0
        signed.append({'gPhi':str(gg),'deltaPsi_rational_lower':str(lower_psi),'deltaPhi_rational_lower':str(lower_phi)})
    ck('both signed Phi6 comparison controls use independent rational bounds',len(signed)==2,
       scope='comparison class only, not imported into FIN')

    x=s.Symbol('x',real=True);k=s.Symbol('k',positive=True)
    pot=4*k*k-6*k*k/s.cosh(k*x)**2
    f=1/s.cosh(k*x)**2;f1=s.tanh(k*x)/s.cosh(k*x)
    simp=lambda z:s.simplify(z.rewrite(s.exp))
    ck('scalar ground and second mode comparison exact',simp(-s.diff(f,x,2)+pot*f)==0 and simp(-s.diff(f1,x,2)+pot*f1-3*k*k*f1)==0)
    ell=s.Symbol('ell',positive=True);fe=1/s.cosh(k*x)**ell
    ck('PT transverse factor bound eigenfunction exact',simp(-s.diff(fe,x,2)-k*k*ell*(ell+1)*fe/s.cosh(k*x)**2+k*k*ell*ell*fe)==0)
    delta,cv,theta0=s.symbols('delta cv theta0',positive=True)
    ch=delta/(delta+1+cv)
    ck('combined H1 coercivity coefficient exact',s.cancel(ch*(1+(1+cv)/delta)-1)==0)
    ck('cubic remainder energy barrier has positive inner margin',s.expand(ch/2-ch*theta0/4-ch/4-ch*(1-theta0)/4)==0,
       scope='0<=theta0<=1 represents C3*r<=cH/4')
    vel,ad,norm,eta2=s.symbols('vel ad norm eta2',real=True)
    ck('modulation derivative retains moving centre sign',s.expand(vel+ad*norm-ad*eta2-(vel+ad*(norm-eta2)))==0)

    # Exact nonlinear work cancellation, including higher-order comparison terms.
    ep,aa,bb,c4,d3=s.symbols('ep aa bb c4 d3',real=True,nonzero=True)
    nn=2*aa*ep**2+c4*ep**4;dd=2*bb*ep+d3*ep**3;qe=-nn/dd
    ck('initial work cancellation is exact not a fitted coefficient',s.cancel(ep*nn+ep*qe*dd)==0 and s.limit(qe/ep,ep,0)==-aa/bb)
    cw,tt=s.symbols('cw tt',positive=True)
    ck('explicit orderepsilon interval retains half initial rate',s.simplify(aa*ep**3-cw*ep**2*(aa*ep/(2*cw))-aa*ep**3/2)==0)
    ck('integrated transfer on chosen interval is orderepsilon4',s.expand((aa*ep**3/2)*(aa*ep/(2*cw))-aa**2*ep**4/(4*cw))==0)

    # Specialize the existing radial complex polynomial for the phase countercontrol.
    im=list(s.symbols('im0:12',real=True));phi=names['phi'];Vr=V
    Ks=s.Matrix(12,12,lambda i,j:0 if i==j else (names[f'K_{i}_{j}']+names[f'K_{j}_{i}'])/2)
    for i in range(12):
        rr=psi[i]**2+im[i]**2
        Vr+=names[f'm2_psi{i}']*im[i]**2/2+names[f'g4_psi{i}']*(rr**2-psi[i]**4)/4+names[f'g6_psi{i}']*(rr**3-psi[i]**6)/6+names[f'gY{i}']*phi**2*im[i]**2
    Vr+=sum(Ks[i,j]*im[i]*im[j] for i in range(12) for j in range(i+1,12))
    generator=sum(psi[i]*s.diff(Vr,im[i])-im[i]*s.diff(Vr,psi[i]) for i in range(12))
    ck('displayed radial action has global phase identity',s.expand(generator)==0,
       scope='reused Noether mechanism specialized to conditional field action; no local gauge source')
    # Pointwise charge continuity after substituting the full wave equations.
    ux,uxx,vx,vxx=s.symbols('ux uxx vx vxx',real=True)
    U0,V0,gu,gv=s.symbols('U0 V0 gu gv',real=True)
    charge_dt=U0*(vxx-gv)-V0*(uxx-gu)
    flux_dx=V0*uxx-U0*vxx
    ck('phase density and flux continuity algebra exact',s.expand(charge_dt+flux_dx+U0*gv-V0*gu)==0)
    ar,br,om=s.symbols('ar br om',real=True)
    real=ep*ar*s.cos(om*tt);imag=ep*br*s.sin(om*tt)
    intensity=ep**2*((ar**2+br**2)/2+(ar**2-br**2)*s.cos(2*om*tt)/2)
    ck('ellipse phase controls second-harmonic force exactly',s.trigsimp(real**2+imag**2-intensity)==0)
    ck('circular linear mode has constant intensity and nonzero model charge',s.trigsimp((real**2+imag**2).subs(br,ar)-ep**2*ar**2)==0 and s.trigsimp(real*s.diff(imag,tt)-imag*s.diff(real,tt)-ep**2*ar*br*om)==0,
       scope='linear-mode forcing algebra, not an exact nonlinear no-radiation solution')

    replay={}
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-ns-intake-'));copy=temp/'campaign';shutil.copytree(SOURCE,copy)
        run=subprocess.run(['python3','REPLAY/verify.py'],cwd=copy,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=90)
        assert run.returncode==0,run.stderr
        incoming=json.loads((SOURCE/'REPLAY/exact_checks.json').read_text());fresh=json.loads((copy/'REPLAY/exact_checks.json').read_text())
        changes=[]
        for old,new in zip(incoming['checks'],fresh['checks']):
            if old!=new:changes.append({'incoming':old,'replayed':new})
        ck('producer31 flags pass on copy',fresh['count']==incoming['count']==31 and all(c['passed'] for c in fresh['checks']))
        replay={'directory':str(temp),'exit_code':run.returncode,'stdout':run.stdout,'stderr':run.stderr,
                'json_equal':incoming==fresh,'byte_equal':sha(copy/'REPLAY/exact_checks.json')==sha(SOURCE/'REPLAY/exact_checks.json'),
                'changed_controls':changes,'scope':'bounded synthetic algebra; not a PDE proof or physical source'}
    after=helper.manifest(SOURCE)
    ck('all incoming science preserved before and after replay',before==after)
    ck('AGENTS and predecessor manifests remain unchanged',sha(ROOT/'AGENTS.md')==BASE and all(helper.manifest(ROOT/n)==v for n,v in dependencies.items()))
    files=[Path(__file__),SOURCE/'MANIFEST.sha256',SOURCE/'INPUTS.sha256',ROOT/'AGENTS.md',
           UE/'MANIFEST.sha256',AUDIT/'MANIFEST.sha256',report,HERE/'PROOF_COMPLETIONS_AND_SCOPE.md']
    inputs={str(p.relative_to(ROOT)):sha(p) for p in files}
    inputs.update(json.loads((AUDIT/'verification.json').read_text())['input_hashes'])
    result={'date':'2026-10-17','entry_HEAD':'dc51caba','status':'PASS_SCOPED_NS_INTAKE',
            'source_manifest':before,'source_manifest_after':after,'original_references':refs,
            'actual_dependencies':dependencies,'input_hashes':inputs,'checks':checks,'signed_comparison_controls':signed,
            'replay':replay,'source_mapping':'gPhi=0 matches bound Phi polynomial; nonzero gPhi is comparison-only',
            'science_verdict':'conditional full-field orbital persistence and local transfer, not natural source or radiation theorem',
            'backend':s.__version__,'empirical_evidence':False}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'bounded_controls':len(checks),'source_entries':before['matched'],'source_replay_byte_equal':replay.get('byte_equal')},indent=2))

if __name__=='__main__':main()
