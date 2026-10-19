#!/usr/bin/env python3
"""Independent RP exact controls and separately labelled numerical quadrature."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_coherent_mode_radiation_campaign_20261017'
NS=ROOT/'fin_nonlinear_stability_transfer_review_20261017'
UE=ROOT/'fin_coupled_evolution_source_review_20261017'
AUDIT=ROOT/'fin_unified_strict_evolution_audit_20261017'
POLICY=ROOT/'fin_agents_update_20261019'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mod(p,name):
    spec=importlib.util.spec_from_file_location(name,p);obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--replay',action='store_true')
    ap.add_argument('--progress',action='store_true');args=ap.parse_args()
    policy=mod(POLICY/'verify_update.py','accepted_instruction_prefix_tools')
    contract=json.loads((POLICY/'UPDATE_CONTRACT.json').read_text())
    agents=(ROOT/'AGENTS.md').read_bytes();tail=policy.check_prefix(agents,contract)
    checks=[]
    def ck(name,ok,**detail):
        assert bool(ok),name;checks.append({'name':name,'passed':True,**detail})
        if args.progress:print(name,file=sys.stderr,flush=True)
    before=policy.manifest(SOURCE)
    ck('all19 source entries have exact hashes and complete coverage',before['entries']==19)
    aliases={'Library:NS-source:HANDOFF(20261017-101344).md':'HANDOFF.md',
             'Library:NS-source:NS_ORBITAL_THEOREM.md':'NS_ORBITAL_THEOREM.md',
             'Library:NS-source:NS_TRANSFER_THEOREM.md':'NS_TRANSFER_THEOREM.md',
             'Library:NS-source:REPORT.md':'REPORT.md'}
    refs=[]
    for row in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        if not row or row.startswith('#'):continue
        h,n=row.split('  ',1);p=ROOT/'fin_nonlinear_stability_transfer_campaign_20261017'/aliases[n]
        assert sha(p)==h,n;refs.append({'reference':n,'local':str(p.relative_to(ROOT)),'sha256':h})
    ck('four Library payload hashes match actual local NS documents',len(refs)==4,
       scope='matching bytes, not Library custody or historical web-access attestation')
    deps={}
    for folder in (NS,UE,AUDIT,POLICY):deps[folder.name]=policy.manifest(folder)
    for folder in (NS,UE,AUDIT):
        saved=json.loads((folder/'verification.json').read_text())
        for rel,h in saved['input_hashes'].items():
            if rel=='AGENTS.md':assert hashlib.sha256(agents[:contract['base_prefix_bytes']]).hexdigest()==h
            else:assert sha(ROOT/rel)==h,rel
    ck('current actual scientific dependencies intact with historical AGENTS distinguished',len(deps)==4)
    ck('live instruction revision prefix independently verified',hashlib.sha256(agents[:contract['updated_prefix_bytes']]).hexdigest()==contract['updated_prefix_sha256'])

    # Exact rational contract arithmetic, no binary64 endpoint assumptions.
    ck('all declared masses negative uniformly',-F(1,4)+F(1,5)<0)
    low=F(399,200);high=F(401,200)
    ck('complement scalar root margin proved by squaring',low*(low+1)<F(299,50) and F(2,25)+low-2==F(3,40))
    ck('full NS comparison gap has strict rational margin',high*(high+1)>F(301,50) and F(21,4)-high**2>F(49,40))
    ck('complement vacuum threshold above active edge',F(1,10)+F(299,50)-F(1,50)==F(303,50) and F(303,50)>6)
    d2,d3,h=s.symbols('d2 d3 h',real=True)
    aa=s.Matrix([[0,h],[h,0]]);yy=s.diag(3+d2,3+d3)
    ck('noncommutation is explicitly in the spectator block',s.expand((aa*yy-yy*aa)[0,1]-h*(d3-d2))==0)

    x=s.Symbol('x',real=True);q=s.Symbol('q',real=True);a=s.Symbol('a',real=True)
    simp=lambda e:s.simplify(e.rewrite(s.exp))
    tt=s.tanh(x);ss=1/s.cosh(x)
    f=s.sqrt(s.Rational(3,4))*ss**2
    ck('active ground mode exact at every parameter value',simp(-s.diff(f,x,2)+(a+6*tt**2)*f-(a+2)*f)==0)
    y=s.Symbol('y',real=True)
    ck('mode normalization exact polynomial integral',s.Rational(3,4)*s.integrate(1-y*y,(y,-1,1))==1)
    ck('quartic source normalization retains DC and halfcosine factor',simp(tt*3*f*f-s.Rational(9,4)*tt*ss**4)==0)
    ck('frequency is strictly above scalar threshold over entire interval',4*(F(5,4))-4==1 and 4*(F(7,4))-4==3)

    D=lambda n,z:s.diff(z,x)+n*tt*z
    Dt=lambda n,z:-s.diff(z,x)+n*tt*z
    T=lambda z:Dt(2,Dt(1,z))
    Tadj=lambda z:D(1,D(2,z))
    test=s.Function('w')(x)
    ck('Darboux normalization polynomial exact',simp(Tadj(T(test))-(s.diff(test,x,4)-5*s.diff(test,x,2)+4*test))==0)
    L=lambda z:-s.diff(z,x,2)+(4-6*ss**2)*z
    ck('intertwining with shifted free operator exact',simp(L(T(test))-T(-s.diff(test,x,2)+4*test))==0)
    plane=s.exp(s.I*q*x);numerator=3*tt**2-3*s.I*q*tt-(q*q+1)
    ck('Jost numerator derived independently by operator composition',simp(T(plane)-plane*numerator)==0)
    psi=plane*numerator/((1-s.I*q)*(2-s.I*q))
    ck('Jost ODE at arbitrary real q exact',simp(L(psi)-(q*q+4)*psi)==0)
    ck('Jost denominator has the required unitary modulus',s.expand((1-s.I*q)*(2-s.I*q)*(1+s.I*q)*(2+s.I*q)-(1+q*q)*(4+q*q))==0)
    f0=ss**2;f1=ss*tt
    ck('both discrete states removed by continuous transform adjoint',simp(Tadj(f0))==0 and simp(Tadj(f1))==0)

    nn=s.Symbol('n',positive=True,integer=True)
    ck('Fourier recurrence follows from exact derivative identity',simp(s.diff(ss**nn,x,2)-nn**2*ss**nn+nn*(nn+1)*ss**(nn+2))==0)
    F4=s.Symbol('F4',real=True)
    F6=(q*q+16)*F4/20
    iq4=-s.I*q*F4/4;iq6=-s.I*q*F6/6
    raw=s.Rational(9,4)*(3*(iq4-iq6)+3*s.I*q*(F4-F6)-(q*q+1)*iq4)
    den=(1+s.I*q)*(2+s.I*q)
    M=s.Rational(9,4)*s.I*q*(q*q+4)*F4/(8*den)
    ck('on-shell amplitude independently reduced from Fourier moments',s.cancel(raw/den-M)==0)
    magnitude=s.cancel(M*s.conjugate(M))
    formula=81*q*(q*q+4)*F4**2/(2048*(1+q*q))
    ck('two spectral roots and positive retarded sign give stated Gamma',s.cancel(magnitude/(2*q)-formula)==0)
    ck('one-root normalization mistake is rejected',s.simplify((magnitude/(4*q)-formula).subs({q:1,F4:1}))!=0)
    pilower=4*sum((F((-1)**j,2*j+1) for j in range(8)),F(0))
    ck('uniform elementary bounds for pi and q',pilower>3 and F(7,4)**2>3 and F(22,7)*F(7,4)/2<3)
    lower_F=F(3*5)/(6*F(27,2));pref=F(81*5,2048*4)
    lower_gamma=pref*lower_F**2
    ck('independent uniform positive Gamma rational certificate',lower_F==F(5,27) and lower_gamma==F(125,73728) and lower_gamma>0,
       scope='analytic exp/pi inequalities paid in the proof completion; not a finite-N or physical rate')

    # Actual source channel invariance, not an invented modal truncation.
    parser=mod(AUDIT/'verify_evolution.py','accepted_safe_archive_parser')
    report=next((ROOT/'material_dowodowy/korpus_qw_pozostaly/raporty_json').glob('report_qw2166_*.json'))
    lit=json.loads(report.read_text())['model']['potential_total']
    names={n.id:s.Symbol(n.id,real=True) for n in ast.walk(ast.parse(lit,mode='eval')) if isinstance(n,ast.Name)}
    V=parser.parse_math(lit,names,{})
    pp=[names[f'psi{i}'] for i in range(12)];phi=names['phi'];zero_other=dict.fromkeys(pp[1:],0)
    sparse={names[f'K_{i}_{j}']:0 for i in range(12) for j in range(12) if i!=j and (i,j) not in ((1,2),(2,1))}
    vactive=s.expand(V.subs(zero_other).subs(sparse))
    ck('active reduced polynomial independent of spectator kernel',s.diff(vactive,names['K_1_2'])==0 and s.diff(vactive,names['K_2_1'])==0)
    ck('all spectator forces vanish on the prepared species subspace',all(s.diff(V,p).subs(zero_other).subs(sparse)==0 for p in pp[1:]))
    active=vactive.subs({names['m2_phi']:-2,names['lam_phi']:2,names['m2_psi0']:a,names['gY0']:3})
    g4=names['g4_psi0'];g6=names['g6_psi0'];p0=pp[0]
    force_phi=2*phi*(phi**2-1)+6*phi*p0**2
    force_psi=(a+6*phi**2)*p0+g4*p0**3+g6*p0**5
    ck('prepared two-field equations derived from actual polynomial',s.expand(s.diff(active,phi)-force_phi)==0 and s.expand(s.diff(active,p0)-force_psi)==0)
    kap,eta=s.symbols('kap eta',real=True)
    phi_error=force_phi.subs(phi,kap+eta)-force_phi.subs({phi:kap,p0:0})-(6*kap**2-2)*eta
    psi_error=force_psi.subs(phi,kap+eta)-(a+6*kap**2)*p0
    ck('exact neutral remainder vertices and signs independently expanded',
       s.expand(phi_error-(6*kap*p0**2+6*eta*p0**2+6*kap*eta**2+2*eta**3))==0 and
       s.expand(psi_error-(12*kap*eta*p0+6*eta**2*p0+g4*p0**3+g6*p0**5))==0)
    ck('Phi6 is still absent and exact neutral parity is preserved',s.diff(V,phi,6)==0 and s.expand(V.subs(phi,-phi)-V)==0 and s.expand(s.diff(V,phi).subs(phi,-phi)+s.diff(V,phi))==0,
       scope='with reflection-symmetric preparation, uniqueness preserves Phi odd and Psi even; not an external pin')

    # Independent finite-interval quadrature, not the producer infinite-line route.
    import mpmath as mp
    mp.mp.dps=45;num=[]
    for qq in (mp.mpf(1),mp.sqrt(2),mp.sqrt(3)):
        dd=(1+1j*qq)*(2+1j*qq)
        def integrand(yy):
            if abs(yy)==1:return mp.mpc(0)
            return mp.mpf(9)/4*yy*(1-yy*yy)*mp.exp(-1j*qq*mp.atanh(yy))*(3*yy*yy+3j*qq*yy-(qq*qq+1))/dd
        numeric=mp.quad(integrand,[-1,-mp.mpf('.5'),0,mp.mpf('.5'),1])
        ff=mp.pi*qq*(qq*qq+4)/(6*mp.sinh(mp.pi*qq/2))
        exact=mp.mpf(9)/4*1j*qq*(qq*qq+4)*ff/(8*dd)
        error=abs(numeric-exact);assert error<mp.mpf('1e-30')
        num.append({'q':mp.nstr(qq,25),'max_formula_error':mp.nstr(error,8),
                    'Gamma_point':mp.nstr(abs(numeric)**2/(2*qq),25)})
    ck('independent three finite-interval quadratures agree numerically',len(num)==3,
       scope='mpmath point diagnostic, not an outward certificate or PDE theorem')

    replay={}
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-rp-intake-'));copy=temp/'campaign';shutil.copytree(SOURCE,copy)
        run=subprocess.run(['python3','REPLAY/verify.py'],cwd=copy,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=90)
        assert run.returncode==0,run.stderr
        ck('producer16 numerical controls pass on temporary copy','PASS 16/16 bounded controls' in run.stdout)
        replay={'directory':str(temp),'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr,
                'scope':'40digit point quadrature/differentiation with tolerances, not exact algebra or outward/PDE proof'}
    after=policy.manifest(SOURCE)
    ck('source and prior scientific packages remain byte-preserved',before==after and all(policy.manifest(ROOT/n)==v for n,v in deps.items()))
    ck('live AGENTS not changed by scientific replay',(ROOT/'AGENTS.md').read_bytes()==agents)
    files=[Path(__file__),SOURCE/'MANIFEST.sha256',SOURCE/'INPUTS.sha256',POLICY/'UPDATE_CONTRACT.json',report,
           HERE/'PROOF_COMPLETIONS_AND_SCOPE.md',ROOT/'FIN_COHERENT_MODE_RADIATION_NEXT_CHAIN_20261017.md']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in files}
    hashes.update({str(p.relative_to(ROOT)):sha(p) for p in [ROOT/r['local'] for r in refs]})
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True,check=True).stdout.strip()
    result={'date':'2026-10-19','entry_HEAD':'4c1faefa','verification_HEAD':head,'status':'PASS_SCOPED_RP_ON_SHELL_INTAKE',
            'source_manifest':before,'source_manifest_after':after,'matched_Library_payloads':refs,
            'dependencies':deps,'input_hashes':hashes,
            'AGENTS_context':{'bytes':len(agents),'sha256':hashlib.sha256(agents).hexdigest(),
                              'prior_verified_prefix_bytes':contract['updated_prefix_bytes'],
                              'prior_verified_prefix_sha256':contract['updated_prefix_sha256'],
                              'later_bytes':tail,'scope':'actual live instruction context for this intake; no old fullhash substitution'},
            'checks':checks,'numeric_crosscheck':num,'replay':replay,
            'uniform_Gamma_lower_rational':str(lower_gamma),'RP004':'OPEN_NONLINEAR_OUTWARD_FLUX',
            'kernel_scope':'active e1 is spectator-kernel-blind in this constrained family',
            'formal_or_empirical_validation':False}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'source_entries':before['entries'],'bounded_checks':len(checks),'RP004':result['RP004']},indent=2))

if __name__=='__main__':main()
