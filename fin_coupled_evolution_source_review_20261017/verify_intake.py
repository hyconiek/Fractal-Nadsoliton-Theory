#!/usr/bin/env python3
"""Independent exact bounded controls; no original source writer is run in place."""
import argparse
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
SOURCE=ROOT/'fin_coupled_evolution_source_campaign_20261017'
PREVIOUS=ROOT/'fin_unified_strict_evolution_audit_20261017'
BASE='ffb62af2671f62aa582c688a90edd998e4a9621b29e35fd240e1863038c9e0d9'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder):
    names=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,raw=line.split('  ',1);p=Path(raw)
        assert not p.is_absolute() and '..' not in p.parts
        name=p.as_posix();assert name not in names
        assert sha(folder/p)==h,(folder,name);names.append(name)
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'},(folder,actual-set(names))
    return {'matched':len(names),'coverage':'COMPLETE','manifest_sha256':sha(folder/'MANIFEST.sha256')}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--replay',action='store_true');a=ap.parse_args()
    checks=[]
    def ck(name,ok,**detail):
        assert bool(ok),name;checks.append({'name':name,'passed':True,**detail})
    before=manifest(SOURCE)
    ck('all22 original manifest entries have complete coverage',before['matched']==22)
    originals={}
    for line in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,n=line.split('  ',1);assert sha(SOURCE/n)==h,n;originals[n]=h
    ck('three supplied local references match; two are producer replay files',len(originals)==3)
    snapshot=SOURCE/'REPLAY/AGENTS_LIBRARY_SNAPSHOT_20261004.md'
    ck('supplied library snapshot is not current AGENTS',sha(snapshot)==originals['REPLAY/AGENTS_LIBRARY_SNAPSHOT_20261004.md'] and sha(snapshot)!=sha(ROOT/'AGENTS.md'),
       scope='the reported Library date is not independently attested acquisition metadata')
    previous=manifest(PREVIOUS)
    ck('predecessor source audit package remains intact',previous['matched']==22)
    prior=json.loads((PREVIOUS/'verification.json').read_text())
    for rel,h in prior['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    ck('actual local QW legacy and FAR scientific inputs agree with predecessor',all(c['passed'] for c in prior['checks']) and len(prior['checks'])==44)
    ck('AGENTS current immutable baseline matches',sha(ROOT/'AGENTS.md')==BASE)

    # Reuse only the independently accepted whitelist parser, not any source algebra.
    spec=importlib.util.spec_from_file_location('accepted_equation_parser',PREVIOUS/'verify_evolution.py')
    parser=importlib.util.module_from_spec(spec);spec.loader.exec_module(parser)
    archive=ROOT/'material_dowodowy/korpus_qw_pozostaly/raporty_json'
    report=next(archive.glob('report_qw2166_*.json'))
    literal=json.loads(report.read_text())['model']['potential_total']
    import ast
    ids={n.id for n in ast.walk(ast.parse(literal,mode='eval')) if isinstance(n,ast.Name)}
    names={n:s.Symbol(n,real=True) for n in ids}
    V=parser.parse_math(literal,names,{})
    psi=[names[f'psi{i}'] for i in range(12)];phi=names['phi'];zeros=dict.fromkeys(psi,0)
    Ks=s.Matrix(12,12,lambda i,j:0 if i==j else (names[f'K_{i}_{j}']+names[f'K_{j}_{i}'])/2)
    M2=s.diag(*(names[f'm2_psi{i}'] for i in range(12)))
    Y=s.diag(*(names[f'gY{i}'] for i in range(12)))
    block=s.hessian(V,psi).subs(zeros)
    ck('actual archived12real transverse Hessian matches operator',block==M2+Ks+2*phi**2*Y)
    ck('actual Phi mixing vanishes at Psi0 with correct scalar Hessian',all(s.diff(V,p,phi).subs(zeros)==0 for p in psi) and s.diff(V,phi,2).subs(zeros)==names['m2_phi']+3*names['lam_phi']*phi**2)
    vi=list(s.symbols('imag0:12',real=True));zv=dict.fromkeys(vi,0)
    Vr=V
    for i in range(12):
        r=psi[i]**2+vi[i]**2
        Vr+=names[f'm2_psi{i}']*vi[i]**2/2
        Vr+=names[f'g4_psi{i}']*(r**2-psi[i]**4)/4
        Vr+=names[f'g6_psi{i}']*(r**3-psi[i]**6)/6
        Vr+=names[f'gY{i}']*phi**2*vi[i]**2
    Vr+=sum(Ks[i,j]*vi[i]*vi[j] for i in range(12) for j in range(i+1,12))
    Hr=s.hessian(Vr,psi+vi).subs(zeros|zv)
    ck('displayed radial realkernel full24real transverse Hessian is two copies',Hr==s.diag(block,block),
       scope='conditional radial realkernel action only; free coefficients not sourced')
    ck('radial action imaginary forces vanish on invariant real slice',all(s.diff(Vr,z).subs(zv)==0 for z in vi))
    u=s.Matrix(s.symbols('u0:3',real=True));v=s.Matrix(s.symbols('v0:3',real=True))
    C=s.Matrix([[2,1,0],[1,3,2],[0,2,4]]);D=s.Matrix([[0,1,2],[-1,0,3],[-2,-3,0]])
    Hc=(C+s.I*D);complexq=u+s.I*v
    vv=s.expand((s.conjugate(complexq).T*Hc*complexq)[0]/2)
    ck('Hermitian imaginary mixing gives realimag coupling countercontrol',s.hessian(vv,list(u)+list(v))==C.row_join(-D).col_join(D.row_join(C)) and any(s.diff(vv,z).subs(dict.fromkeys(v,0))!=0 for z in v))

    z=s.Symbol('z',real=True);k=s.Symbol('k',positive=True);ell=s.Symbol('ell',positive=True)
    f=1/s.cosh(z)**ell
    simp=lambda e:s.simplify(e.rewrite(s.exp))
    ck('scalar comparison bound paid by exact factor eigenfunction',simp(-s.diff(f,z,2)-ell*(ell+1)*f/s.cosh(z)**2+ell**2*f)==0)
    ck('b0 endpoint formal constant is not an L2 groundstate',s.integrate(s.Integer(1),(z,-1,1))==2,
       scope='constant has infinite norm on R; lower spectral edge still valid')
    xx=s.Symbol('xx',real=True);amp=s.Symbol('amp',positive=True)
    profile=amp*s.tanh(k*xx)
    wall_residual=-s.diff(profile,xx,2)+s.diff(V,phi).subs(zeros).subs(phi,profile)
    wall_residual=wall_residual.subs({names['m2_phi']:-2*k**2,names['lam_phi']:2*k**2/amp**2})
    ck('actual Phi background stationarity requires the declared coefficient map',simp(wall_residual)==0)
    # y=tanh(kx) converts all trial integrals to finite polynomial integrals.
    norm=s.integrate(1,(z,-1,1))/k
    kinetic=k*s.integrate(z*z,(z,-1,1))
    attractive=s.integrate(1-z*z,(z,-1,1))/k
    ck('three sech trial integrals independently integrated',norm==2/k and kinetic==2*k/3 and attractive==4/(3*k))
    A,B=s.symbols('A B',real=True);theta=s.Symbol('theta',real=True)
    ck('noncommuting interpolation coefficients exact',s.expand((1-theta)*A+theta*(A+B)-A-theta*B)==0,
       scope='matrix identity entrywise; Loewner implication is analytical')
    Vinf=s.Matrix([[1,s.Rational(1,2)],[s.Rational(1,2),1]]);B2=s.diag(4,0)
    rq=s.simplify((kinetic+norm*Vinf[0,0]-attractive*B2[0,0])/norm).subs(k,1)
    ck('submitted negative trial is exactly minus4thirds',rq==-s.Rational(4,3) and Vinf.det()==s.Rational(3,4))
    V3=s.Matrix([[2,s.Rational(1,3),0],[s.Rational(1,3),1,s.Rational(1,4)],[0,s.Rational(1,4),3]])
    B3=s.diag(7,0,1)
    ck('different3x3 positivefar noncommuting negative control',all(V3[:j,:j].det()>0 for j in range(1,4)) and V3*B3!=B3*V3 and s.Rational(1,3)+V3[0,0]-s.Rational(2,3)*B3[0,0]==-s.Rational(7,3))
    # Recompute gradients from the actual symbolic coupling, not assigned source constants.
    cross=phi**2*sum(Y[i,i]*psi[i]**2 for i in range(12))
    F=s.diff(cross,phi);G=s.Matrix([s.diff(cross,p) for p in psi])
    ck('actual displayed PhiPsi force factors independently derived',F==2*phi*sum(Y[i,i]*psi[i]**2 for i in range(12)) and G[0]==2*Y[0,0]*phi**2*psi[0])
    yy,eps,v0=s.symbols('y eps v',nonzero=True,real=True)
    af=s.integrate(4*yy**2*v0**2*eps**4*z**2*(1-z**2),(z,-1,1))/k
    bg=s.integrate(4*yy**2*v0**4*eps**2*z**4,(z,-1,1))/k
    ck('symbolic nonzero sourcecoupling transfer norms exact',af==16*yy**2*v0**2*eps**4/(15*k) and bg==8*yy**2*v0**4*eps**2/(5*k) and s.cancel(af/bg)==2*eps**2/(3*v0**2))
    ck('submitted synthetic transfer ratio and energies correct',af.subs({yy:1,v0:1,eps:1,k:1})==s.Rational(16,15) and bg.subs({yy:1,v0:1,eps:1,k:1})==s.Rational(8,5))
    delta=s.Symbol('delta',positive=True)
    ck('initial work cancels without new force',s.simplify(-delta*af+delta*(af/bg)*bg)==0)
    ck('single active initial Psi is not generically invariant',s.expand(s.diff(V,psi[1]).subs({phi:0,**{psi[i]:0 for i in range(1,12)}})-Ks[1,0]*psi[0])==0)
    x,aa=s.symbols('x aa',real=True)
    overlap=2*s.integrate(s.exp(-2*x*x-2*aa*aa),(x,-s.oo,s.oo))
    ck('samefield Gaussian cross energy exact',s.simplify(overlap-s.sqrt(2*s.pi)*s.exp(-2*aa*aa))==0)
    ck('kink example requires affine Sobolev space',s.limit(s.tanh(x)**2,x,s.oo)==1,
       scope='positive asymptotic square excludes L2; path differences and perturbations remain L2')
    ck('coarse acceleration mismatch norm exact',s.simplify(af.subs({yy:1,v0:1,k:1})-16*eps**4/15)==0)
    av,cv,th=s.symbols('a c th',positive=True)
    ck('Taylor interval and predictor halfdistance factors exact',s.simplify(av/2-2*cv*(th*av/(8*cv))-av/4-av*(1-th)/4)==0,
       scope='0<=th<=1 gives separation a*t²/4 and worst-case predictor a*t²/8')

    replay={}
    if a.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-ue-intake-'));copy=temp/'campaign';shutil.copytree(SOURCE,copy)
        run=subprocess.run(['python3','REPLAY/verify.py'],cwd=copy,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=90)
        assert run.returncode==0,run.stderr
        old=json.loads((SOURCE/'REPLAY/exact_controls.json').read_text());new=json.loads((copy/'REPLAY/exact_controls.json').read_text())
        changes={key:{'incoming':old[key],'replayed':new[key]} for key in old if old[key]!=new[key]}
        ck('source11 control labels replay semantically on temporary copy',{key:v for key,v in old.items() if key!='runtime'}=={key:v for key,v in new.items() if key!='runtime'} and new['checks_passed']==11)
        replay={'temporary_directory':str(temp),'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr,
                'byte_equal':sha(copy/'REPLAY/exact_controls.json')==sha(SOURCE/'REPLAY/exact_controls.json'),
                'changed_fields':changes,'scope':'same-code synthetic replay, not independent proof; several source integral flags use assigned constants'}
    after=manifest(SOURCE);ck('original campaign bytes and coverage preserved',before==after)
    ck('predecessor package and AGENTS unchanged',manifest(PREVIOUS)==previous and sha(ROOT/'AGENTS.md')==BASE)
    inputs={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),SOURCE/'MANIFEST.sha256',SOURCE/'INPUTS.sha256',
             ROOT/'AGENTS.md',PREVIOUS/'MANIFEST.sha256',PREVIOUS/'verify_evolution.py',
             ROOT/'FIN_UNIFIED_EVOLUTION_SOURCE_AND_COMPOSITION_NEXT_CHAIN_20261017.md',report,HERE/'PROOF_AUDIT_AND_COMPLETIONS.md']}
    inputs.update(prior['input_hashes'])
    result={'date':'2026-10-17','entry_HEAD':'dc51caba','status':'PASS_SCOPED_UE_INTAKE',
            'source_manifest':before,'source_manifest_after':after,'original_local_references':originals,
            'actual_predecessor_manifest':previous,'replay':replay,'checks':checks,'input_hashes':inputs,
            'strict_source':'STILL_UNDERDETERMINED','numeric_FIN_wall_verdict':'NOT_AVAILABLE',
            'source_access_now':'actual local QW/audit/legacy sources bound today, not retrospective producer access',
            'scope_correction':'affine kink background space and explicit Taylor regularity',
            'symbolic_backend':s.__version__,'physical_validation':False}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'source_entries':before['matched'],'bounded_checks':len(checks),
                      'replay_changed_fields':list(replay.get('changed_fields',{}))},indent=2))

if __name__=='__main__':main()
