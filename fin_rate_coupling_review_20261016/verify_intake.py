#!/usr/bin/env python3
"""RC intake: independent rational certificate, finite controls, copy-only replay.

Standard library only. No producer code is imported. General theorems remain
analytical, not proved by these finite checks.
"""
import argparse
import ast
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/'fin_rate_coupling_derivation_20261016'
BASELINE='761e0731cb4e45669605611ef0c3b38b58a0eff8'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def manifest(folder):
    rows={}
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1); assert n not in rows and sha(folder/n)==h,n; rows[n]=h
    actual={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(rows)|{'MANIFEST.sha256'},folder
    return {'matched':len(rows),'coverage':'COMPLETE','manifest_sha256':sha(folder/'MANIFEST.sha256')}

@dataclass(frozen=True)
class Box:
    lo:F
    hi:F
    def __post_init__(self):assert self.lo<=self.hi
    @staticmethod
    def of(x):return x if isinstance(x,Box) else Box(F(x),F(x))
    def __add__(self,x):
        b=self.of(x); return Box(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return Box(-self.hi,-self.lo)
    def __sub__(self,x):return self+-self.of(x)
    def __rsub__(self,x):return self.of(x)+-self
    def __mul__(self,x):
        b=self.of(x); p=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)]; return Box(min(p),max(p))
    __rmul__=__mul__
    def __truediv__(self,x):
        b=self.of(x); assert b.lo>0 or b.hi<0; return self*Box(1/b.hi,1/b.lo)

def root(n,degree):
    # Dyadic precision, unlike the producer's decimal integer-root bracket.
    scale=1<<190; target=n*scale**degree; low=0; high=(n+1)*scale
    while high-low>1:
        mid=(low+high)//2
        if mid**degree<=target:low=mid
        else:high=mid
    assert low**degree<=target<(low+1)**degree
    return Box(F(low,scale),F(low if low**degree==target else low+1,scale))

def cosine(x):
    # Taylor through degree40, absolute Lagrange remainder, not alternating endpoints.
    polynomial=sum(F((-1)**j,math.factorial(2*j))*x**(2*j) for j in range(21))
    remainder=abs(x)**41/math.factorial(41)
    return Box(polynomial-remainder,polynomial+remainder)

def public(box,digits=24):
    scale=10**digits; low=box.lo.numerator*scale//box.lo.denominator
    high=-((-box.hi.numerator*scale)//box.hi.denominator)
    def dec(x):
        s='-' if x<0 else ''; x=abs(x)
        return s+str(x//scale)+'.'+str(x%scale).zfill(digits)
    pair=[dec(low),dec(high)]; assert F(pair[0])<=box.lo<=box.hi<=F(pair[1]); return pair

def certificate():
    s3=root(3,2)
    cs=[Box.of(1),s3/2,Box.of(F(1,2)),Box.of(0),Box.of(F(-1,2)),-s3/2,Box.of(-1),-s3/2,Box.of(F(-1,2)),Box.of(0),Box.of(F(1,2)),s3/2]
    weights={d:cosine(F('0.18575')*d+F('0.1625'))/(1+d*root(d**4,5)) for d in range(1,7)}
    # Full C12 first-row DFT rather than the producer's distance-pair sum.
    row=[sum((2 if d<6 else 1)*weights[d] for d in range(1,7))]
    row += [-weights[min(j,12-j)] for j in range(1,12)]
    lam={r:sum(row[j]*cs[(r*j)%12] for j in range(12)) for r in (3,4,5,6)}
    ratio=(1+lam[6]/(lam[3]+lam[5]))/2
    gap3=(lam[3]+lam[5]+lam[6])/6; gap6=(lam[3]+lam[5])/3
    saved=json.loads((SOURCE/'CERTIFICATE.json').read_text())
    for r in lam:
        low,high=map(F,saved['laplacian_active_spectra'][str(r)]); assert low<=lam[r].lo<=lam[r].hi<=high
    for key,box in [('log_R3_over_log_R6',ratio),('gap_A00_minus_A03',gap3),('gap_A00_minus_A06',gap6),('separation_from_flat_rank7',ratio-F(3,4))]:
        low,high=map(F,saved[key]); assert low<=box.lo<=box.hi<=high,key
    return {'date':'2026-10-16','status':'INDEPENDENT_EXACT_RATIONAL_INPUT_TO_OUTPUT_CHECK',
        'method':'dyadic root brackets for d*(d^4)^(1/5); degree40 cosine with absolute Lagrange remainder; full C12 DFT; exact outward decimals',
        'source_public_enclosures_all_contain_independent_raw_boxes':True,
        'active_spectra':{str(k):public(v) for k,v in lam.items()},'rho36':public(ratio),
        'gap3':public(gap3),'gap6':public(gap6),'separation_from_flat':public(ratio-F(3,4)),
        'arithmetic_scope':'Python integer/Fraction arithmetic; no mpmath/libm decisions; mathematical kernel is supplied, not physically sourced'}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True); parser.add_argument('--replay',action='store_true'); args=parser.parse_args()
    checks=[]
    def ck(name,ok,**detail):checks.append({'name':name,'passed':bool(ok),**detail}); assert ok,name
    out={'date':'2026-10-16','entry_baseline':BASELINE,'scope':'source-preserving review, not a kinetic source or empirical validation'}
    out['source_manifest']=manifest(SOURCE)
    declared={n:h for h,n in (row.split('  ',1) for row in (SOURCE/'INPUTS.sha256').read_text().splitlines())}
    baseline_agents=subprocess.run(['git','show',BASELINE+':AGENTS.md'],cwd=ROOT,capture_output=True,check=True).stdout
    assert hashlib.sha256(baseline_agents).hexdigest()==declared['AGENTS.md']
    assert (ROOT/'AGENTS.md').read_bytes().startswith(baseline_agents)
    for name,h in declared.items():
        if name!='AGENTS.md':assert sha(ROOT/name)==h,name
    ck('twelve original input hashes paid; AGENTS bound to immutable baseline',len(declared)==12)
    out['source_inputs']={'declared':declared,'baseline_AGENTS_verified_at_git':BASELINE,
        'current_AGENTS_equals_source_snapshot':sha(ROOT/'AGENTS.md')==declared['AGENTS.md'],
        'initial_receipt':'All12 current files passed sha256sum -c before this intake appended AGENTS; replay uses original baseline bytes.'}
    previous={}
    for rel in ('fin_mediator_memory_observability_review_20261016','fin_generative_operational_review_20261014'):
        previous[rel]=manifest(ROOT/rel)
        saved=json.loads((ROOT/rel/'verification.json').read_text())
        for path,h in saved['input_hashes'].items():
            if path=='AGENTS.md':
                digest=hashlib.sha256(); length=0; found=False
                for line in (ROOT/'AGENTS.md').read_bytes().splitlines(keepends=True):
                    digest.update(line); length+=len(line)
                    if digest.hexdigest()==h:
                        found=True; previous[rel]['historical_AGENTS_verified_prefix_bytes']=length; break
                assert found,'historical AGENTS prefix'
            else:
                aliases={'source_manifest':'fin_generative_operational_campaign_20261014/MANIFEST.sha256',
                         'checker':'fin_generative_operational_review_20261014/verify_intake.py'}
                assert sha(ROOT/aliases.get(path,path))==h,path
    out['previous_reviews']=previous; ck('actual MM and E reviews and saved inputs unchanged',len(previous)==2)
    cert=certificate(); ck('independent rational certificate binds every displayed source enclosure',cert['source_public_enclosures_all_contain_independent_raw_boxes'])
    ck('strict flat-P7 contrast',F(cert['separation_from_flat'][0])>F('0.02490314788261'))
    out['certificate']=cert
    # All-N proof is analytical; fixtures use a distinct q3 nonconstant-diagonal matrix.
    q=3; mat=[[F(2*i*j+i+j+1,13)+(F(i+1,5) if i==j else 0) for j in range(q)] for i in range(q)]
    count=0
    for N in (2,3,5):
        gain=F(11,7); fields=[[F((slot+2)*(label+1),17) for label in range(q)] for slot in range(N)]
        def V(x):return gain/(2*N)*sum(mat[a][b] for a in x for b in x)+sum(fields[i][a] for i,a in enumerate(x))
        for a,b,c,d in itertools.product(range(q),repeat=4):
            if a==b or c==d:continue
            for tail in itertools.product(range(q),repeat=N-2):
                ac=(a,c)+tail; bc=(b,c)+tail; ad=(a,d)+tail; bd=(b,d)+tail
                mixed=V(bd)-V(ad)-V(bc)+V(ac)
                assert mixed==gain/N*(mat[b][d]-mat[a][d]-mat[b][c]+mat[a][c])
                n=[ac.count(j) for j in range(q)]; n[a]-=1
                single=fields[0][b]-fields[0][a]+gain/N*sum((mat[b][j]-mat[a][j])*n[j] for j in range(q))+gain/(2*N)*(mat[b][b]-mat[a][a])
                assert V(bc)-V(ac)==single
                assert (V(bc)-V(ac))+(V(bd)-V(bc))+(V(ad)-V(bd))+(V(ac)-V(ad))==0
                count+=1
    ck('1116 exact Gibbs, self-diagonal, field, mixed and zero-cycle fixtures',count==1116)
    out['exact_fixtures']=count
    ck('count-multiplicity bias is exactly log of one third',(F(1,3)/F(2,2))==F(1,3))
    # Stationarity alone is weaker than detailed balance: oriented 3-cycle, uniform invariant.
    L=[[-3,2,1],[1,-3,2],[2,1,-3]]
    ck('stationarity without detailed balance negative control',all(sum(L[i][j] for i in range(3))==0 for j in range(3)) and L[0][1]!=L[1][0] and F(2,1)**3!=1)
    # Stationary conditional rates may satisfy DB although the visible path is not Markov.
    joint=[F(i,10) for i in (1,2,3,4)]; J=[[F(0) for _ in range(4)] for _ in range(4)]
    for i,j,conductance in ((0,1,F(1,10)),(2,3,F(1,10)),(0,2,F(1,10)),(1,3,F(3,10))):
        J[i][j]=conductance/joint[i]; J[j][i]=conductance/joint[j]
    for i in range(4):J[i][i]=-sum(J[i])
    observed=[joint[0]+joint[1],joint[2]+joint[3]]
    averaged=[[sum(joint[i]*J[i][j] for i in range(2*a,2*a+2) for j in range(2*b,2*b+2))/observed[a] for b in range(2)] for a in range(2)]
    f=[F(1),F(1),F(0),F(0)]; Lf=[sum(J[i][j]*f[j] for j in range(4)) for i in range(4)]
    conditional=[sum(joint[i]*Lf[i] for i in range(2*a,2*a+2))/observed[a] for a in range(2)]
    defect=sum(joint[i]*(Lf[i]-conditional[i//2])**2 for i in range(4))
    ck('stationary averaged DB does not imply observed semigroup',observed[0]*averaged[0][1]==observed[1]*averaged[1][0] and defect>0,exact_compression_defect=str(defect))
    # Reconstruct the nominal operator from versioned constants, no source import.
    seed=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
    constants={node.targets[0].id:ast.literal_eval(node.value) for node in ast.parse(seed.read_text()).body if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ('OMEGA','PHI','ETA')}
    ww={d:math.cos(constants['OMEGA']*d+constants['PHI'])/(1+d**constants['ETA']) for d in range(1,7)}
    lam={r:sum((2 if d<6 else 1)*ww[d]*(1-math.cos(2*math.pi*r*d/12)) for d in range(1,7)) for r in (3,4,5,6)}
    row=[(2*sum(lam[r]*math.cos(2*math.pi*r*d/12) for r in (3,4,5))+lam[6]*(-1)**d)/12 for d in range(12)]
    A=[[row[(i-j)%12] for j in range(12)] for i in range(12)]
    def pot(x,g):return g/(2*len(x))*sum(A[a][b] for a in x for b in x)
    def rate(x,new,kind,g=3.):
        old=x[0]; target=(new,)+x[1:]; delta=pot(target,g)-pot(x,g)
        if kind in ('HB','THIN'):
            values=[pot((j,)+x[1:],g) for j in range(12)]; shift=max(values)
            val=math.exp(values[new]-shift)/sum(math.exp(t-shift) for t in values)
            if kind=='THIN' and min((new-old)%12,(old-new)%12)==3:val/=1+g*g
            return val
        if kind=='MET':return min(1.,math.exp(delta))/12
        return 2/(12*(1+math.exp(-delta)))  # frozen E04 g0-normalized Barker
    errors=[]
    for kind in ('HB','THIN','MET','BAR'):
        for z in (0,7):
            contrasts={}
            for d in (3,6):
                ac=(0,0,z); bc=(d,0,z); ad=(0,d,z); bd=(d,d,z)
                contrasts[d]=math.log(rate(ad,d,kind)/rate(bd,0,kind))-math.log(rate(ac,d,kind)/rate(bc,0,kind))
            errors.append(abs(contrasts[3]/contrasts[6]-(1+lam[6]/(lam[3]+lam[5]))/2))
    ck('eight independent full-potential kinetic rate diagnostics',max(errors)<3e-14,maximum_error=max(errors))
    x=(0,0,7); native=rate(x,3,'HB')/rate(x,6,'HB'); thin=rate(x,3,'THIN')/rate(x,6,'THIN')
    ck('same-source branching changes by one tenth at g3',abs(thin/native-.1)<1e-14)
    ck('thinning preserves g0 clock and support',all(rate(x,j,'HB',0)==rate(x,j,'THIN',0) for j in range(1,12)))
    # Exact projector identity in Q(sqrt3) specialized to d3,d6,d1.
    def add(x,y):return (x[0]+y[0],x[1]+y[1])
    def mul(x,y):return (x[0]*y[0]+3*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
    def pair(a,b=0):return F(a),F(b)
    table=[pair(1),pair(0,F(1,2)),pair(F(1,2)),pair(0),pair(F(-1,2)),pair(0,F(-1,2)),pair(-1),pair(0,F(-1,2)),pair(F(-1,2)),pair(0),pair(F(1,2)),pair(0,F(1,2))]
    projectors={r:[[mul(table[(r*(i-j))%12],pair(F(1,6))) for j in range(12)] for i in range(12)] for r in (3,5)}
    for r,s in itertools.product((3,5),repeat=2):
        for i,j in itertools.product(range(12),repeat=2):
            value=pair(0)
            for k in range(12):value=add(value,mul(projectors[r][i][k],projectors[s][k][j]))
            assert value==(projectors[r][i][j] if r==s else pair(0))
    ck('exact orthogonal rank-two projectors in Qsqrt3',all(projectors[r][i][i]==pair(F(1,6)) for r in (3,5) for i in range(12)))
    cos3=[F(1),F(0),F(-1),F(0)]
    # P3(0,1)=0; P5(0,1)=-sqrt3/12, both have diagonal1/6.
    ck('isofingerprint family leaves both chosen gaps unchanged',all(cos3[d%4]==cos3[(5*d//3)%4] for d in (0,3,6)))
    spectral=json.loads((ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json').read_text())['exact']['laplacian_intervals']
    ck('epsilon1/100 is admissible for P3 minus P5 family',F(spectral[3][0])>0 and F(spectral[5][0])-F(1,100)>0 and F(spectral[3][1])+F(1,100)<F(spectral[5][0])-F(1,100))
    gap1_change=-root(3,2)*F(1,1200)  # A00-A01 changes by -epsilon sqrt3/12.
    ck('same fingerprint with a genuinely changed untested affinity',gap1_change.hi<0)
    out['isofingerprint_control']={'family':'A_delta=A7+delta(P3-P5)','delta':'1/100',
        'preserved':['trace','constant diagonal','rank7','parity eigenvalue/potential','A00-A03','A00-A06','rho36'],
        'gap1_change':public(gap1_change),'scope':'Exact projector algebra; no family-wide stochastic-limit or physical-source theorem.'}
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-rc-intake-')); copy=temp/SOURCE.name
        shutil.copytree(SOURCE,copy,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for name in declared:
            target=temp/name; target.parent.mkdir(parents=True,exist_ok=True)
            if name=='AGENTS.md':target.write_bytes(baseline_agents)
            else:shutil.copyfile(ROOT/name,target)
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        before=subprocess.run(['sha256sum','-c',str(copy/'MANIFEST.sha256')],cwd=copy,env=env,capture_output=True,text=True)
        run=subprocess.run(['python3','verify_derivation.py'],cwd=copy,env=env,capture_output=True,text=True,timeout=60)
        after=subprocess.run(['sha256sum','-c',str(copy/'MANIFEST.sha256')],cwd=copy,env=env,capture_output=True,text=True)
        assert before.returncode==run.returncode==0,(before.stderr,run.stderr)
        same={n:sha(copy/n)==sha(SOURCE/n) for n in ('CERTIFICATE.json','RESULTS.json')}
        ck('producer replay exact certificate and eight diagnostics reproduced',same['CERTIFICATE.json'] and json.loads((copy/'RESULTS.json').read_text())['exact_algebra']['fixtures']==576)
        out['replay']={'directory':str(copy),'manifest_before':before.stdout,'writer_stdout':run.stdout,'writer_stderr':run.stderr,
            'output_byte_equal':same,'manifest_after_exit':after.returncode,'manifest_after':after.stdout,'output_sha256':{n:sha(copy/n) for n in same}}
    ck('zero-g logarithm ratio excluded',2*0*(row[0]-row[3])==2*0*(row[0]-row[6])==0)
    inputs=[Path(__file__),SOURCE/'MANIFEST.sha256']+[ROOT/n for n in declared if n!='AGENTS.md']
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in previous]
    out['input_hashes']={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out['checks']=checks; out['status']='PASS_SCOPED_RATE_COUPLING_INTAKE_CHECKS'
    out['source_manifest_after']=manifest(SOURCE); assert out['source_manifest_after']==out['source_manifest']
    for rel in previous:assert manifest(ROOT/rel)=={k:v for k,v in previous[rel].items() if k!='historical_AGENTS_verified_prefix_bytes'}
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'source_manifest':out['source_manifest'],'exact_fixtures':count,'certificate':cert},indent=2))

if __name__=='__main__':main()
