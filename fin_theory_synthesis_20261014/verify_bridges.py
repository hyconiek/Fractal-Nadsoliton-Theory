#!/usr/bin/env python3
"""Bounded diagnostics for the intrinsic-theory synthesis, not a campaign.
Exact statements are proved in CORE_DERIVED_BRIDGES.md. Floating checks below
do not turn a proposed stochastic limit into a theorem or physical evidence.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,re
from fractions import Fraction
import numpy as np
from scipy.special import softmax

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p,name):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def manifest(folder):
    entries=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,p=line.split(None,1);assert sha(folder/p)==h,p;entries.append(p)
    return {'files':len(entries),'mismatches':0}

def main():
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--tool-replay',type=Path)
    ns=a.parse_args()
    core=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
    m=load(core,'seed');_,_,lam,X,A=m.build_rank7()
    def H(p,xi,g):
        q=softmax(g*A@p)
        return float(np.dot(p,np.exp(-xi))*np.dot(q,np.exp(xi))-1)
    p=np.arange(1,13,dtype=float);p/=p.sum();xi=np.linspace(-.3,.2,12)
    checks=[]
    for g in (1.,3.7,5.145228719489142):
        q=softmax(g*A@p);force=np.log(p)-g*A@p
        symmetry=abs(H(p,xi,g)-H(p,force-xi,g));hj=abs(H(p,force,g))
        assert symmetry<1e-12 and hj<1e-12
        theta=g*X.T@p;drift=g*(X.T@q-X.T@p)
        gradient=theta/g-X.T@softmax(X@theta)
        graderr=float(np.max(np.abs(drift+g*gradient)));assert graderr<1e-13
        D=np.diag(p)+np.diag(q)-np.outer(p,q)-np.outer(q,p)
        eig=float(np.linalg.eigvalsh(D).min());assert eig>-1e-12
        checks.append({'g':g,'Hamiltonian_reversal_error':symmetry,'Hamilton_Jacobi_zero_error':hj,'retained_gradient_error':graderr,'diffusion_min_eigenvalue_numeric':eig})
    h2=np.array([2,1,-1,-2,-1,1,2,1,-1,-2,-1,1],dtype=int)
    nplus=10+h2;nminus=10-h2;N=120
    assert nplus.sum()==nminus.sum()==N and min(nplus.min(),nminus.min())>0
    pplus=nplus/N;pminus=nminus/N
    means=float(max(np.max(np.abs(X.T@pplus)),np.max(np.abs(X.T@pminus))))
    assert means<1e-14
    # Exact g0 count refresh: q_j=1/12; applied to squared c5 mean.
    c5=X[:,4]
    def square_drift(counts):
        value=0.
        for old in range(12):
            for new in range(12):
                value+=counts[old]/12*((c5[new]-c5[old])/N)**2
        return value
    difference=square_drift(nplus)-square_drift(nminus)
    predicted=lam[5]/(60*N)
    assert abs(difference-predicted)<1e-14 and difference>0
    basis=np.vstack([np.ones(12),X.T]+[(X[:,i]*X[:,j])[None,:] for i in range(7) for j in range(i,7)])
    rank=int(np.linalg.matrix_rank(basis,tol=1e-10));assert rank==12
    text=(ROOT/'GEN001-gen956.md').read_text();headings=[]
    for line_no,line in enumerate(text.splitlines(),1):
        mm=re.match(r'^#{1,6}\s+.*?GEN[- ](\d{1,3})(.*)',line)
        if mm:headings.append({'id':int(mm.group(1)),'line':line_no,'heading':line})
    out={'date':'2026-10-14','baseline':'1b7cd447','source_hashes':{str(core.relative_to(ROOT)):sha(core),'GEN_narrative':sha(ROOT/'GEN001-gen956.md'),'checker':sha(Path(__file__)),'AGENTS.md':sha(ROOT/'AGENTS.md')},
         'F_review_integrity':manifest(ROOT/'fin_generative_calibration_review_20261014'),
         'F_source_integrity':manifest(ROOT/'fin_generative_calibration_campaign_20261014'),
         'READY_integrity':manifest(ROOT/'fin_evidence_intake_tools_20261014'),
         'bounded_numeric_checks':checks,
         'exact_fibre_witness':{'N':N,'n_plus':nplus.tolist(),'n_minus':nminus.tolist(),'mean_max_numeric_residual':means,'generator_square_difference_numeric':difference,'predicted_numeric':predicted,'exact_formula':'Lambda5/(60*120)','gain':0,'coordinate':'retained cosine k5'},
         'moment_lift_rank_numeric':rank,
         'GEN_index':{'lines':len(text.splitlines()),'heading_occurrences':len(headings),'distinct_base_heading_ids':len({h['id'] for h in headings}),'status':'narrative index, not number of proofs'},
         'scope':['proofs are separate analytic finite-model statements','full path LDP and critical stochastic limit are future tasks','no physical records/source/unit result']}
    if ns.tool_replay:
        out['READY_replay']={'directory':str(ns.tool_replay),'unit':json.loads((ns.tool_replay/'UNIT_TESTS.json').read_text()),'independent':json.loads((ns.tool_replay/'INDEPENDENT_VERIFICATION.json').read_text())}
    library=json.loads((Path(__file__).parent/'LITERATURE_REGISTER.json').read_text())
    captures=[]
    for source in library['sources']:
        p=ROOT/'.firecrawl/fin_theory_20261014'/source['cache']
        assert p.is_file(),source['id']
        captures.append({'id':source['id'],'primary_url':source['url'],'cache_sha256':sha(p),'cache_lines':len(p.read_text().splitlines()),'read_scope':source['read_scope']})
    out['literature_access_record']=captures
    ns.output.write_text(json.dumps(out,indent=2)+'\n')
    ns.output.with_name('GEN_HEADINGS_INDEX.json').write_text(json.dumps(headings,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'F_manifest':out['F_review_integrity'],'READY_manifest':out['READY_integrity'],'gain_checks':len(checks),'fibre_difference':difference,'moment_rank':rank,'GEN_index':out['GEN_index']},indent=2))

if __name__=='__main__':main()
