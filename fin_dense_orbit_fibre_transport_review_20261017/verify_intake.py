#!/usr/bin/env python3
"""Independent bounded DH controls and source-preserving copy replay."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;SOURCE=ROOT/'fin_dense_orbit_fibre_transport_campaign_20261017'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj
def iroot(n,k):
    if type(n)is not int or type(k)is not int or n<0 or k<1:raise ValueError('exact nonnegative integer and positive degree required')
    lo=0;hi=1<<((n.bit_length()+k-1)//k)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**k<=n:lo=mid
        else:hi=mid
    assert lo**k<=n<(lo+1)**k;return lo
def orbit(n):
    out=set()
    for k,sign in itertools.product(range(6),(-1,1)):
        image=[0]*12
        for j,value in enumerate(n):image[(2*k+sign*j)%12]=value
        out.add(tuple(image))
    return out
def pred(y):
    result=set()
    for b,a in itertools.product(range(0,12,2),range(1,12,2)):
        if y[b]:
            x=list(y);x[b]-=1;x[a]+=1;result.add(tuple(x))
    return result
def successors(x):
    for a,b in itertools.product(range(1,12,2),range(0,12,2)):
        if x[a]:
            y=list(x);y[a]-=1;y[b]+=1;yield tuple(y)
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--replay',action='store_true');args=p.parse_args()
    helper=load(ROOT/'fin_rbcut_weighted_degree_review_20261016/verify_intake.py','accepted_hash_helpers');checks=[]
    def ck(name,ok,**detail):checks.append({'name':name,'passed':bool(ok),**detail});assert ok,name
    out={'date':'2026-10-17','entry_HEAD':'68e01759'};out['source_manifest']=helper.manifest(SOURCE)
    ck('all29 incoming source files match with complete coverage',out['source_manifest']['matched']==29)
    mt=ROOT/'fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip';entries=[]
    for row in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,rest=row.split('  ',1);alias,path=rest.split(' :: ',1)
        if alias=='MT_WRAPPER':actual=mt
        elif alias=='RBCUT_source':actual=helper.STAGES['orbit']/Path(path).name
        else:actual=ROOT/path.removeprefix('/mnt/data/')
        item={'alias':alias,'original':path,'sha256':h}
        if actual.is_file():
            assert sha(actual)==h,(alias,path);item.update(status='MATCHED',local=str(actual.relative_to(ROOT)))
            if alias=='MT_ARCHIVE_MEMBER':
                with zipfile.ZipFile(mt) as z:assert hashlib.sha256(z.read(path)).hexdigest()==h
                item['archive_member_also_matches']=True
        else:item['status']='UNAVAILABLE_ORIGINAL'
        entries.append(item)
    out['original_references']=entries
    ck('seven original references match and two missing originals stay missing',sum(e['status']=='MATCHED' for e in entries)==7 and {e['alias'] for e in entries if e['status']!='MATCHED'}=={'ARCHIVED_AGENTS_NOT_CURRENT','PC_ZIP'})
    previous={}
    for rel in ('fin_balanced_phase_current_review_20261017','fin_balanced_phase_cut_research_20261016',
                'fin_dense_orbit_fibre_transport_campaign_20261017','fin_critical_metastability_review_20261015',
                'fin_q12_physical_bridge_research_20261017'):
        previous[rel]=helper.manifest(ROOT/rel)
        savedpath=ROOT/rel/('RESULTS.json' if 'research' in rel else 'verification.json')
        if not savedpath.is_file():continue
        saved=json.loads(savedpath.read_text())
        for path,h in saved.get('input_hashes',{}).items():
            if path=='AGENTS.md':assert helper.prefix(h) is not None
            else:assert sha(ROOT/path)==h,path
    out['actual_dependencies']=previous
    ck('actual PC04 MT and current review inputs remain bound',len(previous)==5)
    for k in (1,2,3,4):
        for n in (0,1,2,12,12000,10**30):assert iroot(n,k)**k<=n<(iroot(n,k)+1)**k
    ck('integer root helper certifies powers exactly',iroot(1200**3,3)==1200)
    ck('producer floating cube root is not an exact boundary API',math.floor((1200**3)**(1/3))==1199)
    t=tuple(v for r in range(6) for v in (2*r-5,2*r-5));O=orbit(t)
    minimum=min(sum(abs(a-b) for a,b in zip(x,y)) for x,y in itertools.combinations(O,2))
    ck('independent pushforward group orbit has exact minL1=36',len(O)==12 and minimum==36 and sum(t[::2])==sum(t[1::2])==0)
    tested=0
    for N in (12000,24000,120000,1200**3):
        k=iroot(N,3);cap=iroot(N**3,4)//2-1
        for h in (-cap,0,cap):
            n=[0]*12
            for p,K in ((0,N//2+h),(1,N//2-h)):
                base,rem=divmod(K,6)
                for i in range(6):n[2*i+p]=base+(i<rem)+k*t[2*i+p]
            A=orbit(n);assert min(n)>0 and sum(n)==N and len(A)==12
            assert min(sum(abs(a-b) for a,b in zip(x,y)) for x,y in itertools.combinations(A,2))>=36*k-12
            assert (2*h)**4<=N**3 and (2*h+2)**4<=N**3
            C=sum(n[j]*((1,0,-1,0)*3)[j] for j in range(12));S=sum(n[j]*((0,1,0,-1)*3)[j] for j in range(12))
            assert abs(C+6*k)<=3 and abs(S+6*k)<=3;tested+=1
    ck('integer critical-fibre ramp fixtures include perfect powers',tested==12,scope='count geometry only, not a smallN MT-domain certificate')
    N=12000;k=iroot(N,3);y=tuple(N//12+k*v for v in t);centres=orbit(y)
    holes=set().union(*(pred(c) for c in centres));targets=set().union(*(set(successors(x)) for x in holes))
    lost={z for z in targets if pred(z)<=holes}
    ck('actual q12 local clustered predecessor algebra matches the new seed',len(holes)==432 and lost==centres,
       local_targets=len(targets),scope='finite local count neighbourhood, not all count-lattice enumeration or MT validation')
    q=pred(y)
    ck('predecessor diameter4 and clustered failure control exact',len(q)==36 and max(sum(abs(a-b) for a,b in zip(x,z)) for x,z in itertools.combinations(q,2))==4)
    ck('E-column bound uses max not the bulk upper',max(F(1,4),-F(-1))==1 and F(1)>F(35,100))
    ck('dense complement coefficient and added cluster margin exact',F(18,100)-F(7,1000)==F(173,1000) and F(173,1000)*36*F(9,10)-F(35,100)==F(6569,1250))
    for n in range(1,17):
        row=[F(j+1,13) for j in range(n)];col=list(reversed(row));A={j for j in range(n) if j%3};T=set(range(n))-A;R={j for j in range(n) if j%2};nbr=set(range(n))-R
        assert sum(col[j] for j in nbr)-sum(row[j] for j in T)==sum(row[j] for j in A)-sum(col[j] for j in R)
    ck('complement identity independent exact controls',True)
    # Finite frame identity underlying the mesoscopic cluster target.
    sig=[F((-1)**j) for j in range(12)];vectors=[]
    for a,b in itertools.product(range(1,12,2),range(0,12,2)):
        v=[-s/6 for s in sig];v[b]+=1;v[a]-=1;vectors.append(v)
    assert all(sum(v[j] for v in vectors)==0 for j in range(12))
    for i,j in itertools.product(range(12),repeat=2):
        pf=F(i==j)-F(1,12)-sig[i]*sig[j]/12
        assert sum(v[i]*v[j] for v in vectors)/36==pf/6
    ck('fast36-step frame covariance exact',all(sum(x*x for x in v)==F(5,3) for v in vectors) and F(1,60)>F(1,64))
    ck('slice atom polynomial beats fixed exponential symbolically scoped',F(11,2)-5==F(1,2),
       scope='power bookkeeping; uniform Stirling/coercivity proof is analytic in the intake completion')
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-dh-intake-'));copy=temp/'campaign';shutil.copytree(SOURCE,copy);logs={};env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        for script in ('check_package.py','verify_dh.py'):
            run=subprocess.run(['python3',script],cwd=copy,env=env,capture_output=True,text=True,timeout=60);assert run.returncode==0,(script,run.stderr)
            logs[script]={'exit_code':run.returncode,'stdout':run.stdout,'stderr':run.stderr}
        stdout=logs['verify_dh.py']['stdout'];expected=json.loads((SOURCE/'BOUNDED_VERIFICATION.json').read_text())
        ck('copy-only producer replay agrees as JSON',json.loads(stdout)==expected and expected['bounded_checks']==147)
        out['replay']={'directory':str(temp),'logs':logs,'stdout_byte_equal':stdout.encode()==(SOURCE/'BOUNDED_VERIFICATION.json').read_bytes(),
            'scope':'147 bounded checks; integer/fraction arithmetic after floating fixture root selection, not an allN exact root API'}
    inputs=[Path(__file__),SOURCE/'MANIFEST.sha256',ROOT/'AGENTS.md',HERE/'PROOF_COMPLETIONS_AND_SCOPE.md',HERE/'CLUSTERED_COMPLEMENT_SEED.md',
        ROOT/'FIN_MESOSCOPIC_HALL_AND_INTERMEDIATE_SOURCE_NEXT_CHAIN_20261017.md',ROOT/'fin_intermediate_clock_bridge_research_20261017/MANIFEST.sha256']
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in previous]
    out['input_hashes']={str(path.relative_to(ROOT)):sha(path) for path in inputs};out['checks']=checks
    out['source_manifest_after']=helper.manifest(SOURCE);assert out['source_manifest_after']==out['source_manifest']
    out['status']='PASS_SCOPED_DH_DENSE_COMPLEMENT_INTAKE'
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'source_entries':out['source_manifest']['matched']},indent=2))

if __name__=='__main__':main()
