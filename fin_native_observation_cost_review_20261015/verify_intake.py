#!/usr/bin/env python3
"""NOC intake: source replay and independent small identities, not asymptotic proof."""
import argparse
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
import numpy as np
from scipy.linalg import expm
from scipy.special import softmax

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/"fin_native_observation_cost_campaign_20261015"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder,excluded=()):
    names=[]
    for row in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=row.split("  ",1); assert name not in names and sha(folder/name)==h,name; names.append(name)
    actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}
    extra=actual-set(names)-{"MANIFEST.sha256"}-set(excluded); assert not extra,extra
    return {"matched":len(names),"mismatches":0,"uncovered":sorted(extra),"separately_bound":{n:sha(folder/n) for n in excluded}}
def integer_root(n,power):
    lo=0; hi=1<<((n.bit_length()+power-1)//power)
    while hi-lo>1:
        mid=(hi+lo)//2
        if mid**power<=n:lo=mid
        else:hi=mid
    return lo

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--replay",action="store_true"); args=parser.parse_args()
    checks=[]
    def ck(name,ok,**detail):checks.append({"name":name,"passed":bool(ok),**detail}); assert ok,name
    out={"date":"2026-10-15","baseline":"380e9dc0 plus preserved LI/KO/DKO working-tree additions"}
    out["source_manifest"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"])
    records=json.loads((SOURCE/"SOURCE_INPUTS.json").read_text())
    mapping={}
    for row in records:
        path=row["source_path"]
        if "/fin_dynamic_kernel_observability_campaign_20261015/" in path:
            rel="fin_dynamic_kernel_observability_campaign_20261015/"+path.split("/fin_dynamic_kernel_observability_campaign_20261015/",1)[1]
        elif path.endswith("/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip"):
            rel="fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip"
        else:continue
        assert sha(ROOT/rel)==row["sha256"] and (ROOT/rel).stat().st_size==row["bytes"],rel
        mapping[path]=rel
    expected="".join(r["sha256"]+"  "+r["source_path"]+"\n" for r in records)
    ck("source JSON and absolute-path hash ledger agree",(SOURCE/"SOURCE_INPUTS.sha256").read_text()==expected and (SOURCE/"INPUTS.sha256").read_text()==expected)
    for i in range(1,6):assert (SOURCE/f"NOC-{i:03d}/INPUTS.sha256").read_text()==expected
    out["external_inputs"]={"matched_files_or_wrapper":mapping,"unavailable_original_files":[r["source_path"] for r in records if r["source_path"] not in mapping]}
    out["stage_reference_counts"]={"matched":5*len(mapping),"unavailable_original":5*(len(records)-len(mapping))}
    prior={}
    for rel in ("fin_dynamic_kernel_observability_review_20261015","fin_kernel_observability_review_20261015","fin_critical_metastability_review_20261015"):
        folder=ROOT/rel; item=manifest(folder); old=json.loads((folder/"verification.json").read_text())
        for p,h in old["input_hashes"].items():assert sha(ROOT/p)==h,p
        item["saved_inputs_unchanged"]=len(old["input_hashes"]); prior[rel]=item
    out["actual_review_dependencies"]=prior
    spec=importlib.util.spec_from_file_location("versioned_seed",SEED); seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,_,A=seed.build_rank7(); g=12/lam[6]; gamma=1-lam[3]/lam[6]
    f=np.sqrt(2)*np.array([1.,0.,-1.,0.]*3)
    def direct(labels,k):
        N=len(labels); n=np.bincount(labels,minlength=12); alpha=k/N
        coeff=np.c_[np.full(N,1/math.sqrt(N)),(np.arange(N)<k).astype(float)/math.sqrt(k)-alpha/math.sqrt(k)]
        drift=np.zeros(2); B=np.zeros((2,2)); third=0.
        for a,old in enumerate(labels):
            n0=n.copy(); n0[old]-=1; p=softmax(g/N*(A@n0))
            for j in range(12):
                jump=(f[j]-f[old])*coeff[a]
                drift+=p[j]*jump; B+=p[j]*np.outer(jump,jump); third+=p[j]*np.linalg.norm(jump)**3
        R=float(coeff[:,1]@f[labels]); target=softmax(g/N*(A@n))@f
        forcing=sum(coeff[a,1]*(softmax(g/N*(A@(n-np.eye(12,dtype=int)[old])) )@f-target) for a,old in enumerate(labels))
        return drift,B,third,R,forcing
    rng=np.random.default_rng(20261015); fixtures=[]
    for N,k in ((24,1),(24,8),(120,12),(120,60)):
        labels=rng.integers(0,12,N); drift,B,third,R,forcing=direct(labels,k)
        ck(f"normalized private drift from actual slots N{N} k{k}",abs(drift[1]+R-forcing)<2e-13)
        ck(f"joint native bracket PSD and finite third N{N} k{k}",np.linalg.eigvalsh(B).min()>-1e-13 and math.isfinite(third) and third>0)
        ck(f"alpha-free coefficient cancellation N{N} k{k}",sum(abs(F(int(a<k))-F(k,N)) for a in range(N))==2*F(k)*(1-F(k,N)))
        fixtures.append({"N":N,"k":k,"R_forcing":forcing,"FF_FR_RR":B.tolist(),"third":third,"scope":"finite algebra, not primary critical-window execution"})
    out["actual_generator_fixtures"]=fixtures
    # Exact nonuniform-population conditional private variance (sqrt2 cancels in products).
    population=[0,0,1,2,3,3,4,5,6,7,8,9]; values=[F((1,0,-1,0)[j%4]) for j in population]
    mean=sum(values)/12; variance=sum((x-mean)**2 for x in values)/12
    for k in (1,4,6):
        empirical=sum((sum(values[j] for j in subset)-k*mean)**2*F(2,k) for subset in itertools.combinations(range(12),k))/math.comb(12,k)
        ck(f"exact private second moment without a CLT k{k}",empirical==F(12-k,11)*2*variance)
    # Full finite N2 semigroup / OU Duhamel identity; no Markov assumption for Phi.
    states=list(itertools.product(range(12),repeat=2)); L=np.zeros((144,144)); FF=[]; RR=[]
    for ix,(a,b) in enumerate(states):
        FF.append((f[a]+f[b])/math.sqrt(2)); RR.append((f[a]-f[b])/2)
        for slot,(old,other) in enumerate(((a,b),(b,a))):
            p=softmax(g/2*A[:,other])
            for j in range(12):
                if j==old:continue
                dest=12*j+b if slot==0 else 12*a+j
                L[ix,dest]+=p[j]; L[ix,ix]-=p[j]
    FF=np.array(FF); RR=np.array(RR); alpha=.5; h=.2
    def ou_values(t):
        af=math.sqrt(alpha)*math.exp(-gamma*t); ar=math.exp(-t)
        amp=math.exp(-.5*(alpha/gamma*(1-math.exp(-2*gamma*t))+(1-alpha)*(1-math.exp(-2*t))))
        phase=af*FF+ar*RR
        val=amp*np.cos(phase)
        Gval=amp*((gamma*FF*af+RR*ar)*np.sin(phase)-(af*af+(1-alpha)*ar*ar)*np.cos(phase))
        return val,Gval
    v0,_=ou_values(0); vh,_=ou_values(h); defect=expm(h*L)@v0-vh
    nodes,weights=np.polynomial.legendre.leggauss(14); integ=np.zeros(144)
    for node,weight in zip(nodes,weights):
        t=h*(node+1)/2; v,Gv=ou_values(t)
        integ+=h*weight/2*(expm((h-t)*L)@(L@v-Gv))
    ck("independent finite-state Duhamel diagnostic",np.max(abs(defect-integ))<3e-13 and np.max(abs(defect))>1e-4,residual=float(np.max(abs(defect-integ))))
    ck("full native generator row sums",np.max(abs(L@np.ones(144)))<2e-14)
    # The same(F,R) can have different future observed drift; no quotient is assumed.
    drift_a=float(softmax(g/3*A[:,[1,1]].sum(axis=1))@np.cos(f)-math.cos(f[0]))
    drift_b=float(softmax(g/3*A[:,[1,3]].sum(axis=1))@np.cos(f)-math.cos(f[0]))
    ck("same-chart different hidden-state drift diagnostic",abs(drift_a-drift_b)>.1,
       states=[[0,1,1],[0,1,3]],generator_difference=drift_a-drift_b,
       scope="floating finite witness only, not an independent interval no-go certificate")
    floors=[]
    for N in (120000,1200000,12000000):
        k=integer_root(N**15,16); M=integer_root(N,4)
        ck(f"exact integer floors N{N}",k**16<=N**15<(k+1)**16 and M**4<=N<(M+1)**4)
        floors.append({"N":N,"k_exact":k,"M_exact":M,"k_le_half":k<=N//2})
    out["floor_fixtures"]=floors
    ck("frozen polynomial error exponents",-F(1,4)+2*F(1,16)==-F(1,8) and -F(1,4)+F(1,16)==-F(3,16) and -F(15,32)+F(1,16)==-F(13,32))
    from certify_mean_gap import build
    cert=build(); assert json.loads((HERE/"MEAN_GAP_CERTIFICATE.json").read_text())==cert
    ck("positive gap and explicit corrector-tail certificate",F(cert["positive_linear_mean_gap_coefficient"][0])>0)
    out["certificate"]=cert
    if args.replay:
        tmp=Path(tempfile.mkdtemp(prefix="fin-noc-intake-"))/"campaign"; shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1"); logs={}
        for script in ("verify_package.py","verify_noc.py"):
            run=subprocess.run(["python3",script],cwd=tmp,env=env,capture_output=True,text=True,timeout=60); assert run.returncode==0,(script,run.stdout,run.stderr)
            logs[script]={"exit_code":run.returncode,"stdout":run.stdout,"stderr":run.stderr}
        saved=json.loads((tmp/"INDEPENDENT_VERIFICATION.json").read_text())
        out["replay"]={"directory":str(tmp),"logs":logs,"source_checks":saved["n_checks"],"numerical_output_byte_equal":sha(tmp/"INDEPENDENT_VERIFICATION.json")==sha(SOURCE/"INDEPENDENT_VERIFICATION.json")}
    inputs=[SEED,Path(__file__),HERE/"certify_mean_gap.py",SOURCE/"MANIFEST.sha256",ROOT/"FIN_NATIVE_OBSERVATION_COST_NEXT_CHAIN_20261015.md"]
    inputs += [ROOT/rel/"MANIFEST.sha256" for rel in prior]
    inputs += [ROOT/rel for rel in json.loads((SOURCE/"PROVENANCE.json").read_text())["unavailable_as_separate_review_bytes"]]
    assert all(p.is_file() for p in inputs)
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["source_manifest_after"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"]); assert out["source_manifest_after"]==out["source_manifest"]
    out["source_flag_scope"]="producer point-tail estimate is not a rigorous infinite-tail certificate; floating floor examples are proposals, not exact integer identities or N0"
    out["checks"]=checks; out["status"]="PASS_SCOPED_NOC_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"source":out["source_manifest"],"matched_inputs":len(mapping),"stage_references":out["stage_reference_counts"],"checks":len(checks),"source_checks":out.get("replay",{}).get("source_checks"),"certificate":cert},indent=2))

if __name__=="__main__":main()
