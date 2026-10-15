#!/usr/bin/env python3
"""DKO bounded independent diagnostics and source replay; not stochastic proof."""
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
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import brentq
from scipy.special import softmax

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/"fin_dynamic_kernel_observability_campaign_20261015"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder,excluded=()):
    names=[]; absent_cache=[]; cache=[]; matched=0
    for row in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=row.split(None,1); assert name not in names,name; names.append(name)
        is_cache="__pycache__" in Path(name).parts or name.endswith(".pyc")
        if is_cache:
            cache.append(name)
            if not (folder/name).is_file():absent_cache.append(name); continue
        assert sha(folder/name)==h,name
        matched+=1
    actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}
    extra=actual-set(names)-{"MANIFEST.sha256"}-set(excluded); assert not extra,extra
    return {"matched":matched,"scientific_entries":len(names)-len(cache),"archival_cache_entries":cache,"optional_cache_missing":absent_cache,
            "mismatches":0,"uncovered":sorted(extra),"separately_bound":{x:sha(folder/x) for x in excluded}}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--replay",action="store_true"); args=parser.parse_args()
    checks=[]
    def ck(name,ok,**detail):checks.append({"name":name,"passed":bool(ok),**detail}); assert ok,name
    out={"date":"2026-10-15","baseline":"380e9dc0 plus preserved LI/KO worktree"}
    out["source_manifest"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"])
    declared={name:h for h,name in [row.split("  ",1) for row in (SOURCE/"SOURCE_INPUTS.sha256").read_text().splitlines()]}
    ck("root source input lists and provenance agree",declared==json.loads((SOURCE/"PROVENANCE.json").read_text())["input_files_sha256"] and (SOURCE/"INPUTS.sha256").read_text()==(SOURCE/"SOURCE_INPUTS.sha256").read_text())
    mt=ROOT/"fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip"
    assert sha(mt)==declared["FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip"]
    out["external_inputs"]={"matched_original_archive":[str(mt.relative_to(ROOT))],"unavailable_original_files":sorted(x for x in declared if x!="FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip")}
    prior={}
    for rel in ("fin_kernel_observability_review_20261015","fin_critical_metastability_review_20261015","fin_parity_record_sufficiency_review_20261015"):
        folder=ROOT/rel; row=manifest(folder); v=json.loads((folder/"verification.json").read_text())
        for path,h in v["input_hashes"].items():assert sha(ROOT/path)==h,path
        row["saved_inputs_unchanged"]=len(v["input_hashes"]); prior[rel]=row
    out["actual_prior_reviews"]=prior
    stages={}
    for i in range(1,6):
        matched=[]; missing=[]
        for row in (SOURCE/f"DKO-{i:03d}/INPUTS.sha256").read_text().splitlines():
            h,name=row.split("  ",1)
            if name.startswith("source_KO_archive::"):
                rel=name.split("::",1)[1]; p=ROOT/"fin_kernel_observability_campaign_20261015"/rel
                assert sha(p)==h,name; matched.append(name)
            else:missing.append({"name":name,"reported_sha256":h,"status":"original upload/wrapper not available"})
        stages[f"DKO-{i:03d}"]={"matched_extracted_files":len(matched),"unavailable_original_references":missing}
    out["stage_references"]=stages

    spec=importlib.util.spec_from_file_location("versioned_seed",SEED); seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,_,A=seed.build_rank7(); g=12/lam[6]; rs=(1,3,4)
    features=np.sqrt(2)*np.cos(2*np.pi*np.arange(12)[:,None]*np.array(rs)/12)
    ck("uniform Fourier normalization",np.max(abs(features.T@features/12-np.eye(3)))<2e-14)
    def direct(labels,k):
        N=len(labels); n=np.bincount(labels,minlength=12); alpha=k/N
        c=np.r_[np.full(k,N/k-1),np.full(N-k,-1.)]
        drift=np.zeros(6); bracket=np.zeros((6,6)); private_error=np.zeros(3)
        full_mean=softmax(g/N*(A@n))@features
        for a,old in enumerate(labels):
            deleted=n.copy(); deleted[old]-=1; p=softmax(g/N*(A@deleted))
            private_error+=c[a]*(p@features-full_mean)/math.sqrt(N)
            for j in range(12):
                jump=(features[j]-features[old])/math.sqrt(N)
                alljump=np.r_[jump,c[a]*jump]
                drift+=p[j]*alljump; bracket+=p[j]*np.outer(alljump,alljump)
        FF=features[labels].sum(axis=0)/math.sqrt(N)
        UU=(c[:,None]*features[labels]).sum(axis=0)/math.sqrt(N)
        return drift,bracket,FF,UU,private_error
    fixtures=[]
    rng=np.random.default_rng(20261015)
    for N,k in ((24,8),(120,60),(120,12)):
        labels=rng.integers(0,12,size=N)
        drift,bracket,FF,UU,error=direct(labels,k)
        ck(f"exact private drift identity actual finite jumps N{N} k{k}",np.max(abs(drift[3:]+UU-error))<2e-13)
        ck(f"actual joint bracket positivity N{N} k{k}",np.linalg.eigvalsh(bracket).min()>-2e-13)
        coeff=np.r_[np.full(k,N/k-1),np.full(N-k,-1.)]
        ck(f"uniform absolute coefficient cancellation N{N} k{k}",abs(abs(coeff).sum()-2*(N-k))<1e-12)
        fixtures.append({"N":N,"k":k,"FF":bracket[:3,:3].tolist(),"UU":bracket[3:,3:].tolist(),"FU":bracket[:3,3:].tolist(),"private_residual":error.tolist(),"scope":"algebra only, not a critical-sector trajectory"})
    out["actual_generator_fixtures"]=fixtures
    # Uniform g0 r3 squared jumps depend on the old label:1+2cos3(old)^2.
    c3=[F(1),F(0),F(-1),F(0)]; h=[1+2*c3[j%4]**2 for j in range(12)]
    vals=[]
    for subset in itertools.combinations(range(12),4):
        B=set(subset); c=[F(2) if j in B else F(-1) for j in range(12)]
        vals.append(sum(c[j]*h[j] for j in range(12))/12)
    ck("global-private cross bracket is zero on conditional average, not every block",sum(vals)==0 and max(abs(v) for v in vals)>0)
    ck("general-feature parity shortcut is rejected",set(h)=={F(1),F(3)})
    gamma=1-lam[3]/lam[6]; control=1-(lam[3]+.01)/lam[6]
    inv=[]
    for nu in (.7,1.9):
        q=math.exp(-nu*1.3)
        covariance=math.exp(-gamma*nu*1.3)/gamma+math.exp(-nu*1.3)
        rho=covariance/(1/gamma+1)
        s=2/math.pi*math.asin(rho); r=math.sin(math.pi*s/2)
        recovered=brentq(lambda z:(q**z+z*q)/(1+z)-r,1e-8,.99999999)
        inv.append(recovered)
    ck("common clock removed by null-channel inversion",max(abs(x-gamma) for x in inv)<2e-12)
    pairs=[(2,-3),(-5,-7),(0,3)]
    sign=lambda x:(x>0)-(x<0)
    ck("actual sign-products invariant under fixed nonzero gains",all(sign(a*x)*sign(a*y)==sign(x)*sign(y) for x,y in pairs for a in (-3,2,7)))
    ck("time-changing gain sign is not licensed",sign(2*2)*sign(-2*3)!=sign(2)*sign(3))
    # Independent Cholesky trace/logdet route, not source relative-eigenvalue KL.
    def K(rate,alpha,M):
        d=abs(np.arange(M)[:,None]-np.arange(M)[None,:])
        return alpha/rate*np.exp(-rate*d)+(1-alpha)*np.exp(-d)
    def KL(a,b):
        cb=cho_factor(b,lower=True); ca=cho_factor(a,lower=True)
        return .5*(np.trace(cho_solve(cb,a))-len(a)+2*np.log(np.diag(cb[0])).sum()-2*np.log(np.diag(ca[0])).sum())
    saved=json.loads((SOURCE/"INDEPENDENT_VERIFICATION.json").read_text())
    errors=[]
    for row in saved["gaussian_resource_fixtures"]:
        value=KL(K(gamma,row["alpha"],row["M"]),K(control,row["alpha"],row["M"]))
        errors.append(abs(value-row["KL0_1"]))
    ck("dependent Gaussian KL fixtures independently reproduced",max(errors)<2e-12,max_difference=max(errors))
    from certify_sign_event import build
    certificate=build(); assert json.loads((HERE/"SIGN_EVENT_CERTIFICATE.json").read_text())==certificate
    ck("validated ideal sign-event separation",F(certificate["Bernoulli_TV_gap"][0])>0)
    for name,key in (("probability_A7","p_A7"),("probability_control","p_control"),("Bernoulli_TV_gap","p_difference"),("equal_prior_single_event_risk","equal_prior_risk")):
        ck("source nominal "+key+" consistent with independent enclosure",float(certificate[name][0])<=saved["snapshot_event"][key]<=float(certificate[name][1]))
    out["certificate"]=certificate
    out["source_flag_qualification"]="The amplitude-units source flag repeats its defining formula;32 flags are bounded checks, not32 independent proofs. Intake gain tests are separate exact pointwise checks."
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix="fin-dko-intake-")); archive=temp/"archive"; tmp=temp/"campaign"
        shutil.copytree(SOURCE,archive)
        shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1"); logs={}
        for script in ("verify_package.py","independent_verify.py"):
            target=archive if script=="verify_package.py" else tmp
            run=subprocess.run(["python3",script],cwd=target,env=env,capture_output=True,text=True,timeout=60)
            # The producer lists one historical cache. Its optional clean-checkout
            # absence is packaging-only, never a scientific checksum waiver.
            if script=="verify_package.py" and run.returncode!=0 and out["source_manifest"]["optional_cache_missing"]:
                report=json.loads(run.stdout)
                assert report["errors"] and all("__pycache__" in error or ".pyc" in error for error in report["errors"]),report
            else:assert run.returncode==0,(script,run.stdout,run.stderr)
            logs[script]={"exit_code":run.returncode,"stdout":run.stdout,"stderr":run.stderr}
        result=json.loads((tmp/"INDEPENDENT_VERIFICATION.json").read_text())
        out["replay"]={"directory":str(tmp),"hash_only_archive_copy":str(archive),"historical_cache_executed":False,
                       "logs":logs,"source_checks":result["n_checks"],"numerical_output_byte_equal":sha(tmp/"INDEPENDENT_VERIFICATION.json")==sha(SOURCE/"INDEPENDENT_VERIFICATION.json")}
    inputs=[SEED,Path(__file__),HERE/"certify_sign_event.py",SOURCE/"MANIFEST.sha256",ROOT/"FIN_DYNAMIC_KERNEL_OBSERVABILITY_NEXT_CHAIN_20261015.md"]
    inputs += [ROOT/rel/"MANIFEST.sha256" for rel in prior]
    inputs += [ROOT/"fin_kernel_observability_review_20261015/PROOF_COMPLETIONS_AND_ARITHMETIC.md",ROOT/"fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md"]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["source_manifest_after"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"]); assert out["source_manifest_after"]==out["source_manifest"]
    out["checks"]=checks; out["status"]="PASS_SCOPED_DKO_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"source":out["source_manifest"],"checks":len(checks),"source_checks":out.get("replay",{}).get("source_checks"),"certificate":certificate},indent=2))


if __name__=="__main__":main()
