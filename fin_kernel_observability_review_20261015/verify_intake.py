#!/usr/bin/env python3
"""Bounded KO intake checks; source writer runs only on a temporary copy."""
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

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/"fin_kernel_observability_campaign_20261015"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder,exclude=()):
    names=[]
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=line.split(None,1); assert name not in names and sha(folder/name)==h,name; names.append(name)
    actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}
    extra=actual-set(names)-{"MANIFEST.sha256"}-set(exclude); assert not extra,extra
    return {"matched":len(names),"mismatches":0,"uncovered":sorted(extra),"separately_hashed_exclusions":{x:sha(folder/x) for x in exclude}}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--replay",action="store_true"); args=parser.parse_args()
    checks=[]
    def ck(name,ok,**details):checks.append({"name":name,"passed":bool(ok),**details}); assert ok,name
    out={"date":"2026-10-15","baseline":"380e9dc0 plus preserved LI intake/plan/AGENTS worktree"}
    out["source_manifest"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"])
    provenance=json.loads((SOURCE/"PROVENANCE.json").read_text())
    aliases={name:h for h,name in [line.split(None,1) for line in (SOURCE/"SOURCE_INPUTS.sha256").read_text().splitlines()]}
    ck("source provenance agrees with its external-reference list",aliases==provenance["source_sha256"])
    archives={"UF_ARCHIVE":ROOT/"fin_critical_metastability_campaign_20261014/source_context/UF_PREDECESSOR.zip",
              "MT_ARCHIVE":ROOT/"fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip"}
    for name,p in archives.items():assert sha(p)==aliases[name],name
    out["external_inputs"]={"matched_archives":{k:str(v.relative_to(ROOT)) for k,v in archives.items()},
                            "unavailable_original_bytes":sorted(set(aliases)-set(archives)),"scope":"wrapper absence does not invalidate locally available producer/review files"}
    stages={}
    for i in range(1,6):
        text=(SOURCE/f"KO-{i:03d}/INPUTS.sha256").read_text()
        listed={name:h for h,name in [line.split(None,1) for line in text.splitlines()]}
        assert listed==aliases
        stages[f"KO-{i:03d}"]={"matched_archive_references":2,"unavailable_original_references":4}
    out["stage_references"]=stages
    previous={}
    for rel in ("fin_parity_record_sufficiency_review_20261015","fin_critical_metastability_review_20261015","fin_intrinsic_observation_review_20261015"):
        folder=ROOT/rel; result=manifest(folder)
        old=json.loads((folder/"verification.json").read_text())
        for p,h in old["input_hashes"].items():assert sha(ROOT/p)==h,p
        result["saved_inputs_unchanged"]=len(old["input_hashes"]); previous[rel]=result
    out["actual_review_dependencies"]=previous
    actual_missing=provenance["not_accessed_as_bytes"]
    ck("formerly unavailable review/seed files now actually present",all((ROOT/p).is_file() for p in actual_missing))
    ck("KO004 genuinely deferred in source contract",json.loads((SOURCE/"KO-004/RESULTS.json").read_text())["scientific_verdict"]=="DEFERRED_NOT_REQUIRED")

    spec=importlib.util.spec_from_file_location("versioned_seed",SEED); seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,_,A=seed.build_rank7(); g=12/lam[6]; sigma=(-1.)**np.arange(12)
    P3=np.array([[math.cos(math.pi*(i-j)/2)/6 for j in range(12)] for i in range(12)])
    P4=np.array([[math.cos(2*math.pi*(i-j)/3)/6 for j in range(12)] for i in range(12)])
    Ae=A+.01*(P3-P4)
    ck("versioned operator versus matched comparator trace diagonal",np.max(abs(np.diag(Ae)-np.diag(A)))<1e-14 and abs(np.trace(Ae)-np.trace(A))<1e-13)
    ck("parity eigenvalue and simple top preserved",np.linalg.norm((Ae-A)@sigma)<1e-14 and abs(np.linalg.eigvalsh(Ae)[-1]-lam[6])<1e-13)
    u=np.ones(12)/12
    for m in (-.01,0,.01):
        p=u+m*sigma/12
        ck(f"exact-form parity-line field diagnostic {m}",np.linalg.norm(Ae@p-A@p)<2e-15)
    # Exact independent multivariate hypergeometric enumeration with repeated labels.
    population=[0,0,1,2,3,3,4,5,6,7,8,9]; N=12; k=6
    n=[population.count(j) for j in range(12)]; tally={}
    for subset in itertools.combinations(range(N),k):
        h=tuple(sum(population[a]==j for a in subset) for j in range(12)); tally[h]=tally.get(h,0)+1
    denom=math.comb(N,k)
    for h,count in tally.items():assert F(count,denom)==F(math.prod(math.comb(n[j],h[j]) for j in range(12)),denom)
    ck("exact repeated-label hypergeometric support/normalization",sum(tally.values())==denom)
    c3=[F(1),F(0),F(-1),F(0)]; c4=[F(1),F(-1,2),F(-1,2)]
    features=[[c3[j%4],c4[j%3]] for j in population]
    means=[sum(row[r] for row in features)/N for r in range(2)]
    cov=[[sum((row[r]-means[r])*(row[s]-means[s]) for row in features)/N for s in range(2)] for r in range(2)]
    sampled=[[F(0) for _ in range(2)] for _ in range(2)]
    for h,count in tally.items():
        val=[sum(F(h[j])*([c3[j%4],c4[j%3]][r]) for j in range(12))/k-means[r] for r in range(2)]
        for r in range(2):
            for s in range(2):sampled[r][s]+=F(count*N,denom)*val[r]*val[s]
    factor=F(N*(N-k),k*(N-1))
    ck("exact finite population covariance with nonzero feature means",sampled==[[factor*x for x in row] for row in cov])
    f=np.sqrt(2)*np.cos(2*np.pi*np.arange(12)/12)
    ck("A-null mode is not Laplacian mode1",np.linalg.norm(A@f)<1e-13 and lam[1]>.7)
    ck("uniform feature variance and factor12",abs(np.mean(f*f)-1)<1e-14 and abs(f@f-12)<1e-13)
    from certify_event import build,endpoints
    certificate=build()
    assert json.loads((HERE/"EVENT_CERTIFICATE.json").read_text())==certificate
    ck("positive-series event certificate and strict gap",F(certificate["Bernoulli_TV_gap"][0])>0)
    ck("wrong primitive covariance normalization rejected",abs(1/(12-g*lam[3])-1/(1-lam[3]/lam[6]))>1)
    out["certificate"]=certificate
    out["source_arithmetic_qualification"]="Source lower-gap bookkeeping writes gamma4+.01 instead of gamma4+.01/Lambda6. Its minimum remains gamma5; new exact maximum-mode check pays the intended gap. Long source endpoint strings are not automatically outward public decimals."
    if args.replay:
        tmp=Path(tempfile.mkdtemp(prefix="fin-ko-intake-"))/"campaign"; shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1"); logs={}
        for script in ("verify_package.py","independent_verify.py"):
            run=subprocess.run(["python3",script],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
            assert run.returncode==0,(script,run.stderr)
            logs[script]={"exit_code":run.returncode,"stdout":run.stdout,"stderr":run.stderr}
        saved=json.loads((tmp/"INDEPENDENT_VERIFICATION.json").read_text())
        out["replay"]={"directory":str(tmp),"logs":logs,"source_assertions":saved["assertions_passed"],
                       "verification_byte_equal":sha(tmp/"INDEPENDENT_VERIFICATION.json")==sha(SOURCE/"INDEPENDENT_VERIFICATION.json")}
        consistent=[]
        for own,key in (("probability_A7","probability_A7_interval"),("probability_comparator","probability_comparator_interval"),
                        ("Bernoulli_TV_gap","event_probability_gap_interval"),("equal_prior_single_event_risk","ideal_equal_prior_risk_interval")):
            low=F(saved[key]["lower"].strip("[]").split(",")[0])
            high=F(saved[key]["upper"].strip("[]").split(",")[-1])
            consistent.append(F(certificate[own][0])<=low<=high<=F(certificate[own][1]))
        ck("all four producer event-number pairs consistent with reviewed coarse enclosures",all(consistent))
    inputs=[SEED,Path(__file__),HERE/"certify_event.py",SOURCE/"MANIFEST.sha256",ROOT/"FIN_KERNEL_OBSERVABILITY_NEXT_CHAIN_20261015.md"]
    inputs += [ROOT/rel/"MANIFEST.sha256" for rel in previous]
    inputs += [ROOT/rel for rel in actual_missing]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["source_manifest_after"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"])
    assert out["source_manifest"]==out["source_manifest_after"]
    out["checks"]=checks; out["status"]="PASS_SCOPED_KO_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"source":out["source_manifest"],"checks":len(checks),"source_assertions":out.get("replay",{}).get("source_assertions"),"event_certificate":certificate},indent=2))


if __name__=="__main__":main()
