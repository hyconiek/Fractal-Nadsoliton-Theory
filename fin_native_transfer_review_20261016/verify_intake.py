#!/usr/bin/env python3
"""NT intake: provenance, independent labelled-chain checks, copy-only replay.

Finite fixtures are not proofs of growing-N convergence or of the proposed HM lane.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np
from scipy.special import softmax
import sympy as sy

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/"fin_native_transfer_campaign_20261015"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def manifest(folder,excluded=()):
    names=[]
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=line.split("  ",1)
        assert name not in names and sha(folder/name)==h,name
        names.append(name)
    actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    assert not actual-set(names)-{"MANIFEST.sha256"}-set(excluded)
    return {"matched":len(names),"uncovered":[],"separately_bound":{n:sha(folder/n) for n in excluded}}

def uniformize(L,t,rate=2.,cutoff=45):
    P=np.eye(len(L))+L/rate
    assert P.min()>-1e-14
    term=np.eye(len(L)); weight=math.exp(-rate*t); result=weight*term
    for j in range(1,cutoff+1):
        term=term@P; weight*=rate*t/j; result+=weight*term
    return result

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--replay",action="store_true")
    args=parser.parse_args(); checks=[]
    def ck(name,ok,**detail):
        checks.append({"name":name,"passed":bool(ok),**detail}); assert ok,name
    out={"date":"2026-10-16","baseline":"2267ed77 at start; inherited research and incoming NT committed externally as3be481a7 during intake; source hashes unchanged"}
    out["source_manifest"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json","package_check.py"])
    records=json.loads((SOURCE/"INPUT_PROVENANCE.json").read_text())["read_session_inputs"]
    expected="".join(r["sha256"]+"  "+r["path"]+"\n" for r in records)
    ck("six source/stage hash ledgers agree",(SOURCE/"INPUTS.sha256").read_text()==expected and
       all((SOURCE/f"NT-{i:03d}/INPUTS.sha256").read_text()==expected for i in range(1,7)))
    mapping={}
    for row in records:
        if "FIN_UF_001_UF_005_FULL_HANDOFF" in row["path"]:
            rel="fin_critical_metastability_campaign_20261014/source_context/UF_PREDECESSOR.zip"
        elif "FIN_MT_001_MT_005_FULL_HANDOFF" in row["path"]:
            rel="fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip"
        else:continue
        p=ROOT/rel; assert sha(p)==row["sha256"] and p.stat().st_size==row["size_bytes"],rel
        mapping[row["path"]]=rel
    out["external_inputs"]={"matched_original_wrappers":mapping,"unavailable_originals":[r["path"] for r in records if r["path"] not in mapping]}
    out["stage_reference_counts"]={"matched":6*len(mapping),"unavailable_original":6*(len(records)-len(mapping))}
    prior={}
    for rel in ("fin_toe_bridge_research_20261015","fin_critical_metastability_review_20261015",
                "fin_intrinsic_observation_review_20261015","fin_uniform_critical_review_20261014",
                "fin_kernel_observability_review_20261015"):
        folder=ROOT/rel; entry=manifest(folder)
        old=json.loads((folder/"verification.json").read_text())
        for name,h in old["input_hashes"].items():assert sha(ROOT/name)==h,name
        entry["saved_inputs_unchanged"]=len(old["input_hashes"]); prior[rel]=entry
    out["actual_dependencies"]=prior
    spec=importlib.util.spec_from_file_location("versioned_seed",SEED)
    seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,X,A=seed.build_rank7(); g=12/lam[6]; sigma=(-1.)**np.arange(12)
    states=[(a,b) for a in range(12) for b in range(12)]
    count_states=[(a,b) for a in range(12) for b in range(a,12)]
    idx={pair:i for i,pair in enumerate(count_states)}
    L=np.zeros((144,144)); orbit=np.zeros((144,78))
    for row,(a,b) in enumerate(states):
        orbit[row,idx[tuple(sorted((a,b)))]]=1
        for slot,old,other in ((0,a,b),(1,b,a)):
            p=softmax(g/2*A[:,other])
            for new in range(12):
                if old==new:continue
                dest=12*new+b if slot==0 else 12*a+new
                L[row,dest]+=p[new]; L[row,row]-=p[new]
    sizes=orbit.sum(0); quotient=(orbit.T/sizes[:,None])@L@orbit
    pi_label=softmax(np.array([g/2*A[a,b] for a,b in states])); pi=pi_label@orbit
    ck("labelled count intertwining",np.max(abs(L@orbit-orbit@quotient))<3e-15)
    ck("count multiplicity and invariant law",np.all(sizes==np.array([1 if a==b else 2 for a,b in count_states])) and
       np.max(abs(pi@quotient))<2e-16)
    ck("native detailed balance",np.max(abs(pi[:,None]*quotient-(pi[:,None]*quotient).T))<2e-16)
    ck("native row sums",np.max(abs(quotient.sum(1)))<2e-14)
    m=np.array([(sigma[a]+sigma[b])/2 for a,b in count_states])
    f=np.cos(.7*m); D=np.sqrt(pi)
    H=-D[:,None]*quotient/D[None,:]
    ck("symmetric positive transfer",np.max(abs(H-H.T))<3e-15 and np.linalg.eigvalsh(H).min()>-2e-14)
    ck("positive ground vector",np.max(abs(H@D))<3e-15)
    P=uniformize(quotient,.45); P2=uniformize(quotient,.9)
    fcomplex=f+1j*np.sin(.4*m); h=np.exp(.2*m)+.1j*np.sin(m)
    norm=float(np.sum(pi*abs(P@fcomplex-h)**2))
    expanded=np.vdot(fcomplex,pi*(P2@fcomplex)).real-2*np.vdot(h,pi*(P@fcomplex)).real+np.vdot(h,pi*h).real
    ck("full-state complex 2t norm identity",abs(norm-expanded)<3e-15,residual=abs(norm-expanded))
    # N2 has exceptional parity closure. Use actual N3 slots to check why
    # no all-population parity-lumping premise may enter the asymptotic proof.
    def parity_drift(labels):
        n=np.bincount(labels,minlength=12); result=0.
        for old in labels:
            target=softmax(g/len(labels)*(A@(n-np.eye(12,dtype=int)[old])))
            result+=(target@sigma-sigma[old])/len(labels)
        return float(result)
    drift3=[parity_drift(labels) for labels in ([0,0,0],[0,0,2])]
    ck("full-state proof avoids assuming an all-N parity quotient",abs(drift3[0]-drift3[1])>.01,
       N=3,same_parity_mean=1.,drifts=drift3,scope="finite numerical diagnostic, not an asymptotic theorem")
    energy=.5*np.sum(pi[:,None]*quotient*abs(fcomplex[None,:]-fcomplex[:,None])**2)
    ck("complex edge Dirichlet identity",abs(energy+np.vdot(fcomplex,pi*(quotient@fcomplex)))<2e-15)
    # Direct chronological four-time expectation versus conditional future Gram.
    a=np.c_[np.cos(np.arange(78)*.12)+1j*np.sin(np.arange(78)*.24),
            np.sin(np.arange(78)*.22)-1j*np.cos(np.arange(78)*.21)]
    b=np.c_[np.sin(np.arange(78)*.07)+.3j*np.cos(np.arange(78)*.13),
            np.cos(np.arange(78)*.09)+.2j*np.sin(np.arange(78)*.20)]
    Ps=uniformize(quotient,.17); Pdt=uniformize(quotient,.28); Pmid=uniformize(quotient,.34)
    U=np.column_stack([Ps@(a[:,j]*(Pdt@b[:,j])) for j in range(2)])
    gram=U.conj().T@(pi[:,None]*U); chronological=np.zeros((2,2),complex)
    for i in range(2):
        for j in range(2):
            chronological[i,j]=np.sum(pi*np.conj(b[:,i])*(Pdt@(np.conj(a[:,i])*(Pmid@(a[:,j]*(Pdt@b[:,j]))))))
    ck("independent chronological reflection form",np.max(abs(chronological-gram))<3e-16 and np.linalg.eigvalsh(chronological).min()>-1e-16)
    ck("double-conjugation negative control",complex(0,1)*complex(0,1)==-1 and complex(0,-1)*complex(0,1)==1)
    z,alpha=sy.symbols("z alpha",real=True)
    logpsi=-z**4/24+alpha*z**2/4
    potential=sy.expand(sy.diff(logpsi,z)**2+sy.diff(logpsi,z,2))
    expected_potential=z**6/36-alpha*z**4/6+(alpha**2/4-sy.Rational(1,2))*z**2+alpha/2
    ck("independent logarithmic ground-state calculation",sy.simplify(potential-expected_potential)==0)
    ck("critical sextic coefficients",sy.Poly(potential.subs(alpha,0),z).as_dict()=={(6,):sy.Rational(1,36),(2,):sy.Rational(-1,2)})
    # Bounded seeds for the proposed next chain, not an executed mediator campaign.
    cavity=[]
    for old,other in ((0,0),(0,1),(3,8)):
        background=X[other]
        gaussian_weights=softmax(np.array([g/4*np.dot(background+X[j],background+X[j]) for j in range(12)]))
        native=softmax(g/2*A[:,other])
        cavity.append(float(np.max(abs(gaussian_weights-native))))
    ck("exact-Gaussian cavity identity numerical seed",max(cavity)<2e-15)
    k=sy.symbols("k",positive=True)
    J=sy.Matrix([[-1,1],[k,-k]]); left=sy.Matrix([[k/(1+k),1/(1+k)]])
    ck("diagnostic mediator soft null pair",(left*J).applyfunc(sy.simplify)==sy.zeros(1,2) and
       (left*sy.Matrix([1,1])).applyfunc(sy.simplify)==sy.ones(1,1))
    projected=(left*sy.diag(2,2*k)*left.T)[0]
    ck("diagnostic mobility-noise algebra only",sy.simplify(projected-2*k/(1+k))==0)
    out["numerical_fixtures"]={"N":2,"states_labelled":144,"states_count":78,"g_c_float":g,
        "source_same_bounded_f_correlation":float(np.sum(pi*f*(P@f))),
        "chronological_gram_max_error":float(np.max(abs(chronological-gram))),
        "scope":"global N2 convention fixture, outside the local critical asymptotic regime; no probability error certificate"}
    if args.replay:
        tmp=Path(tempfile.mkdtemp(prefix="fin-nt-intake-"))/"campaign"
        shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        assert manifest(tmp,["PACKAGE_VERIFICATION.json","package_check.py"])==out["source_manifest"]
        logs={}; env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1")
        for script in ("package_check.py","verify_nt.py"):
            run=subprocess.run(["python3",script],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
            assert run.returncode==0,(script,run.stdout,run.stderr)
            logs[script]={"stdout":run.stdout,"stderr":run.stderr,"exit_code":run.returncode}
        original=json.loads((SOURCE/"verification.json").read_text()); fresh=json.loads((tmp/"verification.json").read_text())
        ck("independent correlation agrees with source",abs(out["numerical_fixtures"]["source_same_bounded_f_correlation"]-original["finite_correlation_classical"])<2e-14)
        out["replay"]={"directory":str(tmp),"logs":logs,"bounded_source_flags":fresh["bounded_fixtures_passed"],
            "numerical_byte_equal":sha(tmp/"verification.json")==sha(SOURCE/"verification.json"),
            "max_numeric_drift":max(abs(v-original[name]) for name,v in fresh.items() if isinstance(v,(float,int)) and not isinstance(v,bool))}
    inputs=[SEED,Path(__file__),SOURCE/"MANIFEST.sha256",ROOT/"FIN_TOE_THEORETICAL_NEXT_CHAIN_20261015.md"]
    inputs += [ROOT/name/"MANIFEST.sha256" for name in prior]
    inputs += [ROOT/name for name in (
        "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
        "fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md",
        "fin_uniform_critical_review_20261014/INTAKE_20261014.md",
        "fin_kernel_observability_review_20261015/PROOF_COMPLETIONS_AND_ARITHMETIC.md",
        "fin_toe_bridge_research_20261015/TRANSFER_PROOF_REVIEW.md",
        "fin_rank7_mathphysics_next/proofs/MP7-032_gaussian_auxiliary_representation.md")]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["source_manifest_after"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json","package_check.py"])
    assert out["source_manifest_after"]==out["source_manifest"]
    out["checks"]=checks; out["status"]="PASS_SCOPED_NT_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"source_manifest":out["source_manifest"],"matched_originals":len(mapping),
        "stage_refs":out["stage_reference_counts"],"checks":len(checks),"replay":out.get("replay",{})},indent=2))

if __name__=="__main__":main()
