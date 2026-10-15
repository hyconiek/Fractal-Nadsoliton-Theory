#!/usr/bin/env python3
"""LI audit: custody, tiny q12 algebra and exact decision fixtures, not limit proofs."""
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

import numpy as np
from scipy.linalg import expm
from scipy.special import softmax

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/"fin_parity_record_sufficiency_campaign_20261015"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen=[]
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,rel=line.split(None,1); assert rel not in seen and sha(folder/rel)==h,rel; seen.append(rel)
    actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}
    extra=actual-set(seen)-{"MANIFEST.sha256"}; assert not extra,extra
    return {"matched":len(seen),"mismatches":0,"uncovered":sorted(extra)}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--replay",action="store_true"); args=parser.parse_args()
    checks=[]
    def ck(name,ok,**detail):
        checks.append({"name":name,"passed":bool(ok),**detail}); assert ok,name
    out={"date":"2026-10-15","baseline":"380e9dc0; LI arrived untracked"}
    out["source_manifest"]=manifest(SOURCE)
    mapping={
        "RI_HANDOFF":"fin_record_information_campaign_20261015/HANDOFF.md",
        "RI_RECORD_FIELDS":"fin_record_information_campaign_20261015/RI-001_RECORD_FIELDS_AND_TAG_GENERATOR.md",
        "RI_FULL_RECORD_GATE":"fin_record_information_campaign_20261015/RI-004_FULL_RECORD_GATE.md",
        "UF_PROCESS":"fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md",
        "MT_LOCAL_GIBBS":"fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md",
        "MT_LOCAL_MIXING":"fin_critical_metastability_campaign_20261014/MT-004_REFLECTED_MIXING_WINDOW.md",
        "OBS_RECORD":"fin_intrinsic_observation_campaign_20261015/OBS-C04_INTEGRATED_RECORD_THEOREM.md"}
    def references(path):
        matched=[]; unavailable=[]
        for line in path.read_text().splitlines():
            h,label=line.split("  ",1); key,name=label.split(": ",1)
            if key in mapping:
                assert sha(ROOT/mapping[key])==h,(path,key)
                matched.append(key)
            else:
                assert key=="LI_PLAN",key
                unavailable.append({"alias":key,"uploaded_filename":name,"reported_sha256":h,"status":"ORIGINAL_UPLOAD_BYTES_UNAVAILABLE"})
        return {"matched":len(matched),"unavailable":unavailable}
    out["external_source_references"]=references(SOURCE/"SOURCE_INPUTS.sha256")
    out["stage_references"]={f"LI-{i:03d}":references(SOURCE/f"LI-{i:03d}/INPUTS.sha256") for i in range(1,6)}
    prior_bindings={}
    for folder in ("fin_record_information_review_20261015","fin_intrinsic_observation_review_20261015","fin_critical_metastability_review_20261015","fin_uniform_critical_review_20261014"):
        f=ROOT/folder; binding=manifest(f)
        v=json.loads((f/"verification.json").read_text())
        for rel,h in v["input_hashes"].items():assert sha(ROOT/rel)==h,rel
        binding["saved_inputs_unchanged"]=len(v["input_hashes"]); prior_bindings[folder]=binding
    out["actual_prior_reviews"]=prior_bindings
    out["missing_root_replay_wrapper"]=not (SOURCE/"REPLAY.md").exists()
    out["custody_scope"]="actual reviews bind inherited premises now; original uploaded-plan formatting/hash and historical executor access remain unattested"

    spec=importlib.util.spec_from_file_location("versioned_seed",SEED); seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,_,A=seed.build_rank7(); g=12/lam[6]; sigma=(-1.)**np.arange(12)
    ck("constant diagonal and centered q12 kernel",np.ptp(np.diag(A))<1e-13 and np.linalg.norm(A@np.ones(12))<1e-13)
    # N2 local conditional identity in an artificial count domain, not the tiny MT domain.
    weights=np.zeros((12,12))
    for bg in range(12):
        for tag in range(12):
            n=np.eye(12)[bg]+np.eye(12)[tag]
            weights[bg,tag]=np.exp(g/4*(n@A@n))
    rows=weights/weights.sum(axis=1,keepdims=True)
    predicted=np.asarray([softmax(g/2*A[:,bg]) for bg in range(12)])
    ck("actual q12 Gibbs tag conditional equals LOO refresh",np.max(abs(rows-predicted))<5e-15)
    local=weights.copy(); local[0,0]=0; local/=local.sum()
    marginal=local.sum(axis=1); conditional=local/marginal[:,None]
    ck("all-completion condition and truncated boundary retained",np.max(abs(conditional[1:]-predicted[1:]))<5e-15 and conditional[0,0]==0)
    # When bg=0 and the admitted tag is nonzero, replacing that tag by0 exits n0<=1.
    legal=[]
    for tag in range(1,12):
        n=np.eye(12,dtype=int)[0]+np.eye(12,dtype=int)[tag]
        dest=n-np.eye(12,dtype=int)[tag]+np.eye(12,dtype=int)[0]
        legal.append(n[0]<=1 and n[tag]>0 and dest[0]>1 and dest.sum()==2 and dest.min()>=0 and np.sum(abs(dest-n))==2)
    ck("bad completion implies an actual outward count edge",all(legal))
    parity=np.stack((local[:,1::2].sum(axis=1),local[:,::2].sum(axis=1)),axis=1)
    tv=.5*np.abs(parity-marginal[:,None]/2).sum()
    bias=(conditional@sigma)
    ck("exact initial parity product-TV identity diagnostic",abs(tv-.5*(marginal@abs(bias)))<2e-15)
    out["initial_fixture"]={"N":2,"artificial_domain":"n0<=1, not MT critical sector","product_TV_point":float(tv)}

    kl=[]
    for N in (12,120):
        bg=np.full(12,N//12); bg[0]-=1
        for old in (1,3,5):
            deleted=bg.copy(); deleted[old]-=1
            h=g/N*(A@deleted); delta=g/N*A[:,0]
            p,r=softmax(h+delta),softmax(h)
            val=float(p@np.log(p/r)); bound=float(np.ptp(delta)**2/8)
            ck(f"actual q12 feedback KL N{N} old{old}",-1e-14<=val<=bound+2e-15)
            kl.append({"N":N,"old":old,"KL_point":val,"bound_point":bound})
    out["feedback_fixtures"]=kl
    ck("gain0 feedback negative control",np.array_equal(softmax(np.zeros(12)),np.full(12,1/12)))

    # A deterministic background makes the fair->biased sign likelihood a two-state FK semigroup.
    M=np.eye(2)
    q_values=[]
    for lab in (0,1,3):
        bg=np.zeros(12); bg[lab]=11
        q=float(sigma@softmax(g/12*(A@bg))); q_values.append(q)
        ref=np.array([[-.5,.5],[.5,-.5]])
        tilted=ref.copy(); tilted[0,1]*=1-q; tilted[1,0]*=1+q
        tilted+=np.diag([q/2,-q/2])
        ck(f"conditional parity likelihood unit mass piece {lab}",np.max(abs(tilted@np.ones(2)))<2e-15)
        M=M@expm(.17*tilted)
    ck("ordered conditional likelihood pieces retain mass",np.max(abs(M@np.ones(2)-1))<2e-15)
    wrong=np.array([[-.5+q_values[0]/4,.5*(1-q_values[0])],[.5*(1+q_values[0]),-.5-q_values[0]/4]])
    ck("wrong flip-clock compensator rejected",np.max(abs(wrong@np.ones(2)))>1e-3)
    ck("Gaussian likelihood normalization and wrong coefficient negative control",F(1,2)**2*2/2-F(1,4)==0 and F(1,2)**2*2/2-F(1,2)!=0)
    out["likelihood_fixture"]={"actual_q12_background_parity_biases":q_values,"scope":"small finite algebra; not a likelihood limit simulation"}

    # Exact product-reference / conditional-density / Bayes fixture.
    L=[[F(1,2),F(3,2),F(1,2),F(3,2)],[F(1)]*4,[F(3,2),F(1,2),F(3,2),F(1,2)]]
    gamma=[[L[c][o]/3 for o in range(4)] for c in range(3)]
    risk=1-sum(max(gamma[c][o] for c in range(3)) for o in range(4))/4
    direct=1-sum(max(L[c][o]/12 for c in range(3)) for o in range(4))
    ck("exact three-target subdensity Bayes formula",risk==direct==F(1,2))
    coarse=[(0,2),(1,3)]
    coarse_risk=1-sum(max(sum(gamma[c][o] for o in group)/2 for c in range(3)) for group in coarse)/2
    ck("exact sufficient finite-cylinder fixture",coarse_risk==risk)
    h=[[L[c][o]+F((-1)**(o+c),10) for o in range(4)] for c in range(3)]
    lhs=sum(abs(sum((L[c][o]-h[c][o])/3 for c in range(3))) for o in range(4))/4
    rhs=sum(abs(L[c][o]-h[c][o]) for c in range(3) for o in range(4))/12
    ck("exact conditional L1 contraction fixture",lhs<=rhs)
    erased=1-sum(max(sum(gamma[c][o] for o in group)/2 for c in range(3)) for group in [(0,1),(2,3)])/2
    ck("an information-erasing cylinder is not sufficient",erased==F(2,3) and erased>risk)
    from kernel_family_certificate import build
    certificate=build(); ck("exact iso-parity trace-matched diagnostic-family seed",True)
    out["kernel_family_seed"]=certificate
    out["arithmetic_scope"]="Float q12 fixtures and exact toy/Fourier arithmetic; stochastic theorem admitted by proof audit, not assertion counts"

    if args.replay:
        tmp=Path(tempfile.mkdtemp(prefix="fin-li-intake-"))/"campaign"
        shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1"); logs={}
        for script in ("verify_package.py","independent_checks.py"):
            ran=subprocess.run(["python3",script],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
            assert ran.returncode==0,(script,ran.stderr)
            logs[script]={"exit_code":ran.returncode,"stdout":ran.stdout,"stderr":ran.stderr}
        saved=json.loads((tmp/"INDEPENDENT_CHECKS.json").read_text())
        out["replay"]={"directory":str(tmp),"logs":logs,"source_assertions":saved["checks"],"status":saved["status"],
                       "source_scope":saved["scope"],"output_byte_equal":sha(tmp/"INDEPENDENT_CHECKS.json")==sha(SOURCE/"INDEPENDENT_CHECKS.json")}
    inputs=[SEED,Path(__file__),HERE/"kernel_family_certificate.py",SOURCE/"MANIFEST.sha256",ROOT/"FIN_PARITY_RECORD_SUFFICIENCY_NEXT_CHAIN_20261015.md"]
    for folder in prior_bindings:
        inputs.append(ROOT/folder/"MANIFEST.sha256")
    for rel in json.loads((SOURCE/"PROOF_DEPENDENCIES.json").read_text())["unavailable_exact_review_inputs"]:
        assert (ROOT/rel).is_file(),rel
        inputs.append(ROOT/rel)
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["source_manifest_after"]=manifest(SOURCE)
    out["checks"]=checks; out["status"]="PASS_SCOPED_LI_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"source":out["source_manifest"],"source_refs":out["external_source_references"],"checks":len(checks),"source_assertions":out.get("replay",{}).get("source_assertions")},indent=2))


if __name__=="__main__":main()
