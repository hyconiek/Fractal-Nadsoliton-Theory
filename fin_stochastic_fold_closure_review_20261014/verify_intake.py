#!/usr/bin/env python3
"""Small intake checks. Scripts do not prove weak convergence or energy bounds."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import html
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
import numpy as np
from scipy.optimize import root
from scipy.special import softmax

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"fin_stochastic_fold_closure_campaign_20261014"
THREVIEW=ROOT/"fin_intrinsic_theory_review_20261014"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
FOLD=ROOT/"fin_rank7_followup/certificates/R7P-031_simple_fold.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def manifest(folder, name="MANIFEST.sha256", coverage=False):
    entries=[]
    for line in (folder/name).read_text().splitlines():
        expected, rel=line.split(None,1)
        path=(folder/rel.lstrip("*")).resolve()
        assert sha(path)==expected, path
        entries.append(str(path.relative_to(folder.resolve())) if path.is_relative_to(folder.resolve()) else rel)
    out={"matched":len(entries),"mismatches":0}
    if coverage:
        actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
        extra=actual-set(entries)-{name}
        assert not extra, extra
        out["uncovered"]=sorted(extra)
    return out

def normalized(text):
    return " ".join(re.sub(r"\\([\\`*_{}\[\]()#+.!><:|\-])",r"\1",html.unescape(text)).replace("**","").split())

def exact_exit(cert):
    assert cert["M"]==25 and cert["c_exponent_rational"]=="2/3"
    c=F(2,3); m=cert["M"]
    ip=sum(c**k/F(math.factorial(k)*(3*k+1)) for k in range(m+1))
    im=sum((-c)**k/F(math.factorial(k)*(3*k+1)) for k in range(m+1))
    e=3*c**(m+1)/math.factorial(m+1)
    assert F(cert["I_negative_exact_rational_truncation"])==im
    assert F(cert["I_positive_exact_rational_truncation"])==ip
    assert F(cert["Taylor_remainder_each_integral_strict_upper"])==e
    low=(im-e)/(im-e+ip+e); high=(im+e)/(im+e+ip-e)
    lo,hi=map(F,cert["P_exit_upper_public_outward_interval"])
    ll,lh=map(F,cert["P_exit_lower_public_outward_interval"])
    assert 0<lo<=low<high<=hi<1
    assert ll<=1-high<=1-low<=lh
    return {"upper":cert["P_exit_upper_public_outward_interval"],
            "lower":cert["P_exit_lower_public_outward_interval"],"backend":"stdlib exact Fraction comparison",
            "event":"eventual exit side of canonical diffusion, NOT exit by time 1"}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--replay",action="store_true")
    args=parser.parse_args()
    checks=[]
    def ck(name, condition, **values):
        checks.append({"name":name,"passed":bool(condition),**values})
        assert condition,name
    out={"date":"2026-10-14","baseline":"cf1bf9ae"}
    out["source_manifest"]=manifest(SOURCE,coverage=True)
    out["source_inputs"]=manifest(SOURCE,"SOURCE_INPUTS.sha256")
    out["task_inputs"]={f"SF-{i:03d}":manifest(SOURCE/f"SF-{i:03d}","INPUTS.sha256") for i in range(1,6)}
    out["TH_review_manifest"]=manifest(THREVIEW,coverage=True)
    old=json.loads((THREVIEW/"verification.json").read_text())
    for rel,h in old["input_hashes"].items():
        assert sha(ROOT/rel)==h,rel
    out["inherited_verification_input_hashes_unchanged"]=len(old["input_hashes"])
    plan=ROOT/"FIN_STOCHASTIC_FOLD_CLOSURE_NEXT_CHAIN_20261014.md"
    attached=SOURCE/"source_context/SF_PLAN_USER.md"
    ck("attached plan normalized equality",normalized(plan.read_text())==normalized(attached.read_text()))
    current_intake=THREVIEW/"INTAKE_20261014.md"
    user_intake=SOURCE/"source_context/TH_INTAKE_USER.md"
    ck("attached TH intake normalized equality",normalized(current_intake.read_text())==normalized(user_intake.read_text()))
    out["attached_texts"]={"plan_byte_identical":sha(plan)==sha(attached),
                           "TH_intake_byte_identical":sha(current_intake)==sha(user_intake)}
    ck("actual fold certificate hash",sha(FOLD)=="5537174281996415a20ea3dffce23c33258c51396b3850aedf4168d9280e352e")
    fold=json.loads(FOLD.read_text())
    ck("accepted fold inertia",fold["H7_inertia_at_root"]==[0,1,6] and fold["strict_inclusion"])
    archive=SOURCE/"source_context/FIN_TH_001_TH_005_FULL_HANDOFF_20261014.zip"
    unmatched=[];matched=0
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            if name.endswith("/"):
                continue
            p=ROOT/name
            if not p.exists() or hashlib.sha256(z.read(name)).hexdigest()!=sha(p):
                unmatched.append(name)
            else:
                matched+=1
    ck("embedded TH archive members",not unmatched,matched=matched,unmatched=unmatched)
    out["embedded_TH_archive"]={"matched":matched,"unmatched":unmatched}

    spec=importlib.util.spec_from_file_location("intake_safe_seed",SEED)
    seed=importlib.util.module_from_spec(spec);spec.loader.exec_module(seed)
    _,_,lam,X,A=seed.build_rank7()
    supplied=json.loads((SOURCE/"SF-001/NUMERICAL_CHART.json").read_text())
    C=X[:,[0,2,4,6]]
    def eq(z):
        pp=softmax(C@z[:4]); Sigma=np.diag(pp)-np.outer(pp,pp)
        H=np.eye(4)/z[4]-C.T@Sigma@C
        return np.r_[z[:4]/z[4]-C.T@pp,np.linalg.eigvalsh(H)[0]]
    start=np.r_[np.array(supplied["fold_theta7"])[[0,2,4,6]],supplied["fold_g_nominal"]]
    solved=root(eq,start,tol=1e-11); g=float(solved.x[4])
    theta=np.zeros(7);theta[[0,2,4,6]]=solved.x[:4]
    p=softmax(X@theta); Sigma=np.diag(p)-np.outer(p,p)
    H=np.eye(7)/g-X.T@Sigma@X
    eigen,vec=np.linalg.eigh(H);v=vec[:,0]
    if v[0]<0:v=-v
    ck("versioned independent fold residual",np.max(abs(eq(solved.x)))<3e-13)
    ck("versioned fold vs source",abs(g-supplied["fold_g_nominal"])<1e-12)
    ck("six retained positive directions",eigen[1]>.17)
    r=Sigma@X@v; ell=g*v@X.T; Q=np.eye(12)-np.outer(r,ell)
    J=g*Sigma@A-np.eye(12)
    ck("right and left null normalization",abs(ell@r-1)<5e-14)
    ck("right null identity",np.linalg.norm(J@r)<5e-14)
    ck("left null identity",np.linalg.norm(ell@J)<5e-14)
    ck("oblique projection",np.linalg.norm(Q@Q-Q)<5e-14)
    # Build a different fast basis by eliminating a pivot, then QR.
    T=np.linalg.qr(np.eye(12)[:,:11]-np.eye(12)[:,[11]])[0]
    a=ell@T;pivot=int(np.argmax(abs(a)))
    axes=[i for i in range(11) if i!=pivot]
    U=np.eye(11)[:,axes].copy()
    U[pivot,:]=-a[axes]/a[pivot]
    E=np.linalg.qr(T@U)[0]; B=E.T@J@E
    P=np.linalg.solve(np.kron(np.eye(10),B.T)+np.kron(B.T,np.eye(10)),-np.eye(10).reshape(-1,order="F")).reshape((10,10),order="F")
    Pe=np.linalg.eigvalsh(P)
    ck("fast basis normalized",np.linalg.norm(E.T@E-np.eye(10))<5e-14 and np.linalg.norm(ell@E)<5e-14)
    ck("ten-mode stable restriction",max(np.linalg.eigvals(B).real)<-.5)
    ck("Kronecker Lyapunov residual",np.linalg.norm(B.T@P+P@B+np.eye(10))<5e-14)
    ck("positive derived Lyapunov form",Pe[0]>.1)
    ck("basis invariant P eigenvalues",np.max(abs(Pe[[0,-1]]-supplied["P_eigen_minmax"]))<2e-12)
    ck("actual leading soft noise",abs(2*g*g*v@X.T@Sigma@X@v-2*g)<5e-13)
    # Verify energy formula by actual finite jumps, not state enumeration.
    energy=[]
    for N in (1200,2400):
        target=p+.4*N**(-1/3)*r
        counts=np.floor(N*target).astype(int)
        rem=N-int(counts.sum()); counts[np.argsort(-(N*target-counts))[:rem]]+=1
        pp=counts/N;w=Q@(pp-p);y=E.T@w
        value=float(y@P@y)
        bN=np.zeros(12);direct=quadratic=variance=0.
        for i in range(12):
            if counts[i]==0:continue
            deleted=counts.copy();deleted[i]-=1
            probabilities=softmax((g/N)*A@deleted)
            for j in range(12):
                if j==i:continue
                rate=counts[i]*probabilities[j]
                jump=np.zeros(12);jump[j]=1/N;jump[i]=-1/N
                dy=E.T@Q@jump
                direct+=rate*(float((y+dy)@P@(y+dy))-value)
                quadratic+=rate*float(dy@P@dy)
                bN+=rate*jump
                dz=N**(1/3)*g*float(v@X.T@jump)
                variance+=rate*dz*dz*N**(1/3)
        formula=2*float(y@P@E.T@Q@bN)+quadratic
        ck(f"actual count preparation N{N}",counts.sum()==N and counts.min()>0 and np.max(abs(pp-target))<=1/N+1e-15)
        ck(f"quadratic generator formula N{N}",abs(direct-formula)<1e-14)
        ck(f"positive jump payment N{N}",quadratic>0)
        energy.append({"N":N,"LV_direct":direct,"LV_formula":formula,
                       "quadratic_jump_term":quadratic,"scaled_soft_variance":variance,
                       "scope":"numerical identity, not a certified Lyapunov constant or limit rate"})
    cert=json.loads((SOURCE/"SF-004/EXIT_SIDE_CERTIFICATE.json").read_text())
    out["exact_exit_certificate"]=exact_exit(cert)
    ck("exact fraction exit certificate",True)
    altered=dict(cert);altered["P_exit_upper_public_outward_interval"]=["0.416855651698382284","0.416855651698382284"]
    try:
        exact_exit(altered)
    except AssertionError:
        rejected=True
    else:
        rejected=False
    ck("reject inward exit endpoint",rejected)
    out["finite_jump_energy_diagnostics"]=energy
    out["numeric_chart"]={"g":g,"P_eigen_minmax":Pe[[0,-1]].tolist(),"stable_max_real":float(max(np.linalg.eigvals(B).real))}
    out["checks"]=checks
    inputs=[SEED,FOLD,plan,current_intake,THREVIEW/"verification.json",Path(__file__),SOURCE/"MANIFEST.sha256"]
    out["input_hashes"]={str(path.relative_to(ROOT)):sha(path) for path in inputs}
    prior_agents=subprocess.run(["git","show","cf1bf9ae:AGENTS.md"],cwd=ROOT,capture_output=True,check=True).stdout
    out["guardrail_baseline"]={"commit":"cf1bf9ae","AGENTS_sha256":hashlib.sha256(prior_agents).hexdigest()}
    cap=Path("/tmp/fin-sf-martingale-primary.md")
    out["literature"]= {"url":"https://arxiv.org/html/0712.1929v2", "read_scope":"Theorem 2.1 and full Theorem 5.3 statement/proof",
                        "capture_sha256":sha(cap) if cap.exists() else None,
                        "application":"UI/joint-filtration martingale stability; deterministic-clock FCLT NOT directly applied to random stopped clock"}
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix="fin-sf-intake-"));copy=temp/"campaign"
        shutil.copytree(SOURCE,copy)
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1")
        logs={}
        # Hash check runs BEFORE the writers: changed binary64 outputs need not match.
        for script in ("verify_package.py","check_fold_chart.py","SF-004/first_passage_certificate.py",
                       "SF-004/independent_first_passage_check.py","independent_verify.py"):
            ran=subprocess.run(["python3",script],cwd=copy,env=env,capture_output=True,text=True,timeout=120)
            logs[script]={"exit_code":ran.returncode,"stdout":ran.stdout,"stderr":ran.stderr}
            assert ran.returncode==0,(script,ran.stderr)
        reproduced=json.loads((copy/"INDEPENDENT_VERIFICATION.json").read_text())
        chart=json.loads((copy/"SF-001/NUMERICAL_CHART.json").read_text())
        out["replay"]={"directory":str(copy),"source_checks":[reproduced["pass_count"],reproduced["total"]],"logs":logs,
                        "source_check_tautology":"rate_is_one_per_copy is an unconditional True contractual flag, not an independent check",
                        "chart_fold_gain_drift":abs(chart["fold_g_nominal"]-supplied["fold_g_nominal"]),
                        "byte_identical":{rel:sha(copy/rel)==sha(SOURCE/rel) for rel in
                                          ("SF-001/NUMERICAL_CHART.json","SF-004/EXIT_SIDE_CERTIFICATE.json",
                                           "SF-004/INDEPENDENT_EXIT_CHECK.json","INDEPENDENT_VERIFICATION.json")}}
    out["source_manifest_after"]=manifest(SOURCE,coverage=True)
    out["status"]="PASS_SCOPED_ARITHMETIC_AND_REPLAY"
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"source":out["source_manifest"],"source_inputs":out["source_inputs"],
                      "task_inputs":sum(x["matched"] for x in out["task_inputs"].values()),
                      "TH_archive":out["embedded_TH_archive"],"checks":len(checks),
                      "replay":out.get("replay",{}).get("source_checks"),"exit":out["exact_exit_certificate"]},indent=2))

if __name__=="__main__":main()
