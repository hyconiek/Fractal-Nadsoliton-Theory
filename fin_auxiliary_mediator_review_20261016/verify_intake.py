#!/usr/bin/env python3
"""HM intake: independent finite checks and copy-only source replay.

No finite diagnostic is a proof of the stationary critical process limit.
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
import zipfile
import numpy as np
from scipy.special import softmax
import sympy as sy

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/"fin_auxiliary_mediator_campaign_20261016"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def manifest(folder,excluded=()):
    names=[]
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=line.split("  ",1); assert name not in names and sha(folder/name)==h,name; names.append(name)
    actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    assert actual==set(names)|{"MANIFEST.sha256"}|set(excluded),(folder,actual-set(names))
    return {"matched":len(names),"uncovered":[],"separately_bound":{n:sha(folder/n) for n in excluded}}

def action(L,v,t,rate,cutoff=85):
    P=np.eye(len(L))+L/rate; assert P.min()>-2e-14
    term=v.copy(); weight=math.exp(-rate*t); result=weight*term
    for j in range(1,cutoff+1):
        term=P@term; weight*=rate*t/j; result+=weight*term
    return result

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path,required=True); parser.add_argument("--replay",action="store_true")
    args=parser.parse_args(); checks=[]
    def ck(name,ok,**detail):
        checks.append({"name":name,"passed":bool(ok),**detail}); assert ok,name
    out={"date":"2026-10-16","baseline":"3be481a7 plus preserved NT intake/AGENTS/plan; incoming HM source untracked"}
    out["source_manifest"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"])
    records=json.loads((SOURCE/"INPUT_PROVENANCE.json").read_text())["available"]
    originals=[r for r in records if r.get("sha256")]
    expected="".join(r["sha256"]+"  "+r["path"]+"\n" for r in originals)
    ck("six original inputs and all five stage ledgers agree",(SOURCE/"INPUTS.sha256").read_text()==expected and
       all((SOURCE/f"HM-{i:03d}/INPUTS.sha256").read_text()==expected for i in range(1,6)))
    mapping={}
    for row in originals:
        if row["description"]=="MT archive":rel="fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip"
        elif row["description"]=="UF archive":rel="fin_critical_metastability_campaign_20261014/source_context/UF_PREDECESSOR.zip"
        else:continue
        assert sha(ROOT/rel)==row["sha256"]
        mapping[row["path"]]=rel
    members={}
    for row in records:
        if "member" not in row:continue
        rel=row["member"]; p=ROOT/rel
        assert sha(p)==row["sha256_member_bytes"],rel
        entry={"local_file":rel,"sha256":sha(p),"wrapper_verified":row["archive"] in mapping}
        if entry["wrapper_verified"]:
            with zipfile.ZipFile(ROOT/mapping[row["archive"]]) as z:
                assert hashlib.sha256(z.read(rel)).hexdigest()==sha(p)
        members[rel]=entry
    out["provenance"]={"matched_wrappers":mapping,"unavailable_originals":[r["path"] for r in originals if r["path"] not in mapping],
       "matched_extracted_members":members,"stage_reference_counts":{"matched":5*len(mapping),"unavailable_original":5*(len(originals)-len(mapping))}}
    prior={}
    for rel in ("fin_native_transfer_review_20261016","fin_critical_metastability_review_20261015",
                "fin_kinetic_invariants_review_20261014","fin_dynamic_kernel_observability_review_20261015"):
        folder=ROOT/rel; entry=manifest(folder); saved=json.loads((folder/"verification.json").read_text())
        for name,h in saved["input_hashes"].items():assert sha(ROOT/name)==h,name
        entry["saved_inputs_unchanged"]=len(saved["input_hashes"]); prior[rel]=entry
    out["actual_dependencies"]=prior
    spec=importlib.util.spec_from_file_location("versioned_seed",SEED)
    seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,X,A=seed.build_rank7(); g=12/lam[6]; sigma=(-1.)**np.arange(12)
    rng=np.random.default_rng(20261016)
    for N in (1,2,12,24):
        worst=0.
        for repetition in range(4):
            r=np.bincount(rng.integers(0,12,size=N-1),minlength=12)
            q=softmax(g/N*(A@r)); base=X.T@r
            completed=softmax(np.array([g/(2*N)*np.dot(base+X[j],base+X[j]) for j in range(12)]))
            worst=max(worst,float(np.max(abs(q-completed))))
        ck(f"independent completed-square cavity N{N}",worst<4e-15,maximum_error=worst,scope="finite proposal fixture, not MT")
    # Negative control: a nonconstant diagonal needs the explicit self term.
    B=A.copy(); B[0,0]+=.1; r=np.zeros(12); r[2]=1
    native=softmax(g/2*(B@r)); fixed=softmax(g/2*(B@r)+g/4*np.diag(B))
    completed=softmax(np.array([g/4*(r+np.eye(12)[j])@B@(r+np.eye(12)[j]) for j in range(12)]))
    ck("nonconstant-diagonal shortcut rejected",np.max(abs(native-completed))>1e-4 and np.max(abs(fixed-completed))<2e-15)
    # General Bayes/posterior retention identity, exact rational two-atom fixture.
    p=[F(1,8),F(1,24)]; mean=sum(p)/2; posterior=sum(v*v for v in p)/(2*mean)
    variance=sum((v-mean)**2 for v in p)/2
    ck("exact posterior-retention identity",posterior-mean==variance/mean==F(1,48),scope="two-atom algebra; Gaussian strictness is analytic full-support argument")
    # Joint conditional DB at two actual finite type transitions and field points.
    n=np.zeros(12,dtype=int); n[0]=2
    n1=n.copy(); n1[0]-=1; n1[1]+=1
    for index,theta in enumerate((np.zeros(7),rng.normal(size=7)/3)):
        q=softmax(X@theta)
        w0=math.exp((X.T@n)@theta); w1=2*math.exp((X.T@n1)@theta)
        ck(f"persistent conditional count DB field{index}",abs(w0*2*q[1]-w1*q[0])<3e-15)
    k,Ns,gs,mu,t=sy.symbols("k N g mu t",positive=True)
    density=sy.exp(-Ns*(t-mu)**2/(2*gs))
    residual=-sy.diff(k*(mu-t)*density,t)+k*gs/Ns*sy.diff(density,t,2)
    ck("OU conditional Gaussian stationary identity",sy.simplify(residual)==0)
    # Nonzero reversible compression fixture; producer's commuting toy has zero defect.
    Ltoy=sy.Matrix([[-3,2,1,0],[2,-5,0,3],[1,0,-3,2],[0,3,2,-5]])
    P=sy.diag(sy.ones(2)/2,sy.ones(2)/2); Q=sy.eye(4)-P
    f=sy.Matrix([1,1,-1,-1]); v=Q*Ltoy*f
    defect=(f.T*(P*Ltoy**2*P-(P*Ltoy*P)**2)*f)[0]/4
    ck("exact positive compression defect",defect==(v.T*v)[0]/4==4)
    ck("A7 parity field genuinely nonconstant",abs(sigma@softmax(X@np.r_[np.zeros(6),1.]))>.1)
    # Full 18-dimensional local Jacobian rather than repeated scalar polynomials.
    E=np.linalg.qr(np.eye(12)[:,:11]-np.eye(12)[:,[11]])[0]
    for kappa in (1,4):
        J=np.block([[-np.eye(11),E.T@X/12],[kappa*g*X.T@E,-kappa*np.eye(7)]])
        vals=np.linalg.eigvals(J)
        ck(f"actual full tangent inertia kappa{kappa}",sum(abs(vals)<1e-12)==1 and sum(vals.real< -1e-4)==17)
        left=np.r_[kappa/(1+kappa)*(sigma@E),np.zeros(6),math.sqrt(lam[6]/12)/(1+kappa)]
        right=np.r_[E.T@sigma/12,g*X.T@sigma/12]
        D=np.diag(np.r_[np.full(11,2/12),np.full(7,2*kappa*g)])
        ck(f"18-coordinate soft mobility noise kappa{kappa}",np.max(abs(left@J))<1e-13 and abs(left@right-1)<1e-13 and
           abs(left@D@left-2*kappa/(1+kappa))<2e-13)
    # Exact-native reset lift checked on a finite three-atom residual kernel.
    # This is a first/second-moment algebra fixture, not a Gaussian discretization proof.
    states=seed.compositions(2); L=seed.count_generator(states,g,A).toarray(); pi=seed.count_pi(states,g,A)
    weights=np.array([1/6,2/3,1/6]); residual_values=np.array([-math.sqrt(3),0.,math.sqrt(3)])
    refresh=np.tile(weights,(3,1)); off=L-np.diag(np.diag(L))
    incidence=np.kron(np.eye(len(states)),np.ones((3,1))); jointpi=np.kron(pi,weights)
    hidden=np.tile(residual_values,len(states))*math.sqrt(g/2); reset_cov=[]
    for eta in (0.,2.):
        lifted=np.kron(off,refresh)+np.kron(np.diag(np.diag(L)),np.eye(3))+np.kron(np.eye(len(states)),eta*(refresh-np.eye(3)))
        ck(f"reset count intertwining eta{eta:g}",np.max(abs(lifted@incidence-incidence@L))<4e-14)
        ck(f"reset joint DB eta{eta:g}",np.max(abs(jointpi[:,None]*lifted-(jointpi[:,None]*lifted).T))<3e-16)
        observed=jointpi*hidden@action(lifted,hidden,.37,rate=2+eta)
        predicted=g/2*math.exp(-eta*.37)*np.dot(pi,np.exp(np.diag(L)*.37))
        ck(f"reset hidden covariance eta{eta:g}",abs(observed-predicted)<3e-14,scope="three-atom residual moment fixture")
        reset_cov.append(float(observed))
    ck("hidden covariance differs with identical count law",reset_cov[0]>reset_cov[1]>0)
    # A symbolic seed for the next memory-versus-readout study, not its native theorem.
    a,b,x,y,q=sy.symbols("a b x y q",real=True)
    def C(j):return a*x**j+b*y**j
    def diff(j):return C(j+1)-q*C(j)
    hankel=sy.expand(diff(1)*diff(3)-diff(2)**2)
    expected_hankel=a*b*x*y*(x-q)*(y-q)*(x-y)**2
    ck("reference-annihilated Hankel algebra",sy.expand(hankel-expected_hankel)==0)
    for kappa in (1,4):
        gamma=1-lam[3]/lam[6]; c=1-gamma
        J=np.array([[-1.,1.],[kappa*c,-kappa]])
        covariance=np.array([[1.,c],[c,c]])/gamma
        ck(f"next fast Gaussian stationary covariance seed kappa{kappa}",np.max(abs(J@covariance+covariance@J.T+np.diag([2.,2*kappa*c])))<3e-14)
    if args.replay:
        tmp=Path(tempfile.mkdtemp(prefix="fin-hm-intake-"))/"campaign"
        shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        assert manifest(tmp,["PACKAGE_VERIFICATION.json"])==out["source_manifest"]
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1"); logs={}
        for script in ("package_check.py","verify_hm.py"):
            run=subprocess.run(["python3",script],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
            assert run.returncode==0,(script,run.stdout,run.stderr)
            logs[script]={"exit_code":run.returncode,"stdout":run.stdout,"stderr":run.stderr}
        fresh=json.loads((tmp/"verification.json").read_text())
        original=json.loads((SOURCE/"verification.json").read_text()); changes=[]
        def compare(a,b,path=""):
            assert type(a) is type(b),path
            if isinstance(a,dict):
                assert a.keys()==b.keys(),path
                for key in a:compare(a[key],b[key],path+"/"+key)
            elif isinstance(a,list):
                assert len(a)==len(b),path
                for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+"/"+str(i))
            elif a!=b:
                assert isinstance(a,(float,int)) and not isinstance(a,bool),(path,a,b)
                changes.append({"path":path,"absolute_drift":abs(a-b)})
        compare(original,fresh)
        out["replay"]={"directory":str(tmp),"logs":logs,"source_checks":fresh["checks_count"],
            "numeric_output_byte_equal":sha(tmp/"verification.json")==sha(SOURCE/"verification.json"),
            "numeric_changes":changes,"maximum_absolute_drift":max((r["absolute_drift"] for r in changes),default=0.)}
    inputs=[SEED,Path(__file__),SOURCE/"MANIFEST.sha256",ROOT/"FIN_AUXILIARY_MEDIATOR_NEXT_CHAIN_20261016.md"]
    inputs += [ROOT/rel/"MANIFEST.sha256" for rel in prior]
    inputs += [ROOT/rel for rel in (
        "fin_rank7_mathphysics_next/proofs/MP7-032_gaussian_auxiliary_representation.md",
        "fin_native_transfer_review_20261016/MEDIATOR_BRIDGE_SEED.md",
        "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
        "fin_critical_metastability_review_20261015/INTAKE_20261015.md",
        "fin_dynamic_kernel_observability_review_20261015/PROOF_COMPLETIONS.md")]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["checks"]=checks; out["status"]="PASS_SCOPED_HM_INTAKE_CHECKS"
    out["source_manifest_after"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"])
    assert out["source_manifest_after"]==out["source_manifest"]
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"source":out["source_manifest"],"wrapper_inputs":len(mapping),
        "extracted_members":len(members),"checks":len(checks),"source_checks":out.get("replay",{}).get("source_checks")},indent=2))

if __name__=="__main__":main()
