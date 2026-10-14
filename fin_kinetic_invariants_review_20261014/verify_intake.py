#!/usr/bin/env python3
"""Bounded KI intake checks; analytic stochastic proofs are reviewed separately."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, html, importlib.util, json, math, os, re
import shutil, subprocess, tempfile, zipfile
import numpy as np
from scipy.optimize import root
from scipy.special import softmax, logsumexp, expit

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"fin_kinetic_invariants_campaign_20261014"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def manifest(folder,name="MANIFEST.sha256",coverage=False):
    entries=[]
    for line in (folder/name).read_text().splitlines():
        h,rel=line.split(None,1);rel=str(Path(rel.lstrip("*")))
        assert sha(folder/rel)==h,rel
        entries.append(rel)
    out={"matched":len(entries),"mismatches":0}
    if coverage:
        files={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
        extra=files-set(entries)-{name};assert not extra,extra
        out["uncovered"]=sorted(extra)
    return out
def normalized(text):
    return " ".join(re.sub(r"\\([\\`*_{}\[\]()#+.!><:|\-])",r"\1",html.unescape(text)).replace("**","").split())
def outward(value,places=14,upper=False):
    scale=10**places;num=value*scale
    k=-((-num.numerator)//num.denominator) if upper else num.numerator//num.denominator
    sign="-" if k<0 else "";k=abs(k)
    return f"{sign}{k//scale}.{k%scale:0{places}d}"

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,required=True);parser.add_argument("--replay",action="store_true")
    args=parser.parse_args();checks=[]
    def ck(name,ok,**data):
        checks.append({"name":name,"passed":bool(ok),**data});assert ok,name
    out={"date":"2026-10-14","baseline":"cf1bf9ae plus pre-existing SF integration"}
    out["source_manifest"]=manifest(SOURCE,coverage=True)
    out["source_inputs"]=manifest(SOURCE,"SOURCE_INPUTS.sha256")
    out["task_inputs"]={f"KI-{i:03d}":manifest(SOURCE,f"KI-{i:03d}/INPUTS.sha256") for i in range(1,6)}
    sf=ROOT/"fin_stochastic_fold_closure_review_20261014"
    out["SF_review_manifest"]=manifest(sf,coverage=True)
    previous=json.loads((sf/"verification.json").read_text())
    for rel,h in previous["input_hashes"].items():assert sha(ROOT/rel)==h,rel
    out["inherited_SF_inputs_unchanged"]=len(previous["input_hashes"])
    plan=ROOT/"FIN_KINETIC_INVARIANTS_NEXT_CHAIN_20261014.md"
    attached=SOURCE/"source_context/KI_PLAN_USER.md"
    ck("attached plan normalized text",normalized(plan.read_text())==normalized(attached.read_text()))
    out["plan_binding"]={"repo_sha256":sha(plan),"attached_sha256":sha(attached),"byte_identical":sha(plan)==sha(attached)}
    archive=SOURCE/"source_context/FIN_SF_001_SF_005_FULL_HANDOFF_20261014.zip"
    with zipfile.ZipFile(archive) as z:
        files=[n for n in z.namelist() if not n.endswith("/")]
        mismatch=[n for n in files if not (ROOT/n).exists() or hashlib.sha256(z.read(n)).hexdigest()!=sha(ROOT/n)]
    ck("SF archive members match",not mismatch,matched=len(files),mismatches=mismatch)
    out["SF_archive"]={"matched":len(files),"mismatches":mismatch}
    ck("producer reconstruction core byte bound",sha(SOURCE/"source_context/theory_core.py")==sha(ROOT/"fin_intrinsic_theory_campaign_20261014/theory_core.py"))
    ck("E04 matched-clock contract byte bound",sha(SOURCE/"source_context/E04_KINETIC_CONTRACTS.json")==sha(ROOT/"fin_generative_operational_campaign_20261014/RG-E04/KINETIC_CONTRACTS.json"))
    spec=importlib.util.spec_from_file_location("versioned_seed",SEED);seed=importlib.util.module_from_spec(spec);spec.loader.exec_module(seed)
    _,_,lam,X,A=seed.build_rank7()
    original=json.loads((SOURCE/"KI-002/DERIVED_SOFT_CLOCKS.json").read_text())
    C=X[:,[0,2,4,6]]
    def fold_eq(z):
        p=softmax(C@z[:4]);Sigma=np.diag(p)-np.outer(p,p)
        return np.r_[z[:4]/z[4]-C.T@p,np.linalg.eigvalsh(np.eye(4)/z[4]-C.T@Sigma@C)[0]]
    start=np.array([1.3643114281782407,1.4331027014285425,1.4080457012550374,1.0056820773201627,original["fold_g_nominal"]])
    sol=root(fold_eq,start,tol=1e-11);g=float(sol.x[4]);theta=np.zeros(7);theta[[0,2,4,6]]=sol.x[:4]
    p=softmax(X@theta);Sigma=np.diag(p)-np.outer(p,p)
    H7=np.eye(7)/g-X.T@Sigma@X;_,V=np.linalg.eigh(H7);v=V[:,0]
    if v[0]<0:v=-v
    r=Sigma@X@v;HF=np.diag(1/p)-g*A
    T=np.linalg.qr(np.eye(12)[:,:11]-np.eye(12)[:,[11]])[0]
    He=np.linalg.eigvalsh(T.T@HF@T)
    a=-float(v@theta)/g**2;b=-float(np.sum(p*((X@v)-p@(X@v))**3))
    ck("independent fold residual",np.max(abs(fold_eq(sol.x)))<3e-13)
    ck("primal tangent ten positive eigenvalues",He[1]>12.,eigenvalues=He.tolist())
    ck("primal Hessian null",np.linalg.norm(T.T@HF@r)<1e-13)
    ck("primal dual a identity",abs(-r@A@p-a)<1e-13)
    ck("primal dual b identity",abs(-np.sum(r**3/p**2)-b)<1e-13)
    def rates(pp,gg,kin):
        h=gg*A@pp;d=h[None,:]-h[:,None]
        if kin=="HB":return np.tile(softmax(h),(12,1))
        if kin=="MET":return np.exp(np.minimum(0,d))/12
        return expit(d)/6
    def drift(pp,gg,kin):
        c=pp[:,None]*rates(pp,gg,kin);return c.sum(0)-c.sum(1)
    def matrixK(c,logarithmic=False):
        K=np.zeros((12,12))
        for i in range(12):
            for j in range(i+1,12):
                value=c[i,j]
                if logarithmic:
                    z=math.log(c[i,j]/c[j,i])
                    value=c[j,i]*(math.expm1(z)/z if z else 1.)
                e=np.eye(12)[j]-np.eye(12)[i];K+=value*np.outer(e,e)
        return K
    numeric={};finite_cases=[]
    for kin in ("HB","MET","BAR"):
        c=p[:,None]*rates(p,g,kin);K=matrixK(c);Kt=T.T@K@T
        ki=T@np.linalg.solve(Kt,T.T@r);mob=1/float(r@ki);ell=mob*ki;J=-K@HF
        ck(kin+" strictly positive tangent mobility",np.linalg.eigvalsh(Kt)[0]>.001)
        ck(kin+" equilibrium conductance symmetry",np.max(abs(c-c.T))<1e-13)
        ck(kin+" source soft clock",abs(mob-original["models"][kin]["mobility_m_nominal"])<2e-12)
        D=np.zeros_like(K)
        for i in range(12):
            for j in range(12):
                e=np.eye(12)[j]-np.eye(12)[i];D+=c[i,j]*np.outer(e,e)
        ck(kin+" equilibrium D=2K",np.max(abs(D-2*K))<1e-13)
        ck(kin+" projected variance=2m",abs(ell@D@ell-2*mob)<1e-12)
        aa=ell@T;pivot=int(np.argmax(abs(aa)));axes=[i for i in range(11) if i!=pivot]
        U=np.eye(11)[:,axes].copy();U[pivot,:]=-aa[axes]/aa[pivot]
        E=np.linalg.qr(T@U)[0];B=E.T@J@E
        ck(kin+" actual fast restricted spectrum",max(np.linalg.eigvals(B).real)<-.3)
        h=1e-7
        derivative=np.column_stack([(drift(p+h*T[:,i],g,kin)-drift(p-h*T[:,i],g,kin))/(2*h) for i in range(11)])
        ck(kin+" derivative of actual drift vs -KHF",np.max(abs(derivative-J@T))<2e-5,residual=float(np.max(abs(derivative-J@T))))
        test=p+1e-3*np.linspace(-1,1,12);ct=test[:,None]*rates(test,g,kin)
        KL=matrixK(ct,True);force=np.log(test)-g*A@test
        ck(kin+" nonlinear Onsager identity",np.linalg.norm(drift(test,g,kin)+KL@force)<1e-13)
        g0=rates(p,0.,kin);mask=~np.eye(12,dtype=bool)
        ck(kin+" actual g0 rates",np.max(abs(g0[mask]-1/12))<1e-15)
        # Independent finite jump expansion for V^2: positive jump payment retained.
        P=np.linalg.solve(np.kron(np.eye(10),B.T)+np.kron(B.T,np.eye(10)),-np.eye(10).reshape(-1,order="F")).reshape((10,10),order="F")
        N=1200;target=p+.4*N**(-1/3)*r;counts=np.floor(N*target).astype(int)
        rem=N-int(counts.sum());counts[np.argsort(-(N*target-counts))[:rem]]+=1
        pp=counts/N;Q=np.eye(12)-np.outer(r,ell);y=E.T@Q@(pp-p);value=float(y@P@y)
        LV=LV2=payment=0.
        for i in range(12):
            deleted=counts.copy();deleted[i]-=1;field=(g/N)*A@deleted
            probs=softmax(field) if kin=="HB" else np.exp(np.minimum(0,field-field[i]))/12 if kin=="MET" else expit(field-field[i])/6
            for j in range(12):
                if j==i:continue
                delta=np.zeros(12);delta[j]=1/N;delta[i]=-1/N
                dy=E.T@Q@delta;new=float((y+dy)@P@(y+dy));rate=counts[i]*probs[j]
                LV+=rate*(new-value);LV2+=rate*(new*new-value*value);payment+=rate*(new-value)**2
        ck(kin+" exact quadratic-square generator identity numeric",abs(LV2-2*value*LV-payment)<1e-15 and payment>0)
        numeric[kin]={"m":mob,"m_over_HB":mob/g,"d":(4/b)**(2/3)/(2*mob),
                      "LV2":LV2,"LV2_positive_payment":payment,"scope":"floating diagnostic; no outward moment constants"}
    ck("HB algebraic clock recovery numerical",abs(numeric["HB"]["m"]-g)<2e-13)
    ck("distinct clock hierarchy",numeric["HB"]["m"]>numeric["BAR"]["m"]>numeric["MET"]["m"])
    ck("reject unmatched Barker normalization",abs(expit(0)/12-1/12)>.04)
    # Finite count DB, using pairwise log weight ratios rather than lgamma sums.
    counts=np.full(12,3);counts[0]+=1;counts[1]-=1;N=int(counts.sum())
    def rate(nn,i,j,kin):
        deleted=nn.copy();deleted[i]-=1;field=(g/N)*A@deleted;delta=field[j]-field[i]
        rr=softmax(field)[j] if kin=="HB" else math.exp(min(0.,delta))/12 if kin=="MET" else expit(delta)/6
        return nn[i]*rr
    for kin in ("HB","MET","BAR"):
        errors=[]
        for i in range(12):
            for j in range(12):
                if i==j:continue
                moved=counts.copy();moved[i]-=1;moved[j]+=1
                logratio=math.log(counts[i]/(counts[j]+1))+(g/(2*N))*float(moved@A@moved-counts@A@counts)
                errors.append(abs(math.log(rate(counts,i,j,kin)/rate(moved,j,i,kin))-logratio))
        ck(kin+" finite count detailed balance",max(errors)<1e-12)
        finite_cases.append({"kinetic":kin,"cases":len(errors),"max_log_discrepancy":max(errors)})
    out["nominal_mobilities"]=numeric;out["finite_DB"]=finite_cases
    out["fold_free_energy_nominal"]=float(theta@theta/(2*g)-logsumexp(X@theta)+math.log(12))
    old_spectrum=json.loads((ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json").read_text())["exact"]["laplacian_intervals"]
    lo6,hi6=map(F,old_spectrum[6]);lo5,hi5=map(F,old_spectrum[5])
    gc_bounds=[outward(12/hi6),outward(12/lo6,upper=True)]
    gap=1-hi5/lo6
    ck("uniform critical gain and positive fast gap",gap>0)
    out["next_lane_derived_spectral_bounds"]={"g_uniform_critical":gc_bounds,"retained_fast_rate_lower":outward(gap),
                                             "scope":"exact Fraction propagation of inherited spectral boxes; not a new critical stochastic theorem"}
    out["checks"]=checks
    inputs=[SEED,plan,sf/"MANIFEST.sha256",sf/"STOPPED_LIMIT_COMPLETION.md",SOURCE/"MANIFEST.sha256",Path(__file__)]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["preexisting_guardrail"]={"description":"AGENTS.md before the KI append, including preserved dirty SF section",
                                  "sha256":"aad8636e9953ad538db3ece550822db508eabf43556045fce11839d1ffc12887"}
    capture=Path("/tmp/fin-ki-uniform-benchmark.md")
    out["benchmark"]={"url":"https://arxiv.org/html/0712.0790v1#S1.Thmtheorem2",
                      "read_scope":"update definition and full Theorem2 statement; not a FIN application",
                      "capture_sha256":sha(capture) if capture.exists() else None}
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix="fin-ki-intake-"));copy=temp/"campaign"
        shutil.copytree(SOURCE,copy,ignore=shutil.ignore_patterns("__pycache__"))
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1")
        logs={}
        for script in ("verify_package.py","kinetic_analysis.py","local_regularities.py","independent_verify.py"):
            ran=subprocess.run(["python3",script],cwd=copy,env=env,capture_output=True,text=True,timeout=120)
            logs[script]={"exit_code":ran.returncode,"stdout":ran.stdout,"stderr":ran.stderr};assert ran.returncode==0,(script,ran.stderr)
        ind=json.loads((copy/"INDEPENDENT_VERIFICATION.json").read_text())
        regenerated=json.loads((copy/"KI-002/DERIVED_SOFT_CLOCKS.json").read_text())
        out["replay"]={"directory":str(copy),"source_flags":[ind["passed"],ind["total"]],"logs":logs,
                       "mobility_max_drift":max(abs(original["models"][k]["mobility_m_nominal"]-regenerated["models"][k]["mobility_m_nominal"]) for k in numeric),
                       "byte_identical":{rel:sha(copy/rel)==sha(SOURCE/rel) for rel in ("KI-002/DERIVED_SOFT_CLOCKS.json","KI-003/LOCAL_DRIFT_DIAGNOSTICS.json","INDEPENDENT_VERIFICATION.json")}}
    out["source_manifest_after"]=manifest(SOURCE,coverage=True)
    out["status"]="PASS_SCOPED_KI_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"source":out["source_manifest"],"source_inputs":out["source_inputs"],"stage_inputs":sum(v["matched"] for v in out["task_inputs"].values()),
                      "SF_archive":out["SF_archive"],"checks":len(checks),"source_flags":out.get("replay",{}).get("source_flags"),
                      "mobilities":{k:v["m"] for k,v in numeric.items()},"fold_free_energy":out["fold_free_energy_nominal"],"uniform_critical":out["next_lane_derived_spectral_bounds"]},indent=2))

if __name__=="__main__":main()
