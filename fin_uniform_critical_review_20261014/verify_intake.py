#!/usr/bin/env python3
"""Small independent checks and source replay; analytical proofs are separate."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,html,importlib.util,json,math,os,re,shutil,subprocess,tempfile,zipfile
import numpy as np
from scipy.special import softmax,expit,gamma
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"fin_uniform_critical_campaign_20261014"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def manifest(folder,name="MANIFEST.sha256",coverage=False):
    files=[]
    for line in (folder/name).read_text().splitlines():
        h,rel=line.split(None,1);rel=str(Path(rel.lstrip("*")));assert sha(folder/rel)==h,rel;files.append(rel)
    out={"matched":len(files),"mismatches":0}
    if coverage:
        actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
        extra=actual-set(files)-{name};assert not extra,extra;out["uncovered"]=sorted(extra)
    return out
def norm(text):return " ".join(re.sub(r"\\([\\`*_{}\[\]()#+.!><:|\-])",r"\1",html.unescape(text)).replace("**","").split())

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,required=True);parser.add_argument("--replay",action="store_true");args=parser.parse_args()
    checks=[]
    def ck(name,condition,**data):
        checks.append({"name":name,"passed":bool(condition),**data});assert condition,name
    out={"date":"2026-10-14","audit_start":"cf1bf9ae plus SF/KI working-tree integrations",
         "baseline":"21a875c5, observed during audit; source inputs unchanged",
         "preexisting_AGENTS_sha256":"c260aa6d17bb1f0bfdee5a308dd65ed1b47c2c7957e2b6748f338a737c977059"}
    out["source_manifest"]=manifest(SOURCE,coverage=True)
    out["source_inputs"]=manifest(SOURCE,"SOURCE_INPUTS.sha256")
    out["task_inputs"]={f"UF-{i:03d}":manifest(SOURCE/f"UF-{i:03d}","INPUTS.sha256") for i in range(1,6)}
    kr=ROOT/"fin_kinetic_invariants_review_20261014"
    out["KI_review_manifest"]=manifest(kr,coverage=True)
    previous=json.loads((kr/"verification.json").read_text())
    for rel,h in previous["input_hashes"].items():assert sha(ROOT/rel)==h,rel
    out["KI_verification_inputs_unchanged"]=len(previous["input_hashes"])
    for local,attached in ((ROOT/"FIN_UNIFORM_CRITICAL_NEXT_CHAIN_20261014.md",SOURCE/"PLAN_UF_USER.md"),(kr/"INTAKE_20261014.md",SOURCE/"KI_INTAKE_USER.md")):
        ck(local.name+" normalized attached text",norm(local.read_text())==norm(attached.read_text()),
           repo_sha256=sha(local),attached_sha256=sha(attached),byte_identical=sha(local)==sha(attached))
    with zipfile.ZipFile(SOURCE/"SOURCE_KI_FULL_HANDOFF.zip") as archive:
        names=[n for n in archive.namelist() if not n.endswith("/")]
        bad=[n for n in names if not (ROOT/n).exists() or hashlib.sha256(archive.read(n)).hexdigest()!=sha(ROOT/n)]
    ck("KI archive file members",not bad,matched=len(names),mismatches=bad)
    out["KI_archive"]={"matched":len(names),"mismatches":bad}
    ck("bundled reconstruction code binding",sha(SOURCE/"source_context/theory_core.py")==sha(ROOT/"fin_intrinsic_theory_campaign_20261014/theory_core.py"))
    spec=importlib.util.spec_from_file_location("versioned_seed",SEED);seed=importlib.util.module_from_spec(spec);spec.loader.exec_module(seed)
    _,_,lam,X,A=seed.build_rank7();sig=(-1.)**np.arange(12);u=np.ones(12)/12;chi=lam[6]/12;gc=1/chi
    Q=np.eye(12)-np.ones((12,12))/12-np.outer(sig,sig)/12
    vec=np.linalg.eigh(Q)[1][:,-10:];J=gc*A/12-np.eye(12);Be=vec.T@J@vec
    ck("full ten-fast spectrum",np.linalg.eigvalsh(Be)[-1]<-.018)
    ck("orthogonal fast projector",np.linalg.norm(Q@Q-Q)<1e-13 and np.linalg.norm(Q@sig)<1e-13)
    ck("exact critical parity linear mode numeric",np.linalg.norm(J@sig)<1e-13)
    def rates(pp,g,kind,N=None):
        out=np.zeros((12,12))
        for i in range(12):
            field=g*A@(pp-(np.eye(12)[i]/N if N else 0))
            out[i]=softmax(field) if kind=="HB" else np.exp(np.minimum(0,field-field[i]))/12 if kind=="MET" else expit(field-field[i])/6
            out[i,i]=0
        return out
    def drift(pp,g,kind):
        c=pp[:,None]*rates(pp,g,kind);return c.sum(0)-c.sum(1)
    rule_checks={}
    for kind in ("HB","MET","BAR"):
        ru=rates(u,gc,kind);mask=~np.eye(12,dtype=bool)
        ck(kind+" actual uniform and g0 rates",np.max(abs(ru[mask]-1/12))<1e-14 and np.max(abs(rates(u,0,kind)[mask]-1/12))<1e-14)
        coeff=float(np.sum(u[:,None]*ru*(sig[None,:]-sig[:,None])**2))
        ck(kind+" parity bracket two",abs(coeff-2)<1e-13)
        for m in (-.01,.01):
            pm=u+m*sig/12;actual=float(sig@drift(pm,gc,kind))
            formula=math.tanh(m)-m if kind!="MET" else .5*((1-m)*min(1,math.exp(2*m))-(1+m)*min(1,math.exp(-2*m)))
            ck(kind+f" exact parity-line formula {m}",abs(actual-formula)<2e-15)
            h=1e-6
            response=max(abs(sig@(drift(pm+h*vec[:,i],gc,kind)-drift(pm-h*vec[:,i],gc,kind))/(2*h)) for i in range(10))
            ck(kind+f" fast derivative cancellation {m}",response<3e-8,residual=float(response))
        # Actual finite-jump energy, drift and bracket on a genuine N-type.
        n=np.full(12,120);n[0]+=5;n[2]-=5;N=int(n.sum());pp=n/N;rr=rates(pp,gc,kind,N);w=Q@(pp-u)
        value=float(w@w);LV=payment=0.;bm=qm=0.;bN=np.zeros(12)
        for i in range(12):
            for j in range(12):
                if j==i:continue
                rate=n[i]*rr[i,j];jump=(np.eye(12)[j]-np.eye(12)[i])/N;dw=Q@jump
                LV+=rate*(float((w+dw)@(w+dw))-value);payment+=rate*float(dw@dw);bN+=rate*jump
                dm=float(sig@jump);bm+=rate*dm;qm+=rate*dm*dm
        ck(kind+" real-jump energy identity",abs(LV-2*w@Q@bN-payment)<1e-14 and payment>0)
        supplied=json.loads((SOURCE/"DIAGNOSTICS.json").read_text())["rules"][kind]["finite_loo_diagnostic"]
        ck(kind+" finite LOO drift replay",abs(bm-supplied["drift"])<1e-14)
        ck(kind+" finite LOO bracket replay",abs(qm-supplied["bracket_rate_m"])<1e-15)
        rule_checks[kind]={"bracket_at_uniform":coeff,"LV_actual":LV,"positive_jump_payment":payment}
    # Check gamma-law moments by a separate deterministic quadrature.
    z=quad(lambda x:math.exp(-x**4/12),-np.inf,np.inf)[0]
    moment=quad(lambda x:x*x*math.exp(-x**4/12),-np.inf,np.inf)[0]/z
    exact_formula=math.sqrt(12)*gamma(.75)/gamma(.25)
    ck("quartic variance integral vs Gamma",abs(moment-exact_formula)<1e-11)
    ck("quartic variance strictly positive",exact_formula>1.)
    # The iid parity fourth moment: independent enumeration for small N, no Monte Carlo.
    for N in (3,8):
        fourth=sum(math.comb(N,k)*(2*k-N)**4 for k in range(N+1))/(2**N*N**4)
        ck(f"exact initial parity fourth moment N{N}",F(sum(math.comb(N,k)*(2*k-N)**4 for k in range(N+1)),2**N*N**4)==F(3*N*N-2*N,N**4))
    from check_gibbs_certificate import build,check_public
    cert=build();check_public(json.loads(Path(__file__).with_name("GIBBS_EXCLUSION_CERTIFICATE.json").read_text()),cert)
    ck("exact rational Gibbs exclusion payments",True)
    altered=dict(cert);altered["exclusion_exponent"]="3/5"
    try:check_public(altered,cert)
    except AssertionError:rejected=True
    else:rejected=False
    ck("reject stronger unpaid Gibbs exponent",rejected)
    out["checks"]=checks;out["actual_rule_checks"]=rule_checks
    out["quartic_variance_nominal"]=float(exact_formula)
    out["global_exclusion"]={"bound":cert["bound"],"conditional_on":"accepted exact input spectral boxes"}
    inputs=[SEED,ROOT/"FIN_UNIFORM_CRITICAL_NEXT_CHAIN_20261014.md",kr/"MANIFEST.sha256",ROOT/"fin_stochastic_fold_closure_review_20261014/STOPPED_LIMIT_COMPLETION.md",SOURCE/"MANIFEST.sha256",Path(__file__),Path(__file__).with_name("check_gibbs_certificate.py")]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    preexisting=subprocess.run(["git","show","21a875c5:AGENTS.md"],cwd=ROOT,capture_output=True,check=True).stdout
    assert hashlib.sha256(preexisting).hexdigest()==out["preexisting_AGENTS_sha256"]
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix="fin-uf-intake-"));copy=temp/"campaign";shutil.copytree(SOURCE,copy,ignore=shutil.ignore_patterns("__pycache__"))
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1");logs={}
        for script in ("uf_analysis.py","independent_verify.py"):
            ran=subprocess.run(["python3",script],cwd=copy,env=env,capture_output=True,text=True,timeout=120)
            logs[script]={"exit_code":ran.returncode,"stdout":ran.stdout,"stderr":ran.stderr};assert ran.returncode==0,(script,ran.stderr)
        ind=json.loads((copy/"INDEPENDENT_VERIFICATION.json").read_text())
        byte_equal={rel:sha(copy/rel)==sha(SOURCE/rel) for rel in ("DIAGNOSTICS.json","INDEPENDENT_VERIFICATION.json")}
        # Producer validator is a same-environment byte gate, not an independent theorem check.
        ran=subprocess.run(["python3","validate_package.py"],cwd=copy,env=env,capture_output=True,text=True,timeout=120)
        assert ran.returncode==0,ran.stderr
        logs["validate_package.py"]={"exit_code":ran.returncode,"stdout":ran.stdout,"stderr":ran.stderr,"scope":"fresh numerical instance replay, not original cross-environment byte custody"}
        out["replay"]={"directory":str(copy),"source_flags":[ind["passed"],ind["total"]],"byte_identical_to_incoming":byte_equal,"logs":logs}
    out["source_manifest_after"]=manifest(SOURCE,coverage=True)
    out["status"]="PASS_SCOPED_UF_ARITHMETIC_AND_REPLAY"
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"source":out["source_manifest"],"source_inputs":out["source_inputs"],"task_inputs":sum(v["matched"] for v in out["task_inputs"].values()),"KI_archive":out["KI_archive"],"checks":len(checks),"source_flags":out.get("replay",{}).get("source_flags"),"Gibbs_bound":cert["bound"]},indent=2))

if __name__=="__main__":main()
