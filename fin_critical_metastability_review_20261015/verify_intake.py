#!/usr/bin/env python3
"""Small arithmetic/replay checks; no stochastic theorem or all-type enumeration."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,html,importlib.util,json,math,os,re,shutil,subprocess,tempfile,zipfile
import numpy as np
from scipy.special import softmax,expit

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"fin_critical_metastability_campaign_20261014"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder,name="MANIFEST.sha256",coverage=False):
    entries=[]
    for line in (folder/name).read_text().splitlines():
        h,rel=line.split(None,1);rel=str(Path(rel.lstrip("*")));assert sha(folder/rel)==h,rel;entries.append(rel)
    out={"matched":len(entries),"mismatches":0}
    if coverage:
        actual={str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
        extra=actual-set(entries)-{name};assert not extra,extra;out["uncovered"]=sorted(extra)
    return out
def norm(s):return " ".join(re.sub(r"\\([\\`*_{}\[\]()#+.!><=:|~\-])",r"\1",html.unescape(s)).replace("**","").split())

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,required=True);parser.add_argument("--replay",action="store_true");args=parser.parse_args()
    checks=[]
    def ck(name,ok,**values):
        checks.append({"name":name,"passed":bool(ok),**values});assert ok,name
    out={"date":"2026-10-15","audit_start":"21a875c5 plus preserved UF review/plan/AGENTS working-tree changes"}
    out["source_manifest"]=manifest(SOURCE,coverage=True)
    provenance=json.loads((SOURCE/"PROVENANCE.json").read_text())
    for row in provenance["inputs"]:
        assert sha(SOURCE/row["path"])==row["sha256"] and (SOURCE/row["path"]).stat().st_size==row["bytes"]
    out["provenance_inputs"]={"matched":len(provenance["inputs"]),"mismatches":0}
    out["stage_inputs"]={f"MT-{i:03d}":manifest(SOURCE/f"MT-{i:03d}","INPUTS.sha256") for i in range(1,6)}
    ur=ROOT/"fin_uniform_critical_review_20261014"
    out["UF_review_manifest"]=manifest(ur,coverage=True)
    previous=json.loads((ur/"verification.json").read_text())
    for rel,h in previous["input_hashes"].items():assert sha(ROOT/rel)==h,rel
    out["inherited_UF_inputs_unchanged"]=len(previous["input_hashes"])
    for local,attached in ((ROOT/"FIN_CRITICAL_METASTABILITY_NEXT_CHAIN_20261014.md",SOURCE/"source_context/MT_PLAN_USER.md"),(ur/"INTAKE_20261014.md",SOURCE/"source_context/UF_INTAKE_USER.md")):
        ck(local.name+" normalized input binding",norm(local.read_text())==norm(attached.read_text()),repo_sha256=sha(local),attached_sha256=sha(attached),byte_equal=sha(local)==sha(attached))
    with zipfile.ZipFile(SOURCE/"source_context/UF_PREDECESSOR.zip") as z:
        names=[n for n in z.namelist() if not n.endswith("/")]
        bad=[n for n in names if not (ROOT/n).exists() or hashlib.sha256(z.read(n)).hexdigest()!=sha(ROOT/n)]
    ck("UF archive file members",not bad,matched=len(names),mismatches=bad)
    out["UF_archive"]={"matched":len(names),"mismatches":bad}
    ck("archived reconstruction core",sha(SOURCE/"source_context/theory_core.py")==sha(ROOT/"fin_intrinsic_theory_campaign_20261014/theory_core.py"))
    from check_geometry import build
    geometry=build();assert json.loads(Path(__file__).with_name("GEOMETRY_CERTIFICATE.json").read_text())==geometry
    ck("exact geometric payments",True)
    try:build(F(1,10))
    except AssertionError:reject=True
    else:reject=False
    ck("reject oversized local cylinder",reject)
    spec=importlib.util.spec_from_file_location("versioned_seed",SEED);seed=importlib.util.module_from_spec(spec);spec.loader.exec_module(seed)
    _,_,lam,X,A=seed.build_rank7();gc=12/lam[6];chi=1/gc;sig=(-1.)**np.arange(12);u=np.ones(12)/12
    Q=np.eye(12)-np.ones((12,12))/12-np.outer(sig,sig)/12
    E=np.linalg.eigh(Q)[1][:,-10:];H=E.T@(12*np.eye(12)-gc*A)@E
    ck("full ten-mode fast Hessian strictly positive",np.linalg.eigvalsh(H)[0]>.22)
    ck("critical parity and projection",np.linalg.norm(A@sig-lam[6]*sig)<1e-13 and np.linalg.norm(Q@sig)<1e-13)
    # Correct anisotropic type-lattice volume, from original integer generators.
    lattice=(np.eye(12)[:,:11]-np.eye(12)[:,[11]])
    scale_map=np.vstack([sig,E.T])@lattice
    ck("scaled lattice determinant constant",abs(abs(np.linalg.det(scale_map))-12)<2e-12)
    def FE(p,g):return float(np.sum(p*np.log(12*p))-g/2*p@A@p)
    def grad(p,g):return np.log(12*p)+1-g*A@p
    for m in (-.009,0,.009):
        pm=u+m*sig/12
        ck(f"transverse gradient zero on parity line {m}",np.linalg.norm(E.T@grad(pm,gc))<2e-13)
        ck(f"positive fast Hessian on parity line {m}",np.linalg.eigvalsh(E.T@(np.diag(1/pm)-gc*A)@E)[0]>.08)
    # Central path balancing fixture: true single-copy moves, no tiny-h double certificate.
    N=120;n=np.full(12,10);n[0]+=3;n[2]-=2;n[1]-=1;n[4]+=1;n[5]-=1
    initial=n.copy();path=[n.copy()]
    def balance(parity):
        group=[j for j in range(12) if j%2==parity]
        while max(n[group])-min(n[group])>=2:
            i=max(group,key=lambda k:(n[k],-k));j=min(group,key=lambda k:(n[k],k))
            n[i]-=1;n[j]+=1;path.append(n.copy())
    for parity in (0,1):balance(parity)
    while n[::2].sum()!=N//2:
        old=0 if n[::2].sum()>N//2 else 1;new=1-old
        i=max(range(old,12,2),key=lambda k:n[k]);j=min(range(new,12,2),key=lambda k:n[k])
        n[i]-=1;n[j]+=1;path.append(n.copy())
        balance(0);balance(1)
    ck("central balancing uses only legal count edges",all(np.sum(abs(b-a))==2 and min(b)>=0 and b.sum()==N for a,b in zip(path,path[1:])))
    ck("central balancing reaches uniform type",np.array_equal(n,np.full(12,10)))
    ck("central parity mass never moves away from zero",all(abs(sig@b)<=abs(sig@a) for a,b in zip(path,path[1:])))
    out["central_path_fixture"]={"N":N,"initial":initial.tolist(),"steps":len(path)-1,"scope":"combinatorial algorithm check only; not a tiny-domain all-N energy proof"}
    # Actual total rates, instead of the producer's tautological sum n_i*2<=2N.
    count=np.full(12,800);count[0]+=2;count[1]-=2;NN=int(count.sum());rates={}
    for rule in ("HB","MET","BAR"):
        R=np.zeros((12,12))
        for i in range(12):
            deleted=count.copy();deleted[i]-=1;h=gc*A@deleted/NN
            probs=softmax(h) if rule=="HB" else np.exp(np.minimum(0,h-h[i]))/12 if rule=="MET" else expit(h-h[i])/6
            R[i]=count[i]*probs;R[i,i]=0
        total=float(R.sum());ck(rule+" actual total rate cap",0<total<=2*NN)
        moved=count.copy();moved[0]-=1;moved[1]+=1
        deleted=moved.copy();deleted[1]-=1;hh=gc*A@deleted/NN
        back=softmax(hh)[0] if rule=="HB" else math.exp(min(0,hh[0]-hh[1]))/12 if rule=="MET" else expit(hh[0]-hh[1])/6
        ratio=math.log(count[0]/(count[1]+1))+gc/(2*NN)*float(moved@A@moved-count@A@count)
        ck(rule+" direct finite Gibbs edge ratio",abs(math.log(R[0,1]/(moved[1]*back))-ratio)<1e-12)
        rates[rule]=total
    out["actual_rate_caps"]=rates
    out["geometry"]=geometry
    out["checks"]=checks
    inputs=[SEED,ROOT/"FIN_CRITICAL_METASTABILITY_NEXT_CHAIN_20261014.md",ur/"MANIFEST.sha256",SOURCE/"MANIFEST.sha256",Path(__file__),Path(__file__).with_name("check_geometry.py")]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["preexisting_AGENTS_sha256"]="00a5c189f95608a468c4fdcf89192f3d87bc2e0139e04b0029db550be5a2cf96"
    out["primary_reference"]={"doi":"10.1007/s00440-023-01257-x","publisher_capture_sha256":sha(Path("/tmp/fin-mt-canonical-primary.md")) if Path("/tmp/fin-mt-canonical-primary.md").exists() else None,
                              "scope":"Proposition4.7 inspected; its walk setup has uniform stationary measure; Gibbs-weighted inequality proved directly in review/source"}
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix="fin-mt-intake-"));copy=temp/"campaign";shutil.copytree(SOURCE,copy,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1");logs={}
        for script in ("package_check.py","mt_analysis.py","independent_verify.py","clean_replay.py"):
            ran=subprocess.run(["python3",script],cwd=copy,env=env,capture_output=True,text=True,timeout=120)
            logs[script]={"exit_code":ran.returncode,"stdout":ran.stdout,"stderr":ran.stderr};assert ran.returncode==0,(script,ran.stderr)
        numeric=json.loads((copy/"NUMERICAL_DIAGNOSTICS.json").read_text());ind=json.loads((copy/"INDEPENDENT_VERIFICATION.json").read_text())
        out["replay"]={"directory":str(copy),"producer_flags":[numeric["checks_passed"],numeric["checks_total"]],"independent_flags":[ind["passed"],ind["num_checks"]],"logs":logs,
                       "byte_identical_to_incoming":{rel:sha(copy/rel)==sha(SOURCE/rel) for rel in ("NUMERICAL_DIAGNOSTICS.json","INDEPENDENT_VERIFICATION.json")}}
    out["source_manifest_after"]=manifest(SOURCE,coverage=True)
    out["status"]="PASS_SCOPED_MT_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"source":out["source_manifest"],"inputs":out["provenance_inputs"],"stage_inputs":sum(v["matched"] for v in out["stage_inputs"].values()),"UF_archive":out["UF_archive"],"checks":len(checks),"producer":out.get("replay",{}).get("producer_flags"),"independent":out.get("replay",{}).get("independent_flags")},indent=2))

if __name__=="__main__":main()
