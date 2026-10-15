#!/usr/bin/env python3
"""SRN audit: independent small identities and copy-only replay, not a proof assistant."""
import argparse
from collections import Counter, defaultdict
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
from scipy.special import softmax
from certify_angular_gap import build

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / "fin_scale_robust_native_observation_campaign_20261015"
SEED = ROOT / "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def manifest(folder, excluded=()):
    names = []
    for row in (folder / "MANIFEST.sha256").read_text().splitlines():
        h, name = row.split("  ", 1)
        assert name not in names and sha(folder / name) == h, name
        names.append(name)
    actual = {str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    extra = actual-set(names)-{"MANIFEST.sha256"}-set(excluded)
    assert not extra, extra
    return {"matched": len(names), "mismatches": 0, "uncovered": sorted(extra),
            "separately_bound": {n: sha(folder/n) for n in excluded}}

def integer_root(n, power):
    lo, hi = 0, 1 << ((n.bit_length()+power-1)//power)
    while hi-lo > 1:
        mid = (hi+lo)//2
        if mid**power <= n:
            lo = mid
        else:
            hi = mid
    return lo

# Features divided by sqrt(2), represented in Q(sqrt(3)).
F1 = [(F(1),F(0)), (F(0),F(1,2)), (F(1,2),F(0)), (F(0),F(0)),
      (F(-1,2),F(0)), (F(0),F(-1,2)), (F(-1),F(0)), (F(0),F(-1,2)),
      (F(-1,2),F(0)), (F(0),F(0)), (F(1,2),F(0)), (F(0),F(1,2))]
F3 = [(F((1,0,-1,0)[j % 4]), F(0)) for j in range(12)]

def add(x, y):
    return (x[0]+y[0], x[1]+y[1])

def mul(x, y):
    return (x[0]*y[0]+3*x[1]*y[1], x[0]*y[1]+x[1]*y[0])

def sign(x):
    a,b = x
    if not b:
        return (a > 0)-(a < 0)
    if not a or a*b > 0:
        return (b > 0)-(b < 0)
    d = a*a-3*b*b
    assert d != 0  # sqrt(3) irrational and a,b rational, both nonzero
    return ((a > 0)-(a < 0)) if d > 0 else ((b > 0)-(b < 0))

def angular(x, y):
    yy,xx = mul(y,y),mul(x,x)
    return int(sign((yy[0]-xx[0], yy[1]-xx[1])) > 0)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    checks = []
    def ck(name, condition, **details):
        checks.append({"name": name, "passed": bool(condition), **details})
        assert condition, name
    out = {"date": "2026-10-15", "baseline": "380e9dc0 plus preserved LI/KO/DKO/NOC worktree additions"}
    out["source_manifest"] = manifest(SOURCE, ["PACKAGE_VERIFICATION.json"])
    ledger = (SOURCE/"SOURCE_INPUTS.sha256").read_text()
    rows = [row.split("  ",1) for row in ledger.splitlines()]
    ck("source and five stage ledgers agree", (SOURCE/"INPUTS.sha256").read_text() == ledger and
       all((SOURCE/f"SRN-{i:03d}/INPUTS.sha256").read_text() == ledger for i in range(1,6)))
    matched,missing = {},[]
    for h,rel in rows:
        p = SOURCE/rel
        if p.is_file():
            assert sha(p) == h
            matched[rel] = h
        else:
            missing.append({"path": rel,"claimed_sha256":h})
    out["external_inputs"] = {"matched_originals": matched, "unavailable_originals": missing}
    out["stage_reference_counts"] = {"matched": 5*len(matched), "unavailable_original": 5*len(missing)}
    prior = {}
    for rel in ("fin_native_observation_cost_review_20261015", "fin_dynamic_kernel_observability_review_20261015",
                "fin_kernel_observability_review_20261015", "fin_critical_metastability_review_20261015"):
        folder = ROOT/rel
        item = manifest(folder)
        old = json.loads((folder/"verification.json").read_text())
        for p,h in old["input_hashes"].items():
            assert sha(ROOT/p) == h, p
        item["saved_inputs_unchanged"] = len(old["input_hashes"])
        prior[rel] = item
    out["actual_review_dependencies"] = prior
    # Verify exact gain invariance on genuine Fourier sums, including exact ties.
    total = ties = 0
    for k in (1,3,6):
        for subset in itertools.combinations(range(12),k):
            x = tuple(sum(F1[j][r] for j in subset) for r in (0,1))
            y = tuple(sum(F3[j][r] for j in subset) for r in (0,1))
            base = angular(x,y)
            for gain in (F(2),F(-3,2),F(1,7)):
                assert angular(tuple(gain*v for v in x),tuple(gain*v for v in y)) == base
            if mul(x,x) == mul(y,y):
                ties += 1
                assert base == 0
            total += 1
    ck("exact Fourier common-gain/tie fixtures in Q(sqrt3)", ties > 0, subsets=total, exact_ties=ties)
    ck("exact separate-gain and offset failures", angular((F(2),F(0)),(F(1),F(0))) == 0 and
       angular((F(2),F(0)),(F(3),F(0))) == 1 and angular((F(0),F(0)),(F(0),F(0))) == 0)
    plus = [add(F1[j],F3[j]) for j in range(12)]
    minus = [(F3[j][0]-F1[j][0],F3[j][1]-F1[j][1]) for j in range(12)]
    ck("exact pair spacings with sqrt2 removed", add(plus[0],tuple(-v for v in plus[6])) == (4,0) and
       add(minus[8],tuple(-v for v in minus[2])) == (3,0))
    conditional_cases = 0
    populations = [list(range(12)), [0,0,6,6,2,2,8,8,1,3,5,7], [0]*4+[6]*4+[2,8,1,5]]
    for nfixture,pop in enumerate(populations):
        n = Counter(pop)
        for k in (1,3,6):
            for a,b in ((0,6),(2,8)):
                groups = defaultdict(Counter)
                for subset in itertools.combinations(range(12),k):
                    draw = Counter(pop[j] for j in subset)
                    key = tuple(draw[j] for j in range(12) if j not in (a,b))
                    groups[key][draw[a]] += 1
                for key,masses in groups.items():
                    t = k-sum(key)
                    denom = math.comb(n[a]+n[b],t)
                    assert denom > 0
                    for x in range(max(0,t-n[b]), min(n[a],t)+1):
                        assert F(masses[x],sum(masses.values())) == F(math.comb(n[a],x)*math.comb(n[b],t-x),denom)
                    conditional_cases += 1
            ck(f"exact two-pair conditional support fixture{nfixture} k{k}", True)
    out["exact_conditioning_cases"] = conditional_cases
    # A deliberately broad exact finite grid tests support/max-atom bookkeeping;
    # it is not the proof of the all-population constant.
    atom_cases=0; worst=F(0)
    for A in range(1,13):
        for B in range(1,13):
            H=A+B
            if min(A,B)*8 < H:
                continue
            for t in range(1,3*H//4+1):
                mass=[F(math.comb(A,x)*math.comb(B,t-x), math.comb(H,t))
                      for x in range(max(0,t-B),min(A,t)+1)]
                assert sum(mass)==1
                worst=max(worst,t*max(mass)**2)
                atom_cases+=1
    ck("exact finite maximum-atom grid", worst < 4, cases=atom_cases, max_t_times_atom_squared=str(worst))
    spec = importlib.util.spec_from_file_location("versioned_seed",SEED)
    seed=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(seed)
    _,_,lam,_,A = seed.build_rank7()
    P3=np.array([[(1,0,-1,0)[(i-j)%4]/6 for j in range(12)] for i in range(12)])
    P4=np.array([[(1,-.5,-.5)[(i-j)%3]/6 for j in range(12)] for i in range(12)])
    f=np.sqrt(2)*np.array([[float(a)+float(b)*math.sqrt(3) for a,b in F1], [float(a) for a,b in F3]]).T
    ck("versioned null and fast feature binding", np.max(abs(A@f[:,0])) < 2e-14 and
       np.max(abs(A@f[:,1]-lam[3]*f[:,1])) < 2e-14)
    rng=np.random.default_rng(150202610)
    fixtures=[]
    for N,k in ((24,1),(24,8),(120,12),(120,60)):
        labels=rng.integers(0,12,N)
        n=np.bincount(labels,minlength=12)
        alpha=k/N
        factors=np.c_[np.full(N,1/math.sqrt(N)),((np.arange(N)<k)-alpha)/math.sqrt(k)]
        for model,MX in (("A7",A),("control",A+.01*(P3-P4))):
            drift=np.zeros(4); Bmat=np.zeros((4,4)); third=0.
            for a,old in enumerate(labels):
                probs=softmax((12/lam[6])/N*(MX@(n-np.eye(12,dtype=int)[old])))
                for j in range(12):
                    jump=np.kron(f[j]-f[old],factors[a])
                    drift += probs[j]*jump
                    Bmat += probs[j]*np.outer(jump,jump)
                    third += probs[j]*np.linalg.norm(jump)**3
            x=(f[labels].T@factors).reshape(-1)
            coeff=(np.arange(N)<k).astype(float)-alpha
            qbar=softmax((12/lam[6])/N*(MX@n))@f
            forcing=np.zeros(2)
            for a,old in enumerate(labels):
                qloo=softmax((12/lam[6])/N*(MX@(n-np.eye(12,dtype=int)[old])))@f
                forcing += coeff[a]*(qloo-qbar)/math.sqrt(k)
            ck(f"native four-dimensional slot brackets/drift N{N} k{k} {model}",
               np.max(abs(drift[[1,3]]+x[[1,3]]-forcing)) < 4e-13 and
               np.linalg.eigvalsh(Bmat).min() > -2e-13 and third > 0)
            fixtures.append({"N":N,"k":k,"model":model,"bracket":Bmat.tolist(),"third":third})
    out["generator_fixtures"] = fixtures
    # Native conditional mixed brackets need not vanish at a fixed block.
    ck("nonzero finite off-diagonal bracket retained", max(abs(row["bracket"][i][j]) for row in fixtures
       for i in range(4) for j in range(i)) > .01)
    powers={"martingale":2*F(1,32)-F(1,4)+2*F(1,64),
            "generator_N":3*F(1,32)-F(1,4)+F(1,64),
            "generator_k":3*F(1,32)-F(63,128)+F(1,64),
            "smoothing":-F(1,32)+F(1,64),"lattice":-F(63,128)+F(1,64)}
    ck("all shrinking-bandwidth powers are negative", list(powers.values()) == [F(-5,32),F(-9,64),F(-49,128),F(-1,64),F(-61,128)])
    out["power_ledger"]={k:str(v) for k,v in powers.items()}
    floors=[]
    for N in (1200,12*2**61):
        k=integer_root(N**63,64); M=integer_root(N,4)
        ck(f"exact integer-root floors N{N}", k**64<=N**63<(k+1)**64 and M**4<=N<(M+1)**4)
        floors.append({"N":N,"k":k,"M":M,"k_le_half":k<=N//2})
    out["floor_fixtures"]=floors
    cert=build()
    assert cert == json.loads((HERE/"ANGULAR_GAP_CERTIFICATE.json").read_text())
    ck("independent outward ideal slope certificate", F(cert["d_ang"][0]) > 0)
    out["ideal_certificate"]=cert
    if args.replay:
        tmp=Path(tempfile.mkdtemp(prefix="fin-srn-intake-"))/"campaign"
        shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        # Verify original manifest before the scientific writer; never run finalize_srn.
        assert manifest(tmp,["PACKAGE_VERIFICATION.json"]) == out["source_manifest"]
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1")
        run=subprocess.run(["python3","verify_srn.py"],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
        assert run.returncode==0,(run.stdout,run.stderr)
        saved=json.loads((tmp/"INDEPENDENT_VERIFICATION.json").read_text())
        original=json.loads((SOURCE/"INDEPENDENT_VERIFICATION.json").read_text())
        changes=[]
        def compare(a,b,path=""):
            assert type(a) is type(b), path
            if isinstance(a,dict):
                assert a.keys()==b.keys(), path
                for key in a:
                    compare(a[key],b[key],path+"/"+key)
            elif isinstance(a,list):
                assert len(a)==len(b), path
                for i,(x,y) in enumerate(zip(a,b)):
                    compare(x,y,path+"/"+str(i))
            elif a!=b:
                assert isinstance(a,(float,int)) and not isinstance(a,bool), (path,a,b)
                changes.append({"path":path,"absolute_drift":abs(a-b)})
        compare(original,saved)
        out["replay"]={"directory":str(tmp),"exit_code":run.returncode,"stdout":run.stdout,"stderr":run.stderr,
                       "source_checks":saved["checks_total"],
                       "numerical_output_byte_equal":sha(tmp/"INDEPENDENT_VERIFICATION.json")==sha(SOURCE/"INDEPENDENT_VERIFICATION.json"),
                       "numerical_changes":changes,"maximum_absolute_drift":max((x["absolute_drift"] for x in changes),default=0.)}
    inputs=[Path(__file__), HERE/"certify_angular_gap.py", SEED, SOURCE/"MANIFEST.sha256",
            ROOT/"FIN_SCALE_ROBUST_NATIVE_OBSERVATION_NEXT_CHAIN_20261015.md"]
    inputs += [ROOT/rel/"MANIFEST.sha256" for rel in prior]
    inputs += [ROOT/p for p in (
        "fin_native_observation_cost_review_20261015/PROOF_COMPLETIONS.md",
        "fin_dynamic_kernel_observability_review_20261015/PROOF_COMPLETIONS.md",
        "fin_kernel_observability_review_20261015/PROOF_COMPLETIONS_AND_ARITHMETIC.md",
        "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
        "fin_parity_record_sufficiency_review_20261015/KERNEL_FAMILY_SEED.json")]
    out["input_hashes"]={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["source_manifest_after"]=manifest(SOURCE,["PACKAGE_VERIFICATION.json"])
    assert out["source_manifest_after"]==out["source_manifest"]
    out["checks"]=checks
    out["status"]="PASS_SCOPED_SRN_INTAKE_CHECKS"
    out["not_tested"]=["all-N anti-concentration or generator estimates by enumeration", "mechanized MT/SRN proofs", "native numerical error or N0", "physical instrument"]
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":out["status"],"source":out["source_manifest"],"source_checks":out.get("replay",{}).get("source_checks"),
                      "checks":len(checks),"conditioning_cases":conditional_cases,"missing_originals":len(missing),"ideal":cert},indent=2))

if __name__=="__main__":
    main()
