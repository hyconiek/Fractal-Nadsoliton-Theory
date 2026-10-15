#!/usr/bin/env python3
"""RI integrity/replay and independent small-model checks, not theorem proofs."""
import argparse
from collections import defaultdict
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
from scipy import sparse
from scipy.linalg import expm
from scipy.special import softmax

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT/"fin_record_information_campaign_20261015"
SEED = ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder, name="MANIFEST.sha256", coverage=False):
    entries = []
    for line in (folder/name).read_text().splitlines():
        h, rel = line.split(None, 1)
        assert sha(folder/rel) == h, rel
        assert rel not in entries, rel
        entries.append(rel)
    out = {"matched": len(entries), "mismatches": 0}
    if coverage:
        actual = {str(p.relative_to(folder)) for p in folder.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts}
        out["uncovered"] = sorted(actual-set(entries)-{name})
        assert not out["uncovered"]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    checks = []

    def ck(name, ok, **detail):
        checks.append({"name": name, "passed": bool(ok), **detail})
        assert ok, name

    out = {"date": "2026-10-15", "baseline": "21a875c5 plus existing UF/MT/OBS worktree additions"}
    out["source_manifest"] = manifest(SOURCE, coverage=True)
    out["stage_inputs"] = {f"RI-{i:03d}": manifest(SOURCE/f"RI-{i:03d}", "INPUTS.sha256") for i in range(1, 6)}
    obs = ROOT/"fin_intrinsic_observation_review_20261015"
    out["OBS_review_manifest"] = manifest(obs, coverage=True)
    old = json.loads((obs/"verification.json").read_text())
    for rel, h in old["input_hashes"].items():
        assert sha(ROOT/rel) == h, rel
    out["OBS_saved_inputs_unchanged"] = len(old["input_hashes"])
    declared = json.loads((SOURCE/"SOURCE_HASHES.json").read_text())["inputs"]
    binding = {}
    for i, name in ((1, "LABELLED_TRANSFER"), (2, "SNAPSHOT_CHANNEL_THEOREM"), (3, "RAW_PATH_NO_GO"), (4, "INTEGRATED_RECORD_THEOREM")):
        p = ROOT/f"fin_intrinsic_observation_campaign_20261015/OBS-C0{i}_{name}.md"
        d = declared[f"OBS_C0{i}_PROOF_DEPENDENCY"]
        binding[str(p.relative_to(ROOT))] = {"sha256_match": sha(p) == d["sha256"], "size_match": p.stat().st_size == d["size_bytes"]}
    ck("four declared OBS proof files locally bound", all(all(v.values()) for v in binding.values()))
    out["declared_predecessor_bindings"] = binding
    out["unavailable_historical_inputs"] = ["executor AGENTS snapshot", "uploaded RI plan original", "OBS ZIP wrapper as an original mounted input"]
    out["custody_scope"] = "actual local reviews/plan are bound now; no remote access, old snapshot equality or historical timing attested"

    table_file = json.loads((SOURCE/"JOINT_TABLE_OR_BOUNDS.json").read_text())
    ck("fixed nine-record three-target contract", len(table_file["cells"]) == 27 and json.loads((SOURCE/"FIXED_EXPERIMENT.json").read_text())["record"]["cutoffs"] == [-.5, .5])
    grids = {}
    for label in ("COARSE", "FINE", "REFINED"):
        d = json.loads((SOURCE/f"GRID_DIAGNOSTIC_{label}.json").read_text())
        t = np.asarray(d["table"])
        risk = 1-t.max(axis=2).sum()
        ck(label+" stored grid finite normalization and risk arithmetic", t.shape == (3,3,3) and np.min(t) >= 0 and abs(t.sum()-1)<2e-13 and abs(risk-d["risk"])<2e-14)
        ck(label+" zero-target limiting-grid convention", np.max(abs(t[:,:,1]))==0)
        grids[label] = {"risk": d["risk"], "unnormalized_total": d["unnormalized_total"],
                        "sign_reversal_residual":float(np.max(abs(t[:,:,0]-t[::-1,::-1,2])))}
    refined = np.asarray(json.loads((SOURCE/"GRID_DIAGNOSTIC_REFINED.json").read_text())["table"])
    ck("27 attached numerical cells are the REFINED grid, not certified bounds", np.max(abs(np.asarray([r["grid_approximation_only"] for r in table_file["cells"]]).reshape(3,3,3)-refined))<1e-15)
    out["stored_grid_diagnostics"] = grids
    # Exact bounded-decision identities on a separate rational toy table.
    p = [[F(3,10), F(1,10)], [F(1,10), F(3,10)], [F(1,10), F(1,10)]]
    risk = 1-sum(max(row) for row in p)
    difference = sum(abs(a-b) for a,b in p)
    ck("finite balanced-binary Bayes advantage identity", risk == F(1,2)-difference/2 and risk < F(1,2))

    spec = importlib.util.spec_from_file_location("versioned_seed", SEED)
    seed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(seed)
    _, _, lam, _, A = seed.build_rank7()
    g = 12/lam[6]
    q = 12

    def labelled_generator(N):
        states = list(itertools.product(range(q), repeat=N))
        index = {x:i for i,x in enumerate(states)}
        rows, cols, vals = [], [], []
        for ix, x in enumerate(states):
            total = 0.
            for a, old_label in enumerate(x):
                others = [x[b] for b in range(N) if b != a]
                h = g/N*A[:,others].sum(axis=1)
                prob = softmax(h)
                for j in range(q):
                    if j == old_label:
                        continue
                    dest = list(x)
                    dest[a] = j
                    rows.append(ix); cols.append(index[tuple(dest)]); vals.append(prob[j])
                    total += prob[j]
            rows.append(ix); cols.append(ix); vals.append(-total)
        return states, sparse.csr_matrix((vals,(rows,cols)), shape=(len(states),len(states)))

    states, L = labelled_generator(3)
    keys = [(tuple(x.count(i) for i in range(q)), x[0]) for x in states]
    aug = sorted(set(keys)); idx = {x:i for i,x in enumerate(aug)}
    lift = sparse.csr_matrix((np.ones(len(states)), (range(len(states)), [idx[x] for x in keys])), shape=(len(states),len(aug)))
    rows, cols, vals = [], [], []
    for ix,(n,tag) in enumerate(aug):
        total = 0.
        for i in range(q):
            if not n[i]:
                continue
            deleted = np.array(n); deleted[i] -= 1
            prob = softmax(g/3*(A@deleted))
            for j in range(q):
                if j == i:
                    continue
                nn = list(n); nn[i] -= 1; nn[j] += 1
                for new_tag, mult in ((j, int(i==tag)), (tag, n[i]-int(i==tag))):
                    if mult:
                        rows.append(ix); cols.append(idx[(tuple(nn),new_tag)]); vals.append(mult*prob[j]); total += mult*prob[j]
        rows.append(ix); cols.append(ix); vals.append(-total)
    B = sparse.csr_matrix((vals,(rows,cols)), shape=(len(aug),len(aug)))
    residual = L@lift-lift@B
    err = float(np.max(abs(residual.data))) if residual.nnz else 0.
    ck("actual q12 N3 tag-count generator intertwining", err < 3e-14, labelled_states=len(states), quotient_states=len(aug), residual=err)
    ck("actual N3 quotient rate and row sums", np.max(abs(B@np.ones(len(aug))))<3e-14 and B.diagonal().min()>=-3.0000000000001)
    out["generator_fixture"] = {"N":3,"labelled_states":len(states),"tag_count_states":len(aug),"max_intertwining_residual":err,"scope":"finite algebra diagnostic, not asymptotic critical data"}

    # Exact path-likelihood block convention: full diagonals retain killing.
    xs, L2 = labelled_generator(2)
    Q = L2.toarray()
    a = [i for i,x in enumerate(xs) if x[0]%2==0]
    b = [i for i,x in enumerate(xs) if x[0]%2==1]
    K0, K1, J = Q[np.ix_(a,a)], Q[np.ix_(b,b)], Q[np.ix_(a,b)]
    T = .4
    prior = np.ones(len(a))/len(xs)  # unnormalized initial observed-plus branch
    p0 = float(prior@expm(T*K0)@np.ones(len(a)))
    upper = np.block([[K0,J],[np.zeros_like(J),K1]])
    p1 = float(np.r_[prior,np.zeros(len(b))]@expm(T*upper)@np.r_[np.zeros(len(a)),np.ones(len(b))])
    rate = 2.
    P = np.eye(len(xs))+Q/rate
    v0 = np.zeros(len(xs)); v0[a] = prior
    v1 = np.zeros(len(xs))
    invisible = np.equal.outer(np.array([x[0]%2 for x in xs]),np.array([x[0]%2 for x in xs]))
    P0 = P*invisible; P1 = P*(~invisible)
    weight = math.exp(-rate*T)
    sums = np.array([weight*v0.sum(), 0.])
    for ell in range(1,35):
        v0, v1 = v0@P0, v1@P0+v0@P1
        weight *= rate*T/ell
        sums += weight*np.array([v0.sum(),v1.sum()])
    ck("killed blocks match independent uniformized visible-count laws", max(abs(sums-[p0,p1])) < 2e-14 and abs(prior.sum()-.5)<1e-15)
    ck("conditional prior needs explicit initial-parity factor", abs((2*prior)@expm(T*K0)@np.ones(len(a))-2*p0)<1e-14 and p0>0)
    ck("resetting killed diagonals incorrectly removes no-switch loss", np.max(abs(K0@np.ones(len(a))))>.1 and p0<prior.sum())
    out["likelihood_fixture"] = {"N":2,"initial_branch_mass":float(prior.sum()),"zero_switch_probability":p0,"one_switch_probability":p1,"scope":"finite floating semigroup diagnostic; not a full-operation error certificate"}

    # New analytic reference lemma: diagnostics of its categorical KL payment.
    payments = []
    for N in (12,120):
        n = np.full(q,N//q); tag = 0; bg = n.copy(); bg[tag] -= 1
        for old_label in (0,1,5):
            h = g/N*(A@(bg-np.eye(q,dtype=int)[old_label]))
            delta = g/N*A[:,tag]
            p = softmax(h+delta); r = softmax(h)
            kl = float(p@np.log(p/r))
            bound = float(np.ptp(delta)**2/8)
            ck(f"tag-feedback categorical KL payment N{N} old{old_label}", -1e-14<=kl<=bound+2e-15)
            payments.append({"N":N,"departure":old_label,"KL_point":kl,"range_squared_over8":bound})
    out["feedback_KL_diagnostics"] = payments

    if args.replay:
        tmp = Path(tempfile.mkdtemp(prefix="fin-ri-intake-"))/"campaign"
        shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        env = dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONDONTWRITEBYTECODE="1")
        logs = {}
        for script in ("verify_package.py","independent_checks.py"):
            ran = subprocess.run(["python3",script],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
            assert ran.returncode==0,(script,ran.stderr)
            logs[script] = {"exit_code":ran.returncode,"stdout":ran.stdout,"stderr":ran.stderr}
        drift = {}
        for label,params in (("COARSE",[".10",".10",".025","4.5","3.8"]),("FINE",[".075",".075",".0125","4.5","3.75"]),("REFINED",[".06",".06",".0125","4.5","3.78"])):
            ran = subprocess.run(["python3","deterministic_table.py",*params],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
            assert ran.returncode==0,ran.stderr
            d = json.loads(ran.stdout); original = json.loads((SOURCE/f"GRID_DIAGNOSTIC_{label}.json").read_text())
            error = float(np.max(abs(np.asarray(d["table"])-np.asarray(original["table"]))))
            ck(label+" replay numerical equivalence", error<2e-12 and abs(d["risk"]-original["risk"])<2e-12)
            drift[label] = {"max_coordinate_difference":error,"risk":d["risk"]}
        ind = json.loads((tmp/"INDEPENDENT_CHECKS.json").read_text())
        out["replay"] = {"directory":str(tmp),"logs":logs,"source_flags":[ind["pass"],ind["checks"]],"grid_drift":drift,
                         "independent_output_byte_equal":sha(tmp/"INDEPENDENT_CHECKS.json")==sha(SOURCE/"INDEPENDENT_CHECKS.json")}
    inputs = [SEED,Path(__file__),SOURCE/"MANIFEST.sha256",obs/"MANIFEST.sha256",obs/"PROOF_COMPLETIONS_AND_CORRECTIONS.md",
              ROOT/"FIN_RECORD_INFORMATION_NEXT_CHAIN_20261015.md",ROOT/"fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md"]
    out["input_hashes"] = {str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out["external_reference"] = {"url":"https://arxiv.org/html/2606.09272v1","reading_scope":"sections5.3-5.5,6.1-6.2,8.1-8.2; uniqueness assumptions7.1 contrasted, not invoked for cubic drift",
                                 "capture_sha256":sha(ROOT/".firecrawl/fin-ri-filter-primary-20261015.md") if (ROOT/".firecrawl/fin-ri-filter-primary-20261015.md").is_file() else None}
    out["source_manifest_after"] = manifest(SOURCE,coverage=True)
    out["checks"] = checks
    out["status"] = "PASS_SCOPED_RI_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"source":out["source_manifest"],"stage_entries":sum(d["matched"] for d in out["stage_inputs"].values()),"checks":len(checks),"source_flags":out.get("replay",{}).get("source_flags")},indent=2))


if __name__ == "__main__":
    main()
