#!/usr/bin/env python3
"""Integrity, exact finite combinatorics and independent generator diagnostics.

No finite fixture proves an asymptotic process/observation theorem.
Historical writers run only in a temporary, cache-free copy.
"""
import argparse
from fractions import Fraction as F
import hashlib
import html
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import zipfile

import mpmath as mp
import numpy as np
from scipy.special import softmax
from observation_kernels import hypergeometric, two_blocks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / "fin_intrinsic_observation_campaign_20261015"
SEED = ROOT / "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder, name="MANIFEST.sha256", coverage=False):
    entries = []
    for line in (folder / name).read_text().splitlines():
        h, rel = line.split(None, 1)
        rel = rel.lstrip("*")
        assert sha(folder / rel) == h, rel
        entries.append(rel)
    result = {"matched": len(entries), "mismatches": 0}
    if coverage:
        actual = {str(p.relative_to(folder)) for p in folder.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
        extra = actual - set(entries) - {name}
        assert not extra, extra
        result["uncovered"] = sorted(extra)
    return result


def normalized(s):
    text = re.sub(r"\\([\\`*_{}\[\]()#+.!><=:|~\-])", r"\1", html.unescape(s)).replace("**", "")
    # One observed single-emphasis formatting artifact, not general deletion
    # of multiplication symbols or scientific differences.
    text = text.replace("sqrt(N)<<k*<<N*", "sqrt(N)<<k<<N")
    return " ".join(text.split())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    checks = []

    def ck(name, ok, **detail):
        checks.append({"name": name, "passed": bool(ok), **detail})
        assert ok, name

    out = {"date": "2026-10-15", "baseline": "21a875c5 plus preserved UF/MT working-tree additions"}
    out["source_manifest"] = manifest(SOURCE, coverage=True)
    out["source_inputs"] = manifest(SOURCE, "SOURCE_INPUTS.sha256")
    out["stage_inputs"] = {f"OBS-C0{i}": manifest(SOURCE / f"OBS-C0{i}", "INPUTS.sha256")
                           for i in range(1, 6)}
    mr = ROOT / "fin_critical_metastability_review_20261015"
    out["MT_review_manifest"] = manifest(mr, coverage=True)
    predecessor = json.loads((mr / "verification.json").read_text())
    for rel, h in predecessor["input_hashes"].items():
        assert sha(ROOT / rel) == h, rel
    out["MT_scientific_inputs_unchanged"] = len(predecessor["input_hashes"])
    for local, attached in ((ROOT / "FIN_INTRINSIC_OBSERVATION_NEXT_CHAIN_20261015.md", SOURCE / "source_context/OBS_PLAN_USER.md"),
                            (ROOT / "FIN_CRITICAL_METASTABILITY_NEXT_CHAIN_20261014.md", SOURCE / "source_context/MT_INTAKE_USER.md")):
        ck(local.name + " normalized binding", normalized(local.read_text()) == normalized(attached.read_text()),
           repo_sha256=sha(local), attached_sha256=sha(attached), byte_equal=sha(local) == sha(attached))
    out["provenance_discrepancy"] = {
        "file": "source_context/MT_INTAKE_USER.md", "actual_content": "format-transformed MT research plan, not the MT intake",
        "actual_MT_intake_sha256": sha(mr/"INTAKE_20261015.md"),
        "disposition": "inherited theorem paid by real local MT review and completions; no retroactive producer-access attestation"}
    ck("mislabeled attachment not falsely accepted as MT intake", normalized((SOURCE/"source_context/MT_INTAKE_USER.md").read_text()) != normalized((mr/"INTAKE_20261015.md").read_text()))
    with zipfile.ZipFile(SOURCE / "source_context/FIN_MT_FULL_HANDOFF.zip") as z:
        names = [n for n in z.namelist() if not n.endswith("/")]
        bad = [n for n in names if not (ROOT / n).is_file() or hashlib.sha256(z.read(n)).hexdigest() != sha(ROOT / n)]
    ck("MT archive members match actual predecessor", not bad, members=len(names))
    out["MT_archive"] = {"matched": len(names), "mismatches": bad}
    ck("archived numerical core binding", sha(SOURCE / "source_context/theory_core.py") ==
       sha(ROOT / "fin_intrinsic_theory_campaign_20261014/theory_core.py"))

    # Exact probabilities including impossible first-block outcomes and edge populations.
    cases = 0
    for N in (2, 6, 12):
        for M in range(N + 1):
            for k in range(1, N + 1):
                ps = [hypergeometric(N, M, k, h) for h in range(k + 1)]
                mean = F(k * M, N)
                assert sum(ps) == 1
                assert sum(h * p for h, p in enumerate(ps)) == mean
                assert sum((h - mean)**2 * p for h, p in enumerate(ps)) == F(k * M * (N - M) * (N - k), N*N*(N-1))
                cases += 1
    ck("exact hypergeometric normalization mean and variance", True, cases=cases)
    joint_cases = 0
    for M in range(13):
        for k1, k2 in ((1, 1), (3, 4), (6, 6), (12, 0)):
            rows = [[two_blocks(12, M, k1, k2, h1, h2) for h2 in range(k2 + 1)] for h1 in range(k1 + 1)]
            assert sum(map(sum, rows)) == 1
            assert all(sum(row) == hypergeometric(12, M, k1, h1) for h1, row in enumerate(rows))
            for h1, row in enumerate(rows):
                for h2, p in enumerate(row):
                    # Independent allocation formula across three fixed groups.
                    rem = M - h1 - h2
                    expected = (F(math.comb(k1, h1)*math.comb(k2, h2)*math.comb(12-k1-k2, rem), math.comb(12, M))
                                if 0 <= rem <= 12-k1-k2 else F(0))
                    assert p == expected
            if k1*k2:
                cov = sum((F(h1)-F(k1*M, 12))*(F(h2)-F(k2*M, 12))*rows[h1][h2]
                          for h1 in range(k1+1) for h2 in range(k2+1))
                assert cov == -F(k1*k2*M*(12-M), 12*12*11)
            joint_cases += 1
    ck("exact two-block allocation normalization support covariance", True, cases=joint_cases)
    invalid_rejected = []
    for values in ((12, 6, True, 0), (12, 6, 1.5, 0), (12, 13, 2, 0), (0, 0, 0, 0)):
        try:
            hypergeometric(*values)
        except (TypeError, ValueError):
            invalid_rejected.append(values)
    ck("new exact helper fails closed on invalid contracts", len(invalid_rejected) == 4)

    spec = importlib.util.spec_from_file_location("versioned_seed", SEED)
    seed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(seed)
    _, _, lam, X, A = seed.build_rank7()
    gc = 12 / lam[6]
    sigma = (-1.)**np.arange(12)
    ck("versioned critical eigenvector and centered kernel", np.linalg.norm(A@sigma-lam[6]*sigma) < 1e-13 and np.linalg.norm(A@np.ones(12)) < 1e-13)
    source_out = json.loads((SOURCE / "NUMERICAL_RESULTS.json").read_text())
    ck("critical numerical input binding", abs(gc-source_out["gc_nominal"]) < 1e-13)

    def direct(labels, k, g):
        N = len(labels)
        n = np.bincount(labels, minlength=12)
        c = np.array([(1/k if a < k else 0)-1/N for a in range(N)])
        drift = np.zeros(2)
        bracket = np.zeros((2, 2))
        probs = []
        for a, old in enumerate(labels):
            deleted = n.copy()
            deleted[old] -= 1
            p = softmax(g/N*(A@deleted))
            probs.append(p)
            for j in range(12):
                jump = (sigma[j]-sigma[old])*np.array([1/N, c[a]])
                drift += p[j]*jump
                bracket += p[j]*np.outer(jump, jump)
        q = np.asarray(probs)@sigma
        h = 2*(1-sigma[labels]*q)
        closed = np.array([[h.sum()/N**2, (h*c).sum()/N], [(h*c).sum()/N, (h*c*c).sum()]])
        m, mb = np.mean(sigma[labels]), np.mean(sigma[labels[:k]])
        e = q[:k].mean()-q.mean()
        return drift, bracket, closed, m, mb, e

    rng = np.random.default_rng(20261015)
    metrics = []
    max_q_bound = 2*abs(gc)*float(np.max(np.abs(A)))  # |q_i-q_j| <= 2 g max|A|/N.
    for N, k in ((12, 1), (24, 2), (120, 3), (120, 12)):
        labels = rng.integers(0, 12, N)
        drift, bracket, closed, m, mb, e = direct(labels, k, gc)
        ck(f"actual tagged drift and all brackets N{N} k{k}", abs(drift[1]+mb-m-e) < 3e-14 and np.max(abs(bracket-closed)) < 3e-14)
        ck(f"uniform LOO forcing and cross-payment N{N} k{k}", abs(e)*N <= max_q_bound+1e-13 and abs(bracket[0, 1]) <= 8/N and bracket[0, 0] <= 4/N and bracket[1, 1] <= 4*(1/k-1/N)+1e-14)
        permuted = labels.copy()
        rng.shuffle(permuted)
        d2, b2, *_ = direct(permuted, k, gc)
        ck(f"count drift bracket slot permutation N{N} k{k}", abs(drift[0]-d2[0]) < 3e-14 and abs(bracket[0, 0]-b2[0, 0]) < 3e-14)
        metrics.append({"N": N, "k": k, "drift": drift.tolist(), "bracket": bracket.tolist(), "N_e": N*float(e)})
    out["actual_generator_fixtures"] = metrics
    for N in (120, 1200):
        labels = np.repeat(np.arange(12), N//12)
        _, bracket, _, _, _, _ = direct(labels, 1, gc)
        ck(f"uniform finite-jump noise diagnostic N{N}", abs(N*bracket[0, 0]-2) < 4/N and abs(bracket[1, 1]-2*(1-1/N)) < 4/N)

    # A real labelled/count-invisible swap changes a tag observable, without new primary kinetics.
    labels = np.array([0, 1, 2])
    swapped = labels[[1, 0, 2]]
    ck("outside-primary swap diagnostic", np.array_equal(np.bincount(labels, minlength=12), np.bincount(swapped, minlength=12)) and sigma[labels[0]] != sigma[swapped[0]])
    # Orbit TV isometry for a tiny exact finite model; not the MT large-N theorem.
    states = list(itertools.product(range(2), repeat=3))
    p = [F(1, 10), F(2, 10), F(3, 10), F(4, 10)]
    q = p[::-1]
    lifted = sum(abs(p[sum(x)]-q[sum(x)])/math.comb(3, sum(x)) for x in states)/2
    ck("exact disjoint orbit TV isometry fixture", lifted == sum(abs(a-b) for a, b in zip(p, q))/2)

    mp.mp.dps = 40
    norm = 2*mp.quad(lambda z: mp.exp(-z**4/12), [0, 1, 3, mp.inf])
    risks = {}
    for alpha in ("0.25", "1.0", "4.0"):
        v = 2*mp.quad(lambda z: mp.exp(-z**4/12)*mp.erfc(mp.sqrt(mp.mpf(alpha)/2)*z)/2, [0, 1, 3, mp.inf])/norm
        risks[alpha] = str(v)
    ck("independent mpmath quadrature of bounded sign risks", all(abs(float(risks[a])-source_out["bounded_sign_risk_nominal"]["alpha"][a]) < 2e-12 for a in risks))
    ck("OU limiting increment variance scalar", abs(source_out["critical_private_ou_increment_variance_alpha1_delta_modeltime1"] - 2*(1-math.exp(-1))) < 1e-14)
    out["risk_point_diagnostics_not_certificates"] = risks

    if args.replay:
        tmp = Path(tempfile.mkdtemp(prefix="fin-obs-intake-")) / "campaign"
        shutil.copytree(SOURCE, tmp, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
        logs = {}
        for script in ("verify_package.py", "run_observation.py", "independent_verify.py"):
            ran = subprocess.run(["python3", script], cwd=tmp, env=env, capture_output=True, text=True, timeout=60)
            assert ran.returncode == 0, (script, ran.stderr)
            logs[script] = {"exit_code": ran.returncode, "stdout": ran.stdout, "stderr": ran.stderr}
        ind = json.loads((tmp / "INDEPENDENT_VERIFICATION.json").read_text())
        out["replay"] = {"directory": str(tmp), "logs": logs, "source_flags": [ind["passed"], ind["total"]],
                         "byte_equal": {n: sha(tmp/n) == sha(SOURCE/n) for n in ("NUMERICAL_RESULTS.json", "INDEPENDENT_VERIFICATION.json")}}
        probe = subprocess.run(["python3", "-c", "from observation_core import joint_hypergeom_pmf; print(joint_hypergeom_pmf(12,2,4,3,3,0))"],
                               cwd=tmp, env=env, capture_output=True, text=True, timeout=15)
        ck("source impossible-outcome bug reproduced and isolated", probe.returncode != 0 and "ValueError" in probe.stderr and two_blocks(12, 2, 4, 3, 3, 0) == 0)
        out["source_helper_bug"] = {"input": [12, 2, 4, 3, 3, 0], "expected_exact_probability": "0", "historical_behavior": "raises ValueError instead of zero", "stderr": probe.stderr,
                                    "disposition": "source preserved; reviewed separate exact helper handles zero-support branch"}
    inputs = [SEED, Path(__file__), HERE/"observation_kernels.py", SOURCE/"MANIFEST.sha256", mr/"MANIFEST.sha256",
              mr/"PROOF_COMPLETIONS.md", ROOT/"FIN_INTRINSIC_OBSERVATION_NEXT_CHAIN_20261015.md",
              ROOT/"fin_uniform_critical_review_20261014/INTAKE_20261014.md",
              ROOT/"fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md",
              ROOT/"fin_stochastic_fold_closure_review_20261014/STOPPED_LIMIT_COMPLETION.md"]
    out["input_hashes"] = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    capture = ROOT/".firecrawl/fin-obs-whitt-20261015.md"
    out["primary_reference"] = {"url": "https://arxiv.org/html/0712.1929v2", "capture_sha256": sha(capture) if capture.is_file() else None,
                                "reading_scope": "Theorem2.1(ii), Theorem5.3 and proof; characterization proof in5.1", "Hajek": "not relied on as unread black-box theorem; direct bounded-draw martingale proof in review"}
    out["source_manifest_after"] = manifest(SOURCE, coverage=True)
    out["checks"] = checks
    out["status"] = "PASS_SCOPED_OBS_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"manifest": out["source_manifest"], "source_inputs": out["source_inputs"], "stage_entries": sum(v["matched"] for v in out["stage_inputs"].values()),
                      "archive": out["MT_archive"], "checks": len(checks), "source_flags": out.get("replay", {}).get("source_flags"), "prototype_bug": out.get("source_helper_bug", {}).get("historical_behavior")}, indent=2))


if __name__ == "__main__":
    main()
