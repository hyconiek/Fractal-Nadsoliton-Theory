#!/usr/bin/env python3
"""Provenance and bounded independent checks; not a stochastic-limit proof.

Archived producer scripts run only in a temporary copy. Exact statements
are explained in the review; floating checks do not certify them by sampling.
"""
from pathlib import Path
import argparse
import hashlib
import html
import importlib.util
import json
import re
import shutil
import subprocess
import tempfile
import os
import sys
from fractions import Fraction
import numpy as np
from scipy.optimize import root
from scipy.special import softmax
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "fin_intrinsic_theory_campaign_20261014"
FOLD = ROOT / "fin_rank7_followup/certificates/R7P-031_simple_fold.json"
SEED = ROOT / "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(folder, name="MANIFEST.sha256", complete=False):
    entries = []
    optional_missing = []
    for line in (folder / name).read_text().splitlines():
        digest, rel = line.split(None, 1)
        rel = str(Path(rel.lstrip("*")))
        if not (folder / rel).exists() and "__pycache__" in Path(rel).parts:
            optional_missing.append(rel)
            continue
        assert sha(folder / rel) == digest, rel
        entries.append(rel)
    actual = {str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file()}
    uncovered = sorted(actual - set(entries) - {name})
    scientific_uncovered = [p for p in uncovered if "__pycache__" not in Path(p).parts]
    if complete:
        assert not scientific_uncovered, scientific_uncovered
    return {"matched": len(entries), "mismatches": 0,
            "uncovered": scientific_uncovered if complete else "NOT_A_PACKAGE_COVERAGE_MANIFEST",
            "optional_missing_cache": optional_missing,
            "ignored_cache_files": [p for p in uncovered if p not in scientific_uncovered] if complete else []}


def normalized(text):
    text = html.unescape(text)
    text = re.sub(r"\\([\\`*_{}\[\]()#+.!>|\-])", r"\1", text)
    return " ".join(text.replace("**", "").split())


def interval_endpoints(value, places=14):
    # Convert binary interval endpoints exactly, then round decimal outwards.
    def exact_endpoint(endpoint):
        sign, man, exp, _ = endpoint._mpi_[0]
        number = Fraction((-1 if sign else 1) * man)
        return number * (2 ** exp if exp >= 0 else Fraction(1, 2 ** (-exp)))
    scale = 10 ** places
    lo, hi = exact_endpoint(value.a), exact_endpoint(value.b)
    lower = (lo * scale).numerator // (lo * scale).denominator
    upper = -((-hi * scale).numerator // (-hi * scale).denominator)
    def decimal(k):
        sign = "-" if k < 0 else ""
        k = abs(k)
        return f"{sign}{k // scale}.{k % scale:0{places}d}"
    return [decimal(lower), decimal(upper)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    checks = []
    def ck(name, passed, **data):
        checks.append({"name": name, "passed": bool(passed), **data})
        assert passed, name

    out = {"baseline": "ec8baa15", "date": "2026-10-14"}
    out["source_manifest"] = manifest(SOURCE, complete=True)
    out["source_input_manifest"] = manifest(SOURCE, "SOURCE_INPUTS.sha256")
    out["stage_input_manifests"] = {
        f"TH-{t:03d}": manifest(SOURCE, f"TH-{t:03d}/INPUTS.sha256")
        for t in range(1, 6)
    }
    out["synthesis_manifest"] = manifest(ROOT / "fin_theory_synthesis_20261014", complete=True)
    plan = ROOT / "FIN_INTRINSIC_THEORY_NEXT_CAMPAIGN_20261014.md"
    attached = SOURCE / "source_context/FIN_INTRINSIC_THEORY_NEXT_CAMPAIGN_20261014_USER.md"
    ck("attached plan normalized text", normalized(plan.read_text()) == normalized(attached.read_text()))
    out["plan_binding"] = {"repo_sha256": sha(plan), "attached_sha256": sha(attached),
                           "byte_identical": sha(plan) == sha(attached), "normalized_equal": True}
    fold = json.loads(FOLD.read_text())
    ck("fold certificate actual hash", sha(FOLD) == "5537174281996415a20ea3dffce23c33258c51396b3850aedf4168d9280e352e")
    ck("saved Krawczyk strict inclusion", all(
        Fraction(box[0]) < Fraction(image[0]) < Fraction(image[1]) < Fraction(box[1])
        for box, image in zip(fold["augmented_root_box"], fold["krawczyk_image"])))
    # Small replay of the existing interval root, not a new root census.
    r7p = ROOT / "fin_rank7_followup"
    sys.path.insert(0, str(r7p))
    fold_spec = importlib.util.spec_from_file_location("review_fold_certificate", r7p/"src/fold_certificate.py")
    fold_module = importlib.util.module_from_spec(fold_spec)
    fold_spec.loader.exec_module(fold_module)
    _, reproduced = fold_module.certify()
    ck("replayed fold Krawczyk and inertia", reproduced["strict_inclusion"]
       and reproduced["simple_fold_transversality"] and reproduced["H4_principal_positive"]
       and reproduced["Hsin_positive_definite"] and reproduced["H7_inertia_at_root"] == [0,1,6])
    out["old_fold_interval_replay"] = {"backend":"mpmath.iv; existing implementation, not a new independent proof",
                                       "strict_inclusion":reproduced["strict_inclusion"],
                                       "inertia":reproduced["H7_inertia_at_root"]}

    spec = importlib.util.spec_from_file_location("versioned_fin_seed", SEED)
    seed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(seed)
    _, _, lam, X, A = seed.build_rank7()
    saved = json.loads((SOURCE / "NUMERICAL_REPLAY.json").read_text())
    ck("versioned spectrum vs source", np.max(abs(lam - saved["laplacian_eigenvalues_0_6"])) < 2e-14)
    # Independent fold equation: smallest eigenvalue, not producer determinant.
    E = X[:, [0, 2, 4, 6]]
    def components(th, g):
        p = softmax(E @ th)
        M = E.T @ (np.diag(p) - np.outer(p, p)) @ E
        return p, th/g - E.T@p, np.eye(4)/g - M
    start = np.array([float((Fraction(a)+Fraction(b))/2)
                      for a, b in fold["augmented_root_box"][:5]])
    def equations(z):
        _, grad, H = components(z[:4], z[4])
        return np.r_[grad, np.linalg.eigvalsh(H)[0]]
    sol = root(equations, start, tol=1e-11)
    th, g = sol.x[:4], sol.x[4]
    p, grad, H4 = components(th, g)
    _, evectors = np.linalg.eigh(H4)
    v4 = evectors[:, 0]
    if v4[0] < 0:
        v4 = -v4
    th7, v = np.zeros(7), np.zeros(7)
    th7[[0,2,4,6]], v[[0,2,4,6]] = th, v4
    ck("independent fold residual", np.max(abs(equations(sol.x))) < 2e-13)
    Sigma = np.diag(p) - np.outer(p, p)
    M = X.T @ Sigma @ X
    H7 = np.eye(7)/g - M
    evals = np.linalg.eigvalsh(H7)
    ck("six retained stable directions numerical", evals[1] > .17, eigenvalues=evals.tolist())
    # Full occupation linearization; hidden drift is -identity, not a guess.
    u = np.ones(12)/12
    T = np.linalg.qr(np.eye(12)[:, :11] - np.eye(12)[:, [11]])[0]
    J = g*Sigma@A - np.eye(12)
    full_eigs = np.sort(np.linalg.eigvals(T.T@J@T).real)
    expected = np.sort(np.r_[-g*evals, -np.ones(4)])
    ck("full eleven-mode linear spectrum", np.max(abs(full_eigs-expected)) < 3e-14,
       eigenvalues=full_eigs.tolist())
    B = np.column_stack([np.ones(12), X])
    R = np.eye(12) - B @ np.linalg.solve(B.T@B, B.T)
    ck("hidden drift has exactly unit decay at equal means", np.linalg.norm(R@J@R+R) < 2e-14)
    noise = 2*g*g*M
    ck("microscopic soft noise 2g", abs(v@noise@v - 2*g) < 3e-13)
    projected_correction = -g*g*float(v@M@(X.T@p))
    ck("leave-one-out correction coefficient", abs(projected_correction + v@th7) < 2e-13)
    z = X@v; delta = g*(z[None,:]-z[:,None]); weights = p[:,None]*p[None,:]
    cumulants = {str(k): float(np.sum(weights*delta**k)) for k in range(1,5)}
    ck("fourth soft jump cumulant nonzero", cumulants["4"] > 100)
    j = np.arange(12)
    c = {k:np.cos(2*np.pi*k*j/12) for k in range(1,7)}
    s = {k:np.sin(2*np.pi*k*j/12) for k in range(1,6)}
    Lift = np.vstack([np.ones(12), X.T, c[3]*c[4], s[3]*c[4], c[3]*c[5], s[3]*c[5]])
    ck("eleven occupation coordinate lift", np.linalg.matrix_rank(Lift, tol=1e-12) == 12)
    for nn in (12,120,240):
        raw = np.arange(1,13, dtype=float)
        targets = nn*raw/raw.sum()
        counts = np.floor(targets).astype(int)
        counts[np.argsort(-(targets-counts))[:nn-int(counts.sum())]] += 1
        testp = counts/nn
        ck(f"lift reconstruction actual N{nn} type", np.max(abs(np.linalg.solve(Lift, Lift@testp)-testp)) < 2e-15)
    h = np.array([2,1,-1,-2,-1,1,2,1,-1,-2,-1,1])
    nplus, nminus = 10+h, 10-h
    ck("exact count witness equal means numerical", np.max(abs(X.T@((nplus-nminus)/120))) < 2e-14)
    # Closed g0 expression, independent of producer edge enumeration.
    def square(n):
        pp=n/120; mu=pp@X[:,4]
        return -2*mu*mu + (np.mean(X[:,4]**2)+pp@(X[:,4]**2))/120
    difference = square(nplus)-square(nminus)
    ck("finite witness square generator", abs(difference-lam[5]/7200) < 2e-16)

    # Bind mathematical exact-decimal kernel to accepted rational spectral boxes.
    mp.iv.dps = 70
    iv = mp.iv
    weights_iv = [iv.mpf(0)] + [iv.cos(iv.mpf(".18575")*d+iv.mpf(".1625"))/(1+iv.mpf(d)**iv.mpf("1.8")) for d in range(1,7)]
    old_spectra = json.loads((ROOT / "fin_rank7_followup/inputs/fin_handoff_audit/results.json").read_text())["exact"]["laplacian_intervals"]
    spectral_bounds = {}
    for k in (3,4,5,6):
        value = sum(2*weights_iv[d]*(1-iv.cos(2*iv.pi*k*d/12)) for d in range(1,6))
        value += weights_iv[6]*(1-iv.cos(iv.pi*k))
        endpoints = interval_endpoints(value, 45)
        lo, hi = map(Fraction, endpoints)
        oldlo, oldhi = map(Fraction, old_spectra[k])
        ck(f"exact-decimal spectrum inside accepted lambda{k}", oldlo <= lo <= hi <= oldhi)
        spectral_bounds[str(k)] = endpoints
    gg = iv.mpf(fold["fold_gain_interval"])
    aa = iv.mpf(fold["fold_coeff_v_dot_Fg"])
    bb = iv.mpf(fold["fold_coeff_D3Phi_vvv"])
    cc = iv.exp(iv.ln(4/bb)/3)
    scalars = {"variance_2g": 2*gg, "sigma": iv.sqrt(2*gg), "c": cc,
               "d": cc*cc/(2*gg), "kappa": -aa*cc/2}
    scalar_bounds = {k:interval_endpoints(val) for k,val in scalars.items()}
    out["independent_interval_bounds"] = {"backend": "mpmath.iv assumption; exact rational endpoint serialization",
                                          "spectra": spectral_bounds, "fold_scalars": scalar_bounds}
    out["independent_numeric"] = {"fold_gain":float(g), "fold_theta":th7.tolist(), "fold_v":v.tolist(),
                                  "full_occupation_spectrum":full_eigs.tolist(), "cumulants":cumulants,
                                  "finite_witness_difference":float(difference)}
    out["checks"] = checks
    inputs = [SEED, FOLD, plan, SOURCE/"MANIFEST.sha256", Path(__file__),
              ROOT/"fin_theory_synthesis_20261014/MANIFEST.sha256",
              ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"]
    inputs += [r7p/"src"/name for name in
               ("fold_certificate.py", "coexistence_certificate.py", "model.py", "derivatives.py")]
    inputs += [ROOT/"fin_rank7_mathphysics_next/proofs"/name for name in
               ("MP7-025_simple_fold_normal_form.md", "MP7-026_controlled_fold_scaling.md")]
    out["input_hashes"] = {str(p.relative_to(ROOT)):sha(p) for p in inputs}
    baseline_agents = subprocess.run(["git", "show", "ec8baa15:AGENTS.md"], cwd=ROOT, capture_output=True, check=True).stdout
    out["guardrail_baseline"] = {"commit":"ec8baa15", "AGENTS_sha256":hashlib.sha256(baseline_agents).hexdigest()}
    capture = Path("/tmp/fin-th-intake-kraaij.md")
    out["primary_theorem"] = {"url":"https://arxiv.org/html/1502.05185v3#S2.Thmtheorem3",
                              "scope":"full Theorem 2.3, conditions and initial LDP read manually",
                              "fresh_capture_sha256":sha(capture) if capture.exists() else None}
    if args.replay:
        temp = Path(tempfile.mkdtemp(prefix="fin-th-intake-"))
        copy = temp / "campaign"
        shutil.copytree(SOURCE, copy)
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
        logs = {}
        for script in ("theory_core.py", "fold_diagnostics.py", "independent_verify.py", "verify_package.py"):
            result = subprocess.run(["python3", script], cwd=copy, env=env, capture_output=True, text=True, timeout=120)
            logs[script] = {"exit_code":result.returncode, "stdout":result.stdout,
                            "stderr":result.stderr}
            assert result.returncode == 0, (script,result.stderr)
            if script == "theory_core.py":
                (copy/"NUMERICAL_REPLAY.json").write_text(result.stdout)
        ind = json.loads((copy/"INDEPENDENT_VERIFICATION.json").read_text())
        drifts = {}
        for rel in ("NUMERICAL_REPLAY.json", "TH-005/FOLD_GENERATOR_DIAGNOSTICS.json", "INDEPENDENT_VERIFICATION.json"):
            drifts[rel] = {"byte_identical":sha(copy/rel)==sha(SOURCE/rel)}
        original=json.loads((SOURCE/"TH-005/FOLD_GENERATOR_DIAGNOSTICS.json").read_text())
        replayed=json.loads((copy/"TH-005/FOLD_GENERATOR_DIAGNOSTICS.json").read_text())
        maxdrift=max(abs(a[k]-b[k]) for a,b in zip(original,replayed)
                     for k in a if isinstance(a[k],(int,float)))
        out["replay"] = {"directory":str(copy), "independent_checks":[ind["passed"],ind["total"]],
                         "outputs":drifts, "max_characteristic_coordinate_drift":maxdrift,
                         "logs":logs, "scope":"numerical reproducibility, not interval propagation or theorem proof"}
    out["source_manifest_after"] = manifest(SOURCE, complete=True)
    out["status"] = "PASS_SCOPED_INTAKE_CHECKS"
    args.output.write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps({"source_manifest":out["source_manifest"], "synthesis_manifest":out["synthesis_manifest"],
                      "checks":len(checks), "interval_bounds":scalar_bounds,
                      "replay_checks":out.get("replay",{}).get("independent_checks")},indent=2))


if __name__ == "__main__":
    main()
