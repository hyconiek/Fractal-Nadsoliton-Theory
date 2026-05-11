#!/usr/bin/env python3
"""P1152 strict candidate registry runner.

Runs P1151 pipeline for each candidate JSON listed in a registry file.
Exports aggregate pass/fail table for comparative strict-rigor screening.
"""
from __future__ import annotations
import json
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run_pipeline(candidate_json: Path) -> dict:
    cmd = [sys.executable, str(ROOT / "p1151_strict_selector_pipeline_runner.py"), str(candidate_json)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    failed_step = None
    summary_path = GEN / "p1151_strict_selector_pipeline_runner_summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            failed_step = summary.get("failed_step")
        except Exception:
            failed_step = None
    return {
        "candidate": str(candidate_json),
        "returncode": p.returncode,
        "stdout": p.stdout.strip(),
        "stderr": p.stderr.strip(),
        "pass": p.returncode == 0,
        "failure_stage": failed_step,
    }



def safe_region_pass(candidate_json: Path) -> bool:
    try:
        c = json.loads(candidate_json.read_text(encoding="utf-8"))
        b = json.loads((GEN / "p1172_safe_operating_region_summary.json").read_text(encoding="utf-8"))["safe_bounds"]
        return (b["omega"]["min"] <= c.get("omega_hint", 0) <= b["omega"]["max"]
                and b["phi"]["min"] <= c.get("phi_hint", 0) <= b["phi"]["max"]
                and b["sigma"]["min"] <= c.get("sigma_hint", 0) <= b["sigma"]["max"]
                and b["kappa"]["min"] <= c.get("kappa_hint", 0) <= b["kappa"]["max"])
    except Exception:
        return False

def main() -> None:
    args = sys.argv[1:]
    allow_failures = "--allow-failures" in args
    require_safe_region = "--require-safe-region" in args
    rank_by_safe_margin = "--rank-by-safe-margin" in args
    enforce_shortlist_consistency = "--enforce-shortlist-consistency" in args
    shortlist_output = None
    cleaned_args = []
    i = 0
    while i < len(args):
        token = args[i]
        if token in ("--allow-failures", "--require-safe-region", "--rank-by-safe-margin", "--enforce-shortlist-consistency"):
            i += 1
            continue
        if token == "--shortlist-output" and i + 1 < len(args):
            shortlist_output = Path(args[i + 1]).resolve()
            i += 2
            continue
        cleaned_args.append(token)
        i += 1
    args = cleaned_args
    registry = Path(args[0]).resolve() if len(args) > 0 else (GEN / "p1152_candidate_registry_example.json").resolve()
    items = json.loads(registry.read_text(encoding="utf-8"))
    candidates = [Path(x).resolve() for x in items.get("candidates", [])]

    results = []
    for c in candidates:
        if require_safe_region and (not safe_region_pass(c)):
            results.append({"candidate": str(c), "pass": False, "returncode": 9, "stdout": "", "stderr": "filtered_by_safe_region", "failure_stage": {"index": -1, "cmd": ["safe_region_gate"]}})
            continue
        results.append(run_pipeline(c))
    passed = sum(1 for r in results if r["pass"])
    failed = len(results) - passed

    ranking_triggered = False
    ranking_summary_path = GEN / "p1153_strict_candidate_quality_ranking_summary.json"
    shortlist = []
    shortlist_size = 3
    shortlist_output = shortlist_output or (GEN / "p1152_safe_margin_shortlist.json")
    pre_rank_out = {
        "packet": "P1152",
        "as_of": "2026-05-10",
        "allow_failures": allow_failures,
        "require_safe_region": require_safe_region,
        "rank_by_safe_margin": rank_by_safe_margin,
        "enforce_shortlist_consistency": enforce_shortlist_consistency,
        "ranking_triggered": False,
        "ranking_summary": None,
        "safe_margin_shortlist_size": 0,
        "safe_margin_shortlist": [],
        "registry": str(registry),
        "total": len(results),
        "passed": passed,
        "failed": failed,
        "results": results,
    }
    out_path = GEN / "p1152_strict_candidate_registry_runner_summary.json"
    out_path.write_text(json.dumps(pre_rank_out, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if rank_by_safe_margin:
        subprocess.run([sys.executable, str(ROOT / "p1153_strict_candidate_quality_ranking.py"), str(out_path)], capture_output=True, text=True)
        ranking_triggered = True
        try:
            ranking_payload = json.loads(ranking_summary_path.read_text(encoding="utf-8"))
            ranked = ranking_payload.get("ranking", [])
            ranked = [r for r in ranked if r.get("pipeline_pass") and isinstance(r.get("safe_region_margin"), (int, float))]
            ranked.sort(key=lambda row: (-row["safe_region_margin"], -(row.get("quality_score") or 0.0), row.get("candidate") or ""))
            shortlist = ranked[:shortlist_size]
        except Exception:
            shortlist = []

    out = {
        "packet": "P1152",
        "as_of": "2026-05-10",
        "allow_failures": allow_failures,
        "require_safe_region": require_safe_region,
        "rank_by_safe_margin": rank_by_safe_margin,
        "enforce_shortlist_consistency": enforce_shortlist_consistency,
        "ranking_triggered": ranking_triggered,
        "ranking_summary": str(ranking_summary_path) if ranking_triggered else None,
        "safe_margin_shortlist_size": shortlist_size if ranking_triggered else 0,
        "safe_margin_shortlist": shortlist,
        "safe_margin_shortlist_path": str(shortlist_output) if rank_by_safe_margin else None,
        "shortlist_consistency_enforced": False,
        "registry": str(registry),
        "total": len(results),
        "passed": passed,
        "failed": failed,
        "results": results,
    }

    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if rank_by_safe_margin:
        shortlist_payload = {
            "packet": "P1152_P1177",
            "as_of": out["as_of"],
            "registry": str(registry),
            "source_summary": str(out_path),
            "shortlist_size": shortlist_size,
            "shortlist": shortlist,
            "note": "Methodological shortlist only; no strict-core closure or QW-2191 discharge claim.",
        }
        shortlist_output.write_text(json.dumps(shortlist_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if enforce_shortlist_consistency:
            check_run = subprocess.run([sys.executable, str(ROOT / "p1177_shortlist_consistency_check.py"), str(out_path)], capture_output=True, text=True)
            out["shortlist_consistency_enforced"] = True
            out["shortlist_consistency_returncode"] = check_run.returncode
            out["shortlist_consistency_stdout"] = check_run.stdout.strip()
            out["shortlist_consistency_stderr"] = check_run.stderr.strip()
            out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            if check_run.returncode != 0:
                print("[P1152] shortlist consistency check failed.")
                raise SystemExit(2)
        else:
            out["shortlist_consistency_enforced"] = False
    print(f"[P1152] total={len(results)} passed={passed} failed={failed} wrote {out_path}")
    if rank_by_safe_margin:
        if shortlist:
            print("[P1152] Safe-margin shortlist for next physical tests:")
            for idx, row in enumerate(shortlist, start=1):
                print(f"  {idx}. {row.get('candidate')} margin={row.get('safe_region_margin')}")
        else:
            print("[P1152] Safe-margin shortlist unavailable (no passing candidates with computed margins).")
    if failed > 0 and not allow_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
