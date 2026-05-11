#!/usr/bin/env python3
"""P1169 end-to-end integration test for top P1166 candidate.

Runs P1151 on top candidate and P1152 registry (without allow-failures) to
verify full pipeline admissibility as a working candidate object.
"""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'


def run(cmd):
    p=subprocess.run(cmd,capture_output=True,text=True)
    return {"cmd":cmd,"returncode":p.returncode,"stdout":p.stdout.strip(),"stderr":p.stderr.strip()}


def main():
    args = sys.argv[1:]
    strict_e2e = "--strict-e2e" in args
    require_shortlist_consistency = "--require-shortlist-consistency" in args
    require_safe_region = "--require-safe-region" in args
    require_out_of_locality_robustness = "--require-out-of-locality-robustness" in args
    robustness_threshold = 0.99
    if "--robustness-threshold" in args:
        idx = args.index("--robustness-threshold")
        if idx + 1 < len(args):
            robustness_threshold = float(args[idx + 1])
    if strict_e2e:
        require_shortlist_consistency = True
        require_safe_region = True
    candidate_override = None
    registry_override = None
    if "--candidate" in args:
        i = args.index("--candidate")
        if i + 1 < len(args):
            candidate_override = Path(args[i + 1]).resolve()
    if "--registry" in args:
        i = args.index("--registry")
        if i + 1 < len(args):
            registry_override = Path(args[i + 1]).resolve()
    cand=(candidate_override or (GEN/'p1169_candidate_top_from_p1166.json').resolve())
    reg=(registry_override or (GEN/'p1169_candidate_registry_top.json').resolve())

    step1=run([sys.executable,str(ROOT/'p1151_strict_selector_pipeline_runner.py'),str(cand)])
    s1=json.loads((GEN/'p1151_strict_selector_pipeline_runner_summary.json').read_text())

    step2_cmd = [sys.executable,str(ROOT/'p1152_strict_candidate_registry_runner.py'),str(reg)]
    if require_safe_region:
        step2_cmd.append("--require-safe-region")
    if require_shortlist_consistency:
        step2_cmd.extend(["--rank-by-safe-margin", "--enforce-shortlist-consistency"])
    step2=run(step2_cmd)
    s2=json.loads((GEN/'p1152_strict_candidate_registry_runner_summary.json').read_text())

    shortlist_consistency_ok = True
    if require_shortlist_consistency:
        shortlist_consistency_ok = bool(s2.get("shortlist_consistency_enforced")) and int(s2.get("shortlist_consistency_returncode", 1)) == 0

    all_registry_pass = s2.get('total', 0) > 0 and s2.get('passed', 0) == s2.get('total', -1)
    robustness_ok = True
    robustness_summary = None
    robustness_run = None
    if require_out_of_locality_robustness:
        robustness_run = run([sys.executable, str(ROOT/'p1171_out_of_locality_robustness_probe.py'), str(cand)])
        rsum = json.loads((GEN/'p1171_out_of_locality_robustness_probe_summary.json').read_text())
        robustness_summary = {"robust_cases": rsum.get("robust_cases"), "cases": rsum.get("cases"), "robust_fraction": rsum.get("robust_fraction")}
        robustness_ok = isinstance(rsum.get("robust_fraction"), (int, float)) and float(rsum.get("robust_fraction")) >= robustness_threshold

    out={"packet":"P1169","as_of":"2026-05-10",
         "candidate":str(cand),"registry":str(reg),
         "pipeline_run":step1,"pipeline_summary":{"overall_pass":s1.get('overall_pass'),"failed_step":s1.get('failed_step')},
         "registry_run":step2,"registry_summary":{"total":s2.get('total'),"passed":s2.get('passed'),"failed":s2.get('failed'), "shortlist_consistency_enforced": s2.get('shortlist_consistency_enforced'), "shortlist_consistency_returncode": s2.get('shortlist_consistency_returncode')},
         "strict_e2e": strict_e2e,
         "require_shortlist_consistency": require_shortlist_consistency,
         "require_safe_region": require_safe_region,
         "require_out_of_locality_robustness": require_out_of_locality_robustness,
         "robustness_threshold": robustness_threshold,
         "out_of_locality_robustness_run": robustness_run,
         "out_of_locality_robustness_summary": robustness_summary,
         "integrated_pass": bool(s1.get('overall_pass')) and (s2.get('failed',1)==0) and shortlist_consistency_ok and ((not strict_e2e) or all_registry_pass) and robustness_ok,
         "note":"Integration validation only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1169_e2e_candidate_integration_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1169] integrated_pass={out['integrated_pass']} wrote {outp}")

if __name__=='__main__':
    main()
