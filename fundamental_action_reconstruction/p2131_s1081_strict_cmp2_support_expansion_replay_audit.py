#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2126 = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
IN_2127 = GEN / "p2127_s1077_strict_cmp2_bootstrap_backend_evidence_stresstest.json"
OUT = GEN / "p2131_s1081_strict_cmp2_support_expansion_replay_audit.json"
MD = GEN / "p2131_s1081_strict_cmp2_support_expansion_replay_audit.md"

SCHEMA_VERSION = "p2131_s1081_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
RNG_SEED = 2131
N_BOOT = 400
MODELS = ["M1_nn", "M2_monotone", "M3_nn_tiebreak", "M4_monotone_penalized"]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def ci_width(ci):
    if ci[0] is None or ci[1] is None:
        return None
    return float(ci[1]-ci[0])


def top_model(r):
    w = r.get("posterior_weights_backend_evidence", {}) or {}
    return max(w.items(), key=lambda kv: kv[1])[0]


def iid_boot(rows, rng):
    n=len(rows)
    cov=np.array([1.0 if r.get("posterior_predictive_covered",False) else 0.0 for r in rows])
    top=[top_model(r) for r in rows]
    boot_cov=[]; boot_top={m:[] for m in MODELS}
    for _ in range(N_BOOT):
        idx=rng.integers(0,n,size=n)
        boot_cov.append(float(np.mean(cov[idx])))
        s=[top[i] for i in idx]
        for m in MODELS:
            boot_top[m].append(s.count(m)/n)
    return {
      "coverage_rate_ci95":[float(np.quantile(boot_cov,0.025)), float(np.quantile(boot_cov,0.975))],
      "posterior_top_model_rank_stability":{m:{"ci95":[float(np.quantile(v,0.025)), float(np.quantile(v,0.975))]} for m,v in boot_top.items()}
    }


def block_boot(rows, rng):
    by={}
    for r in rows:
        by.setdefault(float(r.get("backend_s",0.0)), []).append(r)
    keys=sorted(by.keys())
    cov=[]; top={m:[] for m in MODELS}
    for _ in range(N_BOOT):
        sampled=[keys[i] for i in rng.integers(0,len(keys),size=len(keys))]
        sample=[rr for k in sampled for rr in by[k]]
        m=len(sample)
        c=[1.0 if r.get("posterior_predictive_covered",False) else 0.0 for r in sample]
        cov.append(float(np.mean(c)) if m else 0.0)
        tm=[top_model(r) for r in sample]
        for model in MODELS:
            top[model].append(tm.count(model)/m if m else 0.0)
    return {
      "coverage_rate_ci95":[float(np.quantile(cov,0.025)), float(np.quantile(cov,0.975))],
      "posterior_top_model_rank_stability":{m:{"ci95":[float(np.quantile(v,0.025)), float(np.quantile(v,0.975))]} for m,v in top.items()}
    }


def expand_rows(rows, mult, rng):
    n=len(rows)
    idx=rng.integers(0,n,size=n*mult)
    return [rows[i] for i in idx]


def main():
    GEN.mkdir(exist_ok=True)
    p2126=load(IN_2126); p2127=load(IN_2127)
    pre_ready = p2126.get("result_kind") == "PASS_STRICT_CMP2_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_WITH_TRACE"
    iid_ready = p2127.get("result_kind") == "PASS_STRICT_CMP2_BOOTSTRAP_BACKEND_EVIDENCE_STRESSTEST_WITH_TRACE"
    rows=((p2126.get("backend_evidence_weighted_posterior_predictive_calibration_audit",{}) or {}).get("rows") or [])
    results={}
    if pre_ready and iid_ready and rows:
        for mult in [2,3,5,8]:
            rng=np.random.default_rng(RNG_SEED+mult)
            erows=expand_rows(rows,mult,rng)
            iid=iid_boot(erows,rng)
            blk=block_boot(erows,rng)
            iw=ci_width(iid["coverage_rate_ci95"]); bw=ci_width(blk["coverage_rate_ci95"])
            infl=float(bw/iw) if iw and iw>0 else None
            results[f"x{mult}"]={
                "n_rows":len(erows),
                "iid_coverage_rate_ci95":iid["coverage_rate_ci95"],
                "block_coverage_rate_ci95":blk["coverage_rate_ci95"],
                "coverage_ci_width_inflation_block_over_iid":infl,
            }

    payload={
      "schema_version":SCHEMA_VERSION,
      "packet_id":"P2131","stage_id":"S1081","produced_by":Path(__file__).name,
      "timestamp_utc":TIMESTAMP_UTC,
      "status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"PASS_STRICT_CMP2_SUPPORT_EXPANSION_REPLAY_AUDIT_WITH_TRACE" if results else "OPEN_STRICT_CMP2_SUPPORT_EXPANSION_REPLAY_AUDIT_BLOCKED",
      "support_expansion_replay_audit":{
        "n_bootstrap":N_BOOT,"rng_seed":RNG_SEED,"synthetic_expansion_multipliers":list(results.keys()),"results":results,
        "scope_limit":"synthetic resampling replay only; not substitute for new backend acquisition and not theorem-grade closure"
      },
      "recommended_next_honest_step": {"id":"P2132_candidate","goal":"acquire real additional CMP2 backend rows (non-synthetic) and rerun P2127-P2131"},
      "gatekeeper_checks":{"preconditions_ready":pre_ready and iid_ready,"replay_executed":bool(results),"no_toe_closure_claimed":True,
      "full_d3_covariance_transport_proven":False,"full_cutkosky_closure_proven":False,"c3_theorem_proven":False},
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    MD.write_text("# P2131 S1081: strict CMP2 support-expansion replay audit\n\n- Status: `OPEN_PARTIAL_PROGRESS_WITH_TRACE`\n- Synthetic multipliers: `"+str(list(results.keys()))+"`\n\nSynthetic replay only; no theorem-grade closure claim.\n",encoding='utf-8')

if __name__=='__main__':
    main()
