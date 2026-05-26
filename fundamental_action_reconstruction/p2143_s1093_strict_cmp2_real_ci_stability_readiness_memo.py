#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2142 = GEN / "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.json"
OUT = GEN / "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.json"
MD = GEN / "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2142 = load(IN_2142)
    h = (p2142.get("real_data_handoff_and_rerun_executor", {}) or {})
    rk = (h.get("post_rerun_result_kinds", {}) or {})

    ready = bool(h.get("fully_ready_for_real_ci_stability_assessment", False))
    blocked = [k for k in ["p2132", "p2133", "p2134", "p2135"] if not str(rk.get(k, "OPEN")).startswith("PASS_")]

    payload = {
        "schema_version": "p2143_s1093_v1",
        "packet_id": "P2143",
        "stage_id": "S1093",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_REAL_CI_STABILITY_READINESS_MEMO_WITH_TRACE" if ready else "OPEN_STRICT_CMP2_REAL_CI_STABILITY_READINESS_MEMO_BLOCKED",
        "real_ci_stability_readiness_memo": {
            "source": str(IN_2142.relative_to(ROOT)),
            "fully_ready": ready,
            "blocked_stages": blocked,
            "post_rerun_result_kinds": rk,
            "scope_limit": "readiness memo only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2144_candidate",
            "goal": "if blocked: deliver real extension rows and rerun P2142; if ready: freeze final non-synthetic CI stability interpretation",
        },
        "gatekeeper_checks": {
            "readiness_memo_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        f"# P2143 S1093: strict CMP2 real CI stability readiness memo\n\n- Status: `{payload['status']}`\n- Result kind: `{payload['result_kind']}`\n- Fully ready: `{ready}`\n\nNo theorem-grade global closure claim is made.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
