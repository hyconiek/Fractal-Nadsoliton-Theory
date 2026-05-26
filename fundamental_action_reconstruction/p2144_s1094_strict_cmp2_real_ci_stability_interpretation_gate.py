#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2143 = GEN / "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.json"
IN_2135 = GEN / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.json"
OUT = GEN / "p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.json"
MD = GEN / "p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.md"


def load(p: Path) -> dict[str, Any]:
    if not p.exists():
        return {"_missing": str(p)}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2143 = load(IN_2143)
    p2135 = load(IN_2135)

    memo = p2143.get("real_ci_stability_readiness_memo", {}) or {}
    ready = bool(memo.get("fully_ready", False))

    var_obj = (p2135.get("merged_real_block_variant_stability_audit", {}) or {}).get("variants", {}) or {}
    spread = (p2135.get("merged_real_block_variant_stability_audit", {}) or {}).get("inflation_spread_across_variants")

    interpretation = "BLOCKED_PENDING_REAL_EXTENSION_ROWS" if not ready else (
        "READY_WITH_FINITE_SPREAD" if spread is not None else "READY_BUT_SPREAD_UNAVAILABLE"
    )

    payload = {
        "schema_version": "p2144_s1094_v1",
        "packet_id": "P2144",
        "stage_id": "S1094",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_REAL_CI_STABILITY_INTERPRETATION_GATE_WITH_TRACE" if ready else "OPEN_STRICT_CMP2_REAL_CI_STABILITY_INTERPRETATION_GATE_BLOCKED",
        "real_ci_stability_interpretation_gate": {
            "source_memo": str(IN_2143.relative_to(ROOT)),
            "source_variant_audit": str(IN_2135.relative_to(ROOT)),
            "ready": ready,
            "variant_count": len(var_obj),
            "inflation_spread_across_variants": spread,
            "interpretation": interpretation,
            "scope_limit": "interpretation gate only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2145_candidate",
            "goal": "if blocked: provide real extension rows and rerun P2142-P2144; if ready: freeze non-synthetic CI stability interpretation for archival handoff",
        },
        "gatekeeper_checks": {
            "ready_flag": ready,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2144 S1094: strict CMP2 real CI stability interpretation gate\n\n"
        f"- Status: `{payload['status']}`\n"
        f"- Result kind: `{payload['result_kind']}`\n"
        f"- Interpretation: `{interpretation}`\n\n"
        "No theorem-grade global closure claim is made.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
