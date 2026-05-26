#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2149 = GEN / "p2149_s1099_strict_cmp2_real_extension_delivery_checklist_packet.json"
IN_2147 = GEN / "p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.json"
OUT = GEN / "p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.json"
MD = GEN / "p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2149 = load(IN_2149)
    p2147 = load(IN_2147)

    blocked = ((p2147.get("real_data_required_rerun_checkpoint", {}) or {}).get("result_kinds", {}) or {})
    unresolved = [k for k in ["p2132", "p2133", "p2134", "p2135"] if not str(blocked.get(k, "OPEN")).startswith("PASS_")]

    payload = {
        "schema_version": "p2150_s1100_v1",
        "packet_id": "P2150",
        "stage_id": "S1100",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "OPEN_STRICT_CMP2_REAL_DATA_UNLOCK_ATTEMPT_REGISTER" if unresolved else "PASS_STRICT_CMP2_REAL_DATA_UNLOCK_ATTEMPT_REGISTER",
        "real_data_unlock_attempt_register": {
            "source_checklist": str(IN_2149.relative_to(ROOT)),
            "source_checkpoint": str(IN_2147.relative_to(ROOT)),
            "unresolved_stages": unresolved,
            "attempt_note": "external real data delivery still required" if unresolved else "non-synthetic path unlocked",
            "scope_limit": "attempt register only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2151_candidate",
            "goal": "deliver real extension rows and rerun P2147, then confirm p2132-p2135 PASS before interpretation",
        },
        "gatekeeper_checks": {
            "register_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2150 S1100: strict CMP2 real data unlock attempt register\n\n"
        f"- Result kind: `{payload['result_kind']}`\n"
        f"- Unresolved stages: `{unresolved}`\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
