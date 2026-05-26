#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2180 = GEN / "p2180_s1130_strict_qw2191_independent_audit_replay_against_handoff_contract.json"
OUT = GEN / "p2181_s1131_strict_qw2191_replay_certificate_freeze.json"
MD = GEN / "p2181_s1131_strict_qw2191_replay_certificate_freeze.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2180 = load(IN_2180)

    replay = p2180.get("strict_qw2191_independent_audit_replay_against_handoff_contract", {}) or {}
    checks = replay.get("checks", {}) or {}
    all_checks_pass = bool(replay.get("all_checks_pass", False))

    certificate = {
        "certificate_id": "QW2191_REPLAY_CERTIFICATE_V1",
        "source_packet": str(IN_2180.relative_to(ROOT)),
        "checks_snapshot": checks,
        "certificate_status": "CERTIFIED_WITH_TRACE" if all_checks_pass else "NOT_CERTIFIED",
        "scope": "independent audit replay consistency only",
    }

    result_kind = (
        "PASS_STRICT_QW2191_REPLAY_CERTIFICATE_FREEZE_WITH_TRACE"
        if all_checks_pass
        else "OPEN_STRICT_QW2191_REPLAY_CERTIFICATE_FREEZE_BLOCKED"
    )

    payload = {
        "schema_version": "p2181_s1131_v1",
        "packet_id": "P2181",
        "stage_id": "S1131",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_replay_certificate_freeze": {
            "replay_certificate": certificate,
            "scope_limit": "certificate freeze only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2182_candidate",
            "goal": "run contradiction-targeted delta audit between historical release claims and current strict-lane certificates",
        },
        "gatekeeper_checks": {
            "replay_certificate_exported": True,
            "replay_certified": all_checks_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2180.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2180.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2181 S1131: strict QW-2191 replay certificate freeze",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Replay certified: `{all_checks_pass}`",
                f"- Certificate id: `{certificate['certificate_id']}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
