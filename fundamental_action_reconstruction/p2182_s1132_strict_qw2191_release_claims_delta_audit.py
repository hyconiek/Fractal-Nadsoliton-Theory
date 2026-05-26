#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2181 = GEN / "p2181_s1131_strict_qw2191_replay_certificate_freeze.json"
IN_R81 = ROOT.parent / "RELEASE_8_1_TEXTBOOK_EDITION_EN_PL.md"
OUT = GEN / "p2182_s1132_strict_qw2191_release_claims_delta_audit.json"
MD = GEN / "p2182_s1132_strict_qw2191_release_claims_delta_audit.md"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2181 = load_json(IN_2181)
    r81 = load_text(IN_R81)

    # Minimal historical-claim probes (scoped text markers only).
    r81_has_scoped_closed = "global_closure_theorem_status = CLOSED" in r81
    r81_has_scope_boundary = "within declared Release-8 scope" in r81

    replay_cert = (p2181.get("strict_qw2191_replay_certificate_freeze", {}) or {}).get("replay_certificate", {}) or {}
    replay_certified = replay_cert.get("certificate_status") == "CERTIFIED_WITH_TRACE"

    # Delta audit focuses on scope mismatch risk only.
    delta_findings = []
    if r81_has_scoped_closed and not r81_has_scope_boundary:
        delta_findings.append(
            {
                "id": "D1_release_scope_boundary_missing",
                "severity": "high",
                "detail": "Release text has CLOSED marker without explicit scope boundary marker.",
            }
        )

    if replay_certified and r81_has_scoped_closed:
        delta_findings.append(
            {
                "id": "D2_scoped_closed_vs_lane_replay",
                "severity": "info",
                "detail": "Historical scoped closure and current-lane replay certificate coexist; treat as non-identical scopes.",
            }
        )

    blocking_findings = [f for f in delta_findings if f.get("severity") == "high"]
    contradiction_free = len(blocking_findings) == 0

    result_kind = (
        "PASS_STRICT_QW2191_RELEASE_CLAIMS_DELTA_AUDIT_WITH_TRACE"
        if contradiction_free
        else "OPEN_STRICT_QW2191_RELEASE_CLAIMS_DELTA_AUDIT_BLOCKED"
    )

    payload = {
        "schema_version": "p2182_s1132_v1",
        "packet_id": "P2182",
        "stage_id": "S1132",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_release_claims_delta_audit": {
            "source_replay_certificate": str(IN_2181.relative_to(ROOT)),
            "source_release_text": str(IN_R81.relative_to(ROOT.parent)),
            "probes": {
                "release_has_closed_marker": r81_has_scoped_closed,
                "release_has_scope_boundary": r81_has_scope_boundary,
                "replay_certified": replay_certified,
            },
            "delta_findings": delta_findings,
            "blocking_findings": blocking_findings,
            "scope_limit": "scope-delta audit only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2183_candidate",
            "goal": "freeze scope-delta governance note and wire into next contradiction ledger refresh",
        },
        "gatekeeper_checks": {
            "release_delta_audit_exported": True,
            "blocking_findings_absent": contradiction_free,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2181.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2181.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2182 S1132: strict QW-2191 release-claims delta audit",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Blocking findings: `{len(blocking_findings)}`",
                f"- Total findings: `{len(delta_findings)}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
