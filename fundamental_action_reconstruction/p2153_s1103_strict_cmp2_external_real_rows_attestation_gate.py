#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
EXT = GEN / "cmp2_backend_rows_extension_v1.json"
OUT = GEN / "p2153_s1103_strict_cmp2_external_real_rows_attestation_gate.json"
MD = GEN / "p2153_s1103_strict_cmp2_external_real_rows_attestation_gate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    ext = load(EXT)

    att = (ext.get("external_confirmation", {}) if isinstance(ext, dict) else {}) or {}
    required = ["provider", "delivered_at_utc", "independent_confirmation_id", "confirmed"]
    missing = [k for k in required if k not in att]
    has_required = len(missing) == 0
    confirmed = bool(att.get("confirmed") is True)

    provenance_note = str(ext.get("provenance_note", ""))
    provenance_blocked = "constructed from existing repo backend-evidence rows" in provenance_note.lower()

    attested_external_real_rows = has_required and confirmed and not provenance_blocked

    result_kind = (
        "PASS_STRICT_CMP2_EXTERNAL_REAL_ROWS_ATTESTATION_GATE"
        if attested_external_real_rows
        else "OPEN_STRICT_CMP2_EXTERNAL_REAL_ROWS_ATTESTATION_GATE"
    )

    payload = {
        "schema_version": "p2153_s1103_v1",
        "packet_id": "P2153",
        "stage_id": "S1103",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "external_real_rows_attestation_gate": {
            "extension_file": str(EXT.relative_to(ROOT)),
            "external_confirmation": att,
            "required_confirmation_keys": required,
            "missing_confirmation_keys": missing,
            "confirmed": confirmed,
            "provenance_note": provenance_note,
            "provenance_blocked": provenance_blocked,
            "attested_external_real_rows": attested_external_real_rows,
            "scope_limit": "attestation gate only; no theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2154_candidate",
            "goal": "replace extension payload with externally delivered independently confirmed real rows and rerun P2147/P2150/P2144",
        },
        "gatekeeper_checks": {
            "attestation_gate_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2153 S1103: strict CMP2 external real-rows attestation gate",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Attested external real rows: `{attested_external_real_rows}`",
                f"- Missing confirmation keys: `{missing}`",
                "",
                "No theorem-grade closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
