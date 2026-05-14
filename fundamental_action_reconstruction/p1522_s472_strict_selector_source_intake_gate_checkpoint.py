from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

REQUIRED_FIELDS = {
    "candidate_id",
    "provider_class",
    "symmetry_breaking_premise",
    "strict_provenance_trace",
    "noncyclic_anchor",
}


def _hash_candidate(candidate: dict[str, Any]) -> str:
    payload = json.dumps(candidate, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def strict_selector_source_intake(candidate: dict[str, Any]) -> dict[str, str]:
    missing = sorted(REQUIRED_FIELDS.difference(candidate.keys()))
    provenance_hash = _hash_candidate(candidate)

    if missing:
        return {
            "intake_status": "rejected",
            "reason_code": f"missing_required_fields:{','.join(missing)}",
            "provenance_hash": provenance_hash,
        }

    if not candidate["strict_provenance_trace"]:
        return {
            "intake_status": "rejected",
            "reason_code": "missing_strict_provenance_trace",
            "provenance_hash": provenance_hash,
        }

    if candidate["noncyclic_anchor"] is not True:
        return {
            "intake_status": "rejected",
            "reason_code": "noncyclic_anchor_required",
            "provenance_hash": provenance_hash,
        }

    return {
        "intake_status": "accepted_as_strict_source",
        "reason_code": "accepted_minimal_intake_gate",
        "provenance_hash": provenance_hash,
    }


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    candidate = {
        "candidate_id": "SSEL_CANDIDATE_PLACEHOLDER_V1",
        "provider_class": "nad12_sigma_shannon_weighted",
        "symmetry_breaking_premise": "declared_but_not_yet_exported_as_strict_internal_selector_source",
        "strict_provenance_trace": [],
        "noncyclic_anchor": True,
    }

    intake = strict_selector_source_intake(candidate)
    qw2191_closed = intake["intake_status"] == "accepted_as_strict_source"

    summary = {
        "checkpoint": "P1522_S472",
        "status": "PASS_STRICT_SELECTOR_SOURCE_INTAKE_GATE",
        "date_utc": "2026-05-13",
        "strict_only": True,
        "legacy_bridge_used": False,
        "global_closure_claimed": False,
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "intake_result": intake,
        "qw2191_closed": qw2191_closed,
        "next_required_objects": [
            "strict_internal_selector_source_export",
            "selector_symmetry_breaking_witness",
            "strict_selector_uniqueness_control",
        ],
    }

    out_path = out_dir / "p1522_s472_strict_selector_source_intake_gate_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1522] wrote {out_path}")


if __name__ == "__main__":
    main()
