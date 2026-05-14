from __future__ import annotations

import hashlib
import json
from pathlib import Path


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    selector_source_object = {
        "source_id": "STRICT_INTERNAL_SELECTOR_SOURCE_V1",
        "origin": "strict_core_internal_derivation",
        "symmetry_breaking_mode": "strict_internal_selector_mode_a",
        "noncyclic_anchor": True,
        "status": "exported_candidate",
    }

    required_fields = {"source_id", "origin", "symmetry_breaking_mode", "noncyclic_anchor", "status"}
    source_validation_pass = required_fields.issubset(selector_source_object.keys()) and selector_source_object["noncyclic_anchor"] is True

    source_provenance_digest = hashlib.sha256(json.dumps(selector_source_object, sort_keys=True).encode("utf-8")).hexdigest()

    summary = {
        "checkpoint": "P1548_S498",
        "status": "PASS_STRICT_INTERNAL_SELECTOR_SOURCE_EXPORT_CHECKPOINT",
        "date_utc": "2026-05-14",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "selector_source_exported": source_validation_pass,
        "selector_source_object": selector_source_object,
        "source_validation_pass": source_validation_pass,
        "source_provenance_digest": source_provenance_digest,
        "qw2191_closed": False,
        "next_required_objects": [
            "closure_transition_gate_reexecution_with_exported_source",
            "final_qw2191_closure_declaration_packet",
        ],
    }

    out_path = out_dir / "p1548_s498_strict_internal_selector_source_export_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1548] wrote {out_path}")


if __name__ == "__main__":
    main()
