from __future__ import annotations

import json
from pathlib import Path

from p1522_s472_strict_selector_source_intake_gate_checkpoint import strict_selector_source_intake


def build_minimal_strict_provenance_trace() -> list[dict[str, object]]:
    return [
        {
            "step_id": "trace_step_1",
            "artifact_ref": "P1521_S471_contract",
            "operation": "export_osplit_io_contract",
            "strict_guardrail_ok": True,
        },
        {
            "step_id": "trace_step_2",
            "artifact_ref": "P1522_S472_intake_gate",
            "operation": "selector_candidate_schema_check",
            "strict_guardrail_ok": True,
        },
        {
            "step_id": "trace_step_3",
            "artifact_ref": "P1523_S473_disambiguation",
            "operation": "local_closure_not_equal_selector_source",
            "strict_guardrail_ok": True,
        },
    ]


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    trace = build_minimal_strict_provenance_trace()

    candidate = {
        "candidate_id": "SSEL_CANDIDATE_STRICT_TRACE_V1",
        "provider_class": "nad12_sigma_shannon_weighted",
        "symmetry_breaking_premise": "explicit_candidate_premise_exported_for_intake_testing",
        "strict_provenance_trace": trace,
        "noncyclic_anchor": True,
    }

    intake = strict_selector_source_intake(candidate)

    summary = {
        "checkpoint": "P1524_S474",
        "status": "PASS_STRICT_PROVENANCE_TRACE_BUILDER",
        "date_utc": "2026-05-13",
        "route": "F_nadsoliton_to_LSM_plus_LGR",
        "strict_only": True,
        "legacy_bridge_used": False,
        "intake_result": intake,
        "intake_pass": intake["intake_status"] == "accepted_as_strict_source",
        "selector_uniqueness_witness_present": False,
        "qw2191_closed": False,
        "closure_note": "intake_pass_is_not_selector_closure",
        "next_required_objects": [
            "strict_selector_uniqueness_witness",
            "selector_resolution_theorem_or_internal_source_upgrade",
        ],
    }

    out_path = out_dir / "p1524_s474_strict_provenance_trace_builder_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1524] wrote {out_path}")


if __name__ == "__main__":
    main()
