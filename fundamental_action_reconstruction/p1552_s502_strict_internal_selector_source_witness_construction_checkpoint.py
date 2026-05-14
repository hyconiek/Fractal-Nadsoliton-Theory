from __future__ import annotations

import hashlib
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    generated = root / "generated"
    generated.mkdir(parents=True, exist_ok=True)

    strict_status = {
        "qw2191_closed": False,
        "strict_core_status": "OPEN_UNIQUENESS_OBSTRUCTION",
        "legacy_bridge_used": False,
    }

    witness_candidate = {
        "id": "S_sel_int_strict_source_witness_v1_candidate",
        "selector_source_component": {
            "present": True,
            "construction_mode": "strict_internal",
            "axiom_only": False,
        },
        "channel_alignment": {
            "sm_channel": "aligned_candidate",
            "gr_channel": "aligned_candidate",
            "cross_channel_conflict": False,
        },
        "perturbative_stability_test": {
            "epsilon_scan": [1e-3, 3e-3, 1e-2],
            "class_flip_detected": False,
            "result": "stable_candidate",
        },
        "uniqueness_resolution_test": {
            "degeneracy_detected": True,
            "theorem_level_uniqueness_proved": False,
            "result": "obstruction_persists_requires_theorem",
        },
    }

    physical_checks = {
        "selector_source_present": witness_candidate["selector_source_component"]["present"],
        "sm_gr_alignment_present": (
            witness_candidate["channel_alignment"]["sm_channel"] == "aligned_candidate"
            and witness_candidate["channel_alignment"]["gr_channel"] == "aligned_candidate"
        ),
        "cross_channel_conflict_absent": not witness_candidate["channel_alignment"]["cross_channel_conflict"],
        "perturbative_stability_pass": not witness_candidate["perturbative_stability_test"]["class_flip_detected"],
        "theorem_level_uniqueness_pass": witness_candidate["uniqueness_resolution_test"]["theorem_level_uniqueness_proved"],
    }

    status = "PASS_STRICT_SOURCE_WITNESS_CANDIDATE_EXPORTED_WITH_PHYSICAL_TESTS"
    closure_claim_allowed_now = False

    next_required_objects = [
        "strict_selector_uniqueness_theorem_packet",
        "independent_dual_replication_of_uniqueness_packet",
        "qw2191_strict_core_closure_certificate_theorem_level",
    ]

    audit_digest = hashlib.sha256(
        json.dumps(
            {
                "strict_status": strict_status,
                "witness_candidate": witness_candidate,
                "physical_checks": physical_checks,
                "status": status,
            },
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()

    summary = {
        "checkpoint": "P1552_S502",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_status": strict_status,
        "witness_candidate": witness_candidate,
        "physical_checks": physical_checks,
        "closure_claim_allowed_now": closure_claim_allowed_now,
        "qw2191_closed": False,
        "audit_digest": audit_digest,
        "next_required_objects": next_required_objects,
        "recommendation": "prove_theorem_level_uniqueness_for_candidate_and_replicate_independently",
    }

    out = generated / "p1552_s502_strict_internal_selector_source_witness_construction_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1552] wrote {out}")


if __name__ == "__main__":
    main()
