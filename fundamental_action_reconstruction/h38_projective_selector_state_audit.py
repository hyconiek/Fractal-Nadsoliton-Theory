#!/usr/bin/env python3
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    out = root / "generated"
    out.mkdir(parents=True, exist_ok=True)

    data = {
        "id": "H38",
        "date": "2026-03-16",
        "status": "PASS_PARTIAL_PROJECTIVE_STATE_ONLY__GLOBAL_PROJECTIVE_STATE_OBJECT_EXPORTED_BUT_NO_SIGN_SENSITIVE_ORIENTATION_DATUM",
        "result": "strict_core_supports_at_most_a_local_projective_or_ray_level_selector_representative_on_pair1_and_does_not_furnish_a_physically_individuated_directed_selector_state",
        "frontier": "H38_B1",
        "frontier_text": "strict core supports at most a local projective or ray-level selector representative on pair1 and does not furnish a physically individuated directed selector state",
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "no claim that the projective/ray-level state is already the physically correct selector state",
            "no claim that a sign-sensitive directed selector state or observable is exported",
            "no claim that QW-2191 is discharged",
        ],
    }

    summary = {
        "id": data["id"],
        "status": data["status"],
        "frontier": data["frontier"],
        "result": data["result"],
    }

    (out / "h38_projective_selector_state_audit.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8"
    )
    (out / "h38_projective_selector_state_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
