#!/usr/bin/env python3
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    out = root / "generated"
    out.mkdir(parents=True, exist_ok=True)

    data = {
        "id": "H36",
        "date": "2026-03-06",
        "status": "PASS_PARTIAL_NO_STRICT_DIRECTED_AXIS_SELECTION",
        "result": "strict_core_supports_only_a_coordinate_level_undirected_axis_representative_u_psi0_pair1_inside_pair1_and_not_a_strict_directed_orientation_selection",
        "frontier": "H36_B1",
        "frontier_text": "strict core supports only a coordinate-level undirected axis representative u_psi0_pair1 inside pair1 and contains no strict argument selecting a directed orientation on that axis",
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "no claim that psi0 selects a physically directed axis",
            "no claim that u_psi0_pair1 and -u_psi0_pair1 are physically distinguished by strict core",
            "no claim that QW-2191 is discharged",
        ],
    }

    summary = {
        "id": data["id"],
        "status": data["status"],
        "frontier": data["frontier"],
        "result": data["result"],
    }

    (out / "h36_directed_axis_orientation_audit.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8"
    )
    (out / "h36_directed_axis_orientation_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
