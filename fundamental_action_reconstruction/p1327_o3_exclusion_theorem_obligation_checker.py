"""P1327: theorem-level O3-EXCLUSION obligation checker with v4 integration."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1327_o3_exclusion_theorem_obligation_checker_report_v1.json"


def load(name: str) -> dict:
    path = GEN / name
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p1326 = load("p1326_v4_replay_adversarial_reintegration_report_v1.json")
    p1328 = load("p1328_formal_export_registry_report_v1.json")
    exports = p1328.get("exports", {})

    obligations = {
        "V4_replay_pass": p1326.get("replay_pass") is True,
        "V4_adversarial_pass": p1326.get("adversarial_pass") is True,
        "V4_reintegration_ready": p1326.get("o3_reintegration_ready") is True,
        "theorem_V4_selector_consistency_formally_exported": exports.get(
            "theorem_V4_selector_consistency_formally_exported", False
        ),
        "residual_loophole_elimination_formally_exported": exports.get(
            "residual_loophole_elimination_formally_exported", False
        ),
    }

    theorem_ready = all(obligations.values())

    payload = {
        "packet_id": "P1327_O3_EXCLUSION_THEOREM_OBLIGATION_CHECKER_REPORT_V1",
        "date_utc": "2026-05-12",
        "obligations": obligations,
        "theorem_ready": theorem_ready,
        "missing_obligations": [k for k, v in obligations.items() if not v],
        "qw2191_strict_status": "CLOSED" if theorem_ready else "NOT_CLOSED",
    }

    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
