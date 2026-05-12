"""P1337: full O3 checker refresh after P1336 internal-source export."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1337_full_o3_checker_refresh_report_v1.json"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1326 = load("p1326_v4_replay_adversarial_reintegration_report_v1.json")
    p1328 = load("p1328_formal_export_registry_report_v1.json")
    p1336 = load("p1336_internal_source_v2_tiebreak_family_search_report_v1.json")

    obligations = {
        "v4_replay_pass": p1326.get("replay_pass") is True,
        "v4_adversarial_pass": p1326.get("adversarial_pass") is True,
        "v4_consistency_l1_exported": p1328.get("exports", {}).get("theorem_V4_selector_consistency_formally_exported", False),
        "strict_internal_source_exported": p1336.get("strict_internal_source_exportable", False),
        "global_l2_exported": p1328.get("exports", {}).get("residual_loophole_elimination_formally_exported", False),
    }

    pass_count = sum(1 for v in obligations.values() if v)
    ready = all(obligations.values())

    payload = {
        "packet_id": "P1337_FULL_O3_CHECKER_REFRESH_REPORT_V1",
        "date_utc": "2026-05-12",
        "obligations": obligations,
        "pass_count": pass_count,
        "total_obligations": len(obligations),
        "o3_strict_ready": ready,
        "missing": [k for k, v in obligations.items() if not v],
        "qw2191_strict_status": "CLOSED" if ready else "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
