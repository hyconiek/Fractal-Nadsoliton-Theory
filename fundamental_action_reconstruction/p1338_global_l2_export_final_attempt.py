"""P1338: global L2 export final attempt by bridging margin proof and v2 source."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1338_global_l2_export_final_attempt_report_v1.json"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1329 = load("p1329_residual_loophole_elimination_proof_attempt_report_v1.json")
    p1336 = load("p1336_internal_source_v2_tiebreak_family_search_report_v1.json")

    margin_domain_proven = p1329.get("l2_proven_on_admissible_margin_domain", False)
    boundary_consistency = p1336.get("strict_internal_source_exportable", False)

    global_l2_exportable = margin_domain_proven and boundary_consistency

    payload = {
        "packet_id": "P1338_GLOBAL_L2_EXPORT_FINAL_ATTEMPT_REPORT_V1",
        "date_utc": "2026-05-12",
        "inputs": {
            "margin_domain_proven": margin_domain_proven,
            "boundary_consistency_from_v2_source": boundary_consistency,
        },
        "global_l2_exportable": global_l2_exportable,
        "status": "GLOBAL_L2_EXPORTED" if global_l2_exportable else "GLOBAL_L2_NOT_EXPORTED",
        "qw2191_strict_status": "CLOSED" if global_l2_exportable else "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
