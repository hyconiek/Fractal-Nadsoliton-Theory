"""P1328: formal export registry update for O3 theorem obligations."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1328_formal_export_registry_report_v1.json"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1329 = load("p1329_residual_loophole_elimination_proof_attempt_report_v1.json")
    p1338 = load("p1338_global_l2_export_final_attempt_report_v1.json")

    l2_exported = p1338.get("global_l2_exportable", False)

    payload = {
        "packet_id": "P1328_FORMAL_EXPORT_REGISTRY_REPORT_V1",
        "date_utc": "2026-05-12",
        "exports": {
            "theorem_V4_selector_consistency_formally_exported": True,
            "residual_loophole_elimination_formally_exported": l2_exported,
        },
        "conditional_support": {
            "residual_loophole_elimination_margin_domain_support": p1329.get(
                "l2_proven_on_admissible_margin_domain", False
            )
        },
        "formal_export_progress": "2_of_2" if l2_exported else "1_of_2",
        "theorem_ready": l2_exported,
        "qw2191_strict_status": "CLOSED" if l2_exported else "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
