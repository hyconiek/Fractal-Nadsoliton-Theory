#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_P1375 = GEN / "p1375_l_b1_02_scale_scheme_transport_invariance_attempt_summary.json"
OUT = GEN / "p1376_full_ci_export_for_ci_coefficient_class_summary.json"
OUT_SCHEMA = GEN / "p1376_ci_coefficient_class_schema_v1.json"


def main() -> None:
    p1375 = json.loads(IN_P1375.read_text(encoding="utf-8"))

    schema = {
        "schema_id": "p1376_ci_coefficient_class_schema_v1",
        "as_of": "2026-05-12",
        "coefficients": ["c_g3", "c_g2", "c_g1", "c_mix"],
        "semantics": {
            "c_g3": "SU(3) gauge slot coefficient",
            "c_g2": "SU(2) gauge slot coefficient",
            "c_g1": "U(1) gauge slot coefficient",
            "c_mix": "gauge-mixing slot coefficient"
        },
        "population_status": "NOT_YET"
    }

    summary = {
        "artifact": "p1376_full_ci_export_for_ci_coefficient_class_summary",
        "as_of": "2026-05-12",
        "input_dependency": IN_P1375.name,
        "input_l_b1_02_status": p1375.get("l_b1_02_status"),
        "schema_export": OUT_SCHEMA.name,
        "ci_class_defined": True,
        "ci_class_population_complete": False,
        "l_b1_02_readiness": "PARTIAL",
        "forbidden_claim": "do_not_claim_transport_invariance_before_full_population",
        "next_packet": "P1377_CI_CLASS_POPULATION_AND_FIRST_TRANSPORT_DRIFT_RUN"
    }

    GEN.mkdir(parents=True, exist_ok=True)
    OUT_SCHEMA.write_text(json.dumps(schema, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    OUT.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[p1376] wrote {OUT_SCHEMA}")
    print(f"[p1376] wrote {OUT}")


if __name__ == "__main__":
    main()
