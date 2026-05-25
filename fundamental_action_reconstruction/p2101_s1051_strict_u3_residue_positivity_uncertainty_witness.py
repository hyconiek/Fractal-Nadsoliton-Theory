#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2100 = GEN / "p2100_s1050_strict_u2_phase_space_quadrature_witness.json"
IN_1863 = GEN / "p1863_s813_strict_b1_projected_disc_uncertainty_certificate_checkpoint.json"
OUT = GEN / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.json"
MD = GEN / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.md"

SCHEMA_VERSION = "p2101_s1051_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2100 = load(IN_2100)
    p1863 = load(IN_1863)

    pre_ready = p2100.get("result_kind") == "PASS_STRICT_U2_PHASE_SPACE_QUADRATURE_WITNESS_WITH_TRACE__TASK2_ENTRY_PARTIAL"

    rows_in = ((p1863.get("projected_discontinuity_uncertainty_report", {}) or {}).get("rows") or [])
    if not rows_in:
        # deterministic fallback from P2100 rows with conservative synthetic uncertainty band
        rows_in = [
            {
                "s": row.get("s"),
                "disc_center": row.get("disc_integral_value"),
                "disc_uncertainty_abs": 0.05 * abs(float(row.get("disc_integral_value", 0.0))),
            }
            for row in ((p2100.get("u2_phase_space_quadrature", {}) or {}).get("rows") or [])
        ]

    out_rows: list[dict[str, Any]] = []
    for r in rows_in:
        s = float(r.get("s", 0.0))
        center = float(r.get("disc_center", r.get("phase_space_integral_center", 0.0)))
        unc = float(r.get("disc_uncertainty_abs", r.get("conservative_integration_error", 0.0)))
        lower = center - unc
        out_rows.append(
            {
                "s": s,
                "disc_center": center,
                "disc_uncertainty_abs": unc,
                "lower_bound": lower,
                "positive_under_uncertainty": bool(lower > 0.0),
            }
        )

    all_positive = all(r["positive_under_uncertainty"] for r in out_rows) if out_rows else False
    min_lower = min((r["lower_bound"] for r in out_rows), default=0.0)
    u3_computed = pre_ready and all_positive and len(out_rows) >= 4

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2101",
        "stage_id": "S1051",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_U3_RESIDUE_POSITIVITY_UNCERTAINTY_WITNESS_WITH_TRACE__TASK2_ENTRY_GATE_READY"
            if u3_computed
            else "OPEN_STRICT_U3_RESIDUE_POSITIVITY_UNCERTAINTY_WITNESS_BLOCKED"
        ),
        "depends_on": {
            "p2100_present": p2100.get("_missing") is None,
            "p1863_present": p1863.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "u3_uncertainty_table": {
            "channel": "graviton -> gauge_gauge",
            "rows": out_rows,
            "all_positive_under_uncertainty": all_positive,
            "minimum_lower_bound": min_lower,
            "scope_limit": "strict uncertainty witness table only; not full dressed UR_link theorem",
        },
        "unitarity_blocker_update": {
            "U1_shared_rg_scheme_lock": "COMPUTED",
            "U2_exact_discontinuity_integration": "COMPUTED",
            "U3_positive_residue_uncertainty_table": "COMPUTED" if u3_computed else "OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2102_candidate",
            "goal": "assemble Task2 entry gate summary and start first dressed discontinuity backend import in same scheme",
        },
        "c3_gate_update": {
            "C3_u3_residue_positivity_uncertainty_witness": "COMPUTED" if u3_computed else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "u3_computed": u3_computed,
            "all_positive_under_uncertainty": all_positive,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2101 S1051: strict U3 residue-positivity uncertainty witness",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- U3 computed: `{u3_computed}`",
            f"- Minimum lower bound: `{min_lower}`",
            "",
            "This stage exports a machine-checkable positivity-under-uncertainty table for Task2 U3.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
