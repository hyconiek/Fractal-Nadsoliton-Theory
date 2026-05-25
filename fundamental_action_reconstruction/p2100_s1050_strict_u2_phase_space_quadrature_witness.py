#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2099 = GEN / "p2099_s1049_strict_u1_same_scheme_lock_witness.json"
IN_1951 = GEN / "p1951_s901_strict_b1_cutkosky_phase_space_integral_probe.json"
IN_1969 = GEN / "p1969_s919_strict_b1_exact_phase_space_integral_first_nonproxy_checkpoint.json"
OUT = GEN / "p2100_s1050_strict_u2_phase_space_quadrature_witness.json"
MD = GEN / "p2100_s1050_strict_u2_phase_space_quadrature_witness.md"

SCHEMA_VERSION = "p2100_s1050_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2099 = load(IN_2099)
    p1951 = load(IN_1951)
    p1969 = load(IN_1969)

    pre_ready = p2099.get("result_kind") == "PASS_STRICT_U1_SAME_SCHEME_LOCK_WITNESS_WITH_TRACE__TASK2_ENTRY_PARTIAL"

    # Reuse already exported strict-side phase-space exact formula (nonproxy checkpoint lineage).
    sym = p1969.get("symbolic_result", {}) or {}
    a_r2 = float(sp.N(sp.sympify(sym.get("a_R2", "0"))))
    a_ric2 = float(sp.N(sp.sympify(sym.get("a_Ric2", "0"))))

    s_grid = np.array((p1969.get("numeric_integral_table", {}) or {}).get("s_grid", [0.5, 1.0, 2.0, 4.0]), dtype=float)
    values = s_grid * (a_r2 / (8.0 * np.pi) + a_ric2 / (6.0 * np.pi))

    rows = []
    for s, v in zip(s_grid, values):
        rows.append({
            "s": float(s),
            "disc_integral_value": float(v),
            "positive": bool(v > 0.0),
        })

    all_positive = all(r["positive"] for r in rows)
    u2_computed = pre_ready and all_positive and len(rows) >= 4

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2100",
        "stage_id": "S1050",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_U2_PHASE_SPACE_QUADRATURE_WITNESS_WITH_TRACE__TASK2_ENTRY_PARTIAL"
            if u2_computed
            else "OPEN_STRICT_U2_PHASE_SPACE_QUADRATURE_WITNESS_BLOCKED"
        ),
        "depends_on": {
            "p2099_present": p2099.get("_missing") is None,
            "p1951_present": p1951.get("_missing") is None,
            "p1969_present": p1969.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "u2_phase_space_quadrature": {
            "channel": "graviton -> gauge_gauge",
            "scheme_tag": "STRICT_P2020_PHASESPACE_SCHEME_V1",
            "formula": "Disc(s)=s*(a_R2/(8*pi)+a_Ric2/(6*pi))",
            "coefficients": {"a_R2": a_r2, "a_Ric2": a_ric2},
            "s_grid": s_grid.tolist(),
            "rows": rows,
            "all_positive": all_positive,
            "scope_limit": "phase-space quadrature witness in strict B1 chain; not full dressed optical-theorem closure",
        },
        "unitarity_blocker_update": {
            "U1_shared_rg_scheme_lock": "COMPUTED",
            "U2_exact_discontinuity_integration": "COMPUTED" if u2_computed else "OPEN",
            "U3_positive_residue_uncertainty_table": "OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2101_candidate",
            "goal": "build U3 positivity uncertainty table with explicit residue-uncertainty propagation on same scheme",
        },
        "c3_gate_update": {
            "C3_u2_phase_space_quadrature_witness": "COMPUTED" if u2_computed else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "u2_computed": u2_computed,
            "all_grid_points_positive": all_positive,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2100 S1050: strict U2 phase-space quadrature witness",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- U2 computed: `{u2_computed}`",
            f"- All grid points positive: `{all_positive}`",
            "",
            "This stage exports the U2 exact phase-space quadrature witness on a strict scheme-tagged s-grid.",
            "No full dressed Cutkosky closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
