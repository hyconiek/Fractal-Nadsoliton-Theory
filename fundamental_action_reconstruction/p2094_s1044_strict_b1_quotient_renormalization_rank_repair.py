#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
INP = GEN / "p1950_s900_strict_renormalization_exact_integration_probe.json"
OUT = GEN / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.json"
MD = GEN / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.md"

SCHEMA_VERSION = "p2094_s1044_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1950 = load(INP)
    coeffs = p1950.get("aggregate_metrics", {})
    rows = [r for r in p1950.get("computed_rows", []) if "row_operator" in r]

    ready = len(rows) == 4 and all(k in coeffs for k in ["a_R2", "a_Ric2", "a_Riem2", "a_GB"])

    # Build reduced 3x3 quotient system by removing GB as independent channel.
    # Use rows/channels: R2, Ric2, Riem2 only.
    row_map = {r["channel"]: r for r in rows}
    channels = ["R2", "Ric2", "Riem2"]
    A3 = np.array(
        [[float(row_map[ch]["row_operator"][k]) for k in ["a_R2", "a_Ric2", "a_Riem2"]] for ch in channels],
        dtype=float,
    ) if ready else np.zeros((3, 3), dtype=float)
    b3 = np.array([float(row_map[ch]["rhs_divergence_target"]) for ch in channels], dtype=float) if ready else np.zeros(3, dtype=float)

    x3 = np.linalg.solve(A3, b3) if ready else np.zeros(3, dtype=float)
    residual3 = A3 @ x3 - b3
    residual3_l2 = float(np.linalg.norm(residual3, ord=2))
    residual3_max = float(np.max(np.abs(residual3))) if residual3.size else 0.0

    # Check GB row is reconstructed by topological identity combination in p1950 profile export.
    gb_row = row_map.get("GB", {})
    gb_lhs_from_x3 = float(
        gb_row.get("row_operator", {}).get("a_R2", 0.0) * x3[0]
        + gb_row.get("row_operator", {}).get("a_Ric2", 0.0) * x3[1]
        + gb_row.get("row_operator", {}).get("a_Riem2", 0.0) * x3[2]
    ) if ready else 0.0
    gb_rhs = float(gb_row.get("rhs_divergence_target", 0.0)) if ready else 0.0
    gb_residual = gb_lhs_from_x3 - gb_rhs

    tol = 1e-9
    quotient_pass = ready and residual3_max <= tol
    gb_consistent = ready and abs(gb_residual) <= 1e-6

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2094",
        "stage_id": "S1044",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_B1_QUOTIENT_RENORMALIZATION_RANK_REPAIR_WITH_TRACE__NO_FULL_GB_INDEPENDENCE"
            if quotient_pass and gb_consistent
            else "OPEN_STRICT_B1_QUOTIENT_RENORMALIZATION_RANK_REPAIR_BLOCKED"
        ),
        "depends_on": {
            "p1950_present": p1950.get("_missing") is None,
            "preconditions_ready": ready,
            "p1950_rank_defect_seen": p1950.get("fail_trace", "").startswith("operator_profile_gram_rank=3"),
        },
        "quotient_rank_repair": {
            "independent_channels": channels,
            "dropped_independent_channel": "GB",
            "rationale": "GB handled as topological/derived channel; not independent in current strict B1 projection matrix.",
            "A3": A3.tolist(),
            "b3": b3.tolist(),
            "solution": {
                "a_R2": float(x3[0]),
                "a_Ric2": float(x3[1]),
                "a_Riem2": float(x3[2]),
            },
            "residual_l2": residual3_l2,
            "residual_abs_max": residual3_max,
            "gb_row_reconstruction": {
                "lhs_from_3channel_solution": gb_lhs_from_x3,
                "rhs_target": gb_rhs,
                "signed_residual": gb_residual,
            },
        },
        "c3_gate_update": {
            "C3_strict_b1_quotient_renormalization_rank_repair": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "quotient_3x3_residual_small": quotient_pass,
            "gb_row_consistent_as_derived_channel": gb_consistent,
            "full_4channel_independence_proven": False,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2094 S1044: strict B1 quotient renormalization rank-repair",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- 3x3 quotient residual max: `{residual3_max}`",
            f"- GB derived-row residual: `{gb_residual}`",
            "",
            "This stage repairs the rank-defect path by quotienting out independent GB and solving strict 3-channel renormalization rows.",
            "No full 4-channel independence claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
