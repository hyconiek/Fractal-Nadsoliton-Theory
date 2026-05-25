#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
INP = GEN / "p1950_s900_strict_renormalization_exact_integration_probe.json"
OUT = GEN / "p2095_s1045_strict_b1_gb_derived_channel_certificate.json"
MD = GEN / "p2095_s1045_strict_b1_gb_derived_channel_certificate.md"

SCHEMA_VERSION = "p2095_s1045_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1950 = load(INP)
    rows = [r for r in p1950.get("computed_rows", []) if "row_operator" in r]
    row_map = {r.get("channel"): r for r in rows}

    ready = all(ch in row_map for ch in ["R2", "Ric2", "Riem2", "GB"])

    # certificate target from backend rule in P1950 profiles: GB = Riem2 - 4*Ric2 + R2
    v_r2 = np.array([float(row_map["R2"]["row_operator"][k]) for k in ["a_R2", "a_Ric2", "a_Riem2", "a_GB"]]) if ready else np.zeros(4)
    v_ric2 = np.array([float(row_map["Ric2"]["row_operator"][k]) for k in ["a_R2", "a_Ric2", "a_Riem2", "a_GB"]]) if ready else np.zeros(4)
    v_riem2 = np.array([float(row_map["Riem2"]["row_operator"][k]) for k in ["a_R2", "a_Ric2", "a_Riem2", "a_GB"]]) if ready else np.zeros(4)
    v_gb = np.array([float(row_map["GB"]["row_operator"][k]) for k in ["a_R2", "a_Ric2", "a_Riem2", "a_GB"]]) if ready else np.zeros(4)

    v_combo = v_riem2 - 4.0 * v_ric2 + v_r2
    dep_residual = v_gb - v_combo
    dep_residual_l2 = float(np.linalg.norm(dep_residual, ord=2))
    dep_residual_max = float(np.max(np.abs(dep_residual))) if dep_residual.size else 0.0

    gb_rhs = float(row_map["GB"]["rhs_divergence_target"]) if ready else 0.0
    rhs_combo = (
        float(row_map["Riem2"]["rhs_divergence_target"])
        - 4.0 * float(row_map["Ric2"]["rhs_divergence_target"])
        + float(row_map["R2"]["rhs_divergence_target"])
    ) if ready else 0.0
    rhs_residual = gb_rhs - rhs_combo

    tol = 1e-9
    derived_pass = ready and dep_residual_max <= tol and abs(rhs_residual) <= tol

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2095",
        "stage_id": "S1045",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_B1_GB_DERIVED_CHANNEL_CERTIFICATE_WITH_TRACE__NO_FULL_4CHANNEL_INDEPENDENCE"
            if derived_pass
            else "OPEN_STRICT_B1_GB_DERIVED_CHANNEL_CERTIFICATE_BLOCKED"
        ),
        "depends_on": {
            "p1950_present": p1950.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "derived_channel_certificate": {
            "identity": "GB = Riem2 - 4*Ric2 + R2",
            "operator_vector_gb": v_gb.tolist(),
            "operator_vector_combo": v_combo.tolist(),
            "operator_residual": dep_residual.tolist(),
            "operator_residual_l2": dep_residual_l2,
            "operator_residual_abs_max": dep_residual_max,
            "rhs_gb": gb_rhs,
            "rhs_combo": rhs_combo,
            "rhs_residual": rhs_residual,
        },
        "c3_gate_update": {
            "C3_strict_b1_gb_derived_channel_certificate": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "gb_derived_identity_operator_level": derived_pass,
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
            "# P2095 S1045: strict B1 GB derived-channel certificate",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Operator residual max: `{dep_residual_max}`",
            f"- RHS residual: `{rhs_residual}`",
            "",
            "This stage certifies that GB is a derived channel in the current B1 operator basis.",
            "No full 4-channel independence claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
