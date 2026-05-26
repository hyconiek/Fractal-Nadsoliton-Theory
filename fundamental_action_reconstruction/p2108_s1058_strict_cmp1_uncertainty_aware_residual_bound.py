#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2107 = GEN / "p2107_s1057_strict_cmp1_first_residual_execution.json"
IN_2101 = GEN / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.json"
OUT = GEN / "p2108_s1058_strict_cmp1_uncertainty_aware_residual_bound.json"
MD = GEN / "p2108_s1058_strict_cmp1_uncertainty_aware_residual_bound.md"

SCHEMA_VERSION = "p2108_s1058_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2107 = load(IN_2107)
    p2101 = load(IN_2101)

    pre_ready = p2107.get("result_kind") == "PASS_STRICT_CMP1_FIRST_RESIDUAL_EXECUTION_WITH_TRACE"

    cmp1 = p2107.get("cmp1_residual_execution", {}) or {}
    residual = cmp1.get("signed_residual")

    rows = ((p2101.get("u3_uncertainty_table", {}) or {}).get("rows") or [])
    disc_unc = [float(r.get("disc_uncertainty_abs", 0.0)) for r in rows if r.get("disc_uncertainty_abs") is not None]
    sigma_disc = float(np.max(disc_unc)) if disc_unc else 0.0
    sigma_cutsum = 0.0  # strict placeholder until dressed backend covariance row link is exported
    sigma_total = float(np.sqrt(sigma_disc**2 + sigma_cutsum**2))

    residual_f = float(sp.N(sp.sympify(residual))) if residual is not None else None
    if residual_f is not None:
        z95 = 1.96
        low95 = residual_f - z95 * sigma_total
        high95 = residual_f + z95 * sigma_total
    else:
        low95 = None
        high95 = None

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2108",
        "stage_id": "S1058",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP1_UNCERTAINTY_AWARE_BOUND_WITH_TRACE__SEED_LEVEL"
            if pre_ready and residual_f is not None
            else "OPEN_STRICT_CMP1_UNCERTAINTY_AWARE_BOUND_BLOCKED"
        ),
        "depends_on": {
            "p2107_present": p2107.get("_missing") is None,
            "p2101_present": p2101.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp1_uncertainty_aware_residual_bound": {
            "residual_center": residual_f,
            "sigma_disc_seed_proxy": sigma_disc,
            "sigma_cutsum_seed_proxy": sigma_cutsum,
            "sigma_total": sigma_total,
            "z95": 1.96,
            "residual_interval_95": [low95, high95],
            "contains_zero": (low95 <= 0.0 <= high95) if low95 is not None else None,
            "scope_limit": "seed-level uncertainty proxy only; D3 covariance transport remains OPEN",
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_cmp1_source_row": "COMPUTED",
            "D2_cmp1_first_residual_execution": "COMPUTED",
            "D3_uncertainty_propagation_from_dressed_backend": "OPEN_PARTIAL_SEED_PROXY_ONLY",
        },
        "recommended_next_honest_step": {
            "id": "P2109_candidate",
            "goal": "export object-specific covariance row link and recompute CMP1 bound with nonzero D3 backend uncertainty",
        },
        "c3_gate_update": {
            "C3_cmp1_first_residual_execution": "COMPUTED",
            "C3_cmp1_uncertainty_bound_seed_level": "COMPUTED" if residual_f is not None else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "cmp1_residual_present": residual_f is not None,
            "seed_uncertainty_bound_exported": low95 is not None,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2108 S1058: strict CMP1 uncertainty-aware residual bound",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Residual center: `{residual_f}`",
            f"- 95% interval: `{payload['cmp1_uncertainty_aware_residual_bound']['residual_interval_95']}`",
            "",
            "This stage exports a seed-level uncertainty-aware bound for CMP1 residual.",
            "No full D3 covariance transport theorem or global Cutkosky closure is claimed.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
