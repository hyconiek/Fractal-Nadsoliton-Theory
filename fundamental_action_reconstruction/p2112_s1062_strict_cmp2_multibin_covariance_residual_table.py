#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2111 = GEN / "p2111_s1061_strict_cmp2_covariance_row_export_and_residual_execution.json"
IN_1862 = GEN / "p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json"
IN_2101 = GEN / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.json"
OUT = GEN / "p2112_s1062_strict_cmp2_multibin_covariance_residual_table.json"
MD = GEN / "p2112_s1062_strict_cmp2_multibin_covariance_residual_table.md"

SCHEMA_VERSION = "p2112_s1062_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2111 = load(IN_2111)
    p1862 = load(IN_1862)
    p2101 = load(IN_2101)

    pre_ready = p2111.get("result_kind") == "PASS_STRICT_CMP2_NUMERIC_COVARIANCE_ROW_EXPORT_AND_FIRST_RESIDUAL_WITH_TRACE"

    disc_rows = ((p1862.get("projected_discontinuity_seed_evaluation", {}) or {}).get("rows") or [])
    u3_rows = ((p2101.get("u3_uncertainty_table", {}) or {}).get("rows") or [])

    s_bins = sorted({float(r.get("s")) for r in disc_rows if r.get("s") is not None})

    cmp2_rows = []
    for s in s_bins:
        vals = [float(r.get("disc_projected_seed")) for r in disc_rows if float(r.get("s")) == s and r.get("disc_projected_seed") is not None]
        dressed_cov = float(np.var(vals, ddof=0)) if vals else None

        matching_u3 = [r for r in u3_rows if r.get("s") is not None and float(r.get("s")) == s]
        cutsum_unc = float(sp.N(sp.sympify(matching_u3[0].get("disc_uncertainty_abs", 0.0)))) if matching_u3 else None

        sigma_proxy = float(sp.N(sp.sympify(matching_u3[0].get("disc_uncertainty_abs", 0.0)))) if matching_u3 else 0.0
        residual = (dressed_cov - cutsum_unc) if (dressed_cov is not None and cutsum_unc is not None) else None
        z95 = 1.96
        interval = [residual - z95 * sigma_proxy, residual + z95 * sigma_proxy] if residual is not None else [None, None]

        cmp2_rows.append(
            {
                "s": s,
                "dressed_covariance_row_value": dressed_cov,
                "cutsum_uncertainty_row_value": cutsum_unc,
                "signed_residual": residual,
                "abs_residual": abs(residual) if residual is not None else None,
                "sigma_proxy": sigma_proxy,
                "residual_interval_95": interval,
                "contains_zero": (interval[0] <= 0.0 <= interval[1]) if residual is not None else None,
            }
        )

    residuals = [r["signed_residual"] for r in cmp2_rows if r["signed_residual"] is not None]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2112",
        "stage_id": "S1062",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_MULTIBIN_RESIDUAL_TABLE_WITH_TRACE"
            if pre_ready and len(cmp2_rows) > 0 and all(r["signed_residual"] is not None for r in cmp2_rows)
            else "OPEN_STRICT_CMP2_MULTIBIN_RESIDUAL_TABLE_BLOCKED"
        ),
        "depends_on": {
            "p2111_present": p2111.get("_missing") is None,
            "p1862_present": p1862.get("_missing") is None,
            "p2101_present": p2101.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp2_multibin_residual_table": {
            "row_id": "CMP2",
            "scheme_tag": "STRICT_P2020_PHASESPACE_SCHEME_V1",
            "rows": cmp2_rows,
            "residual_summary": {
                "n_bins": len(cmp2_rows),
                "mean_signed_residual": float(np.mean(residuals)) if residuals else None,
                "max_abs_residual": float(np.max(np.abs(residuals))) if residuals else None,
            },
            "scope_limit": "multibin operational table only; not global covariance transport theorem",
        },
        "backend_import_blocker_update": {
            "D1_dressed_pole_residue_backend_on_common_basis": "COMPUTED",
            "D2_cmp1_first_residual_execution": "COMPUTED",
            "D2_cmp2_residual_execution": "COMPUTED_MULTIBIN_TABLE",
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_MULTIBIN_PROXY_TABLE_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2113_candidate",
            "goal": "replace CMP2 sigma_proxy by explicit dressed-backend covariance transport per s-bin and re-evaluate interval table",
        },
        "c3_gate_update": {
            "C3_cmp2_first_residual_execution": "COMPUTED",
            "C3_cmp2_multibin_residual_table": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "cmp2_multibin_rows_exported": len(cmp2_rows) > 0,
            "cmp2_multibin_residuals_computed": all(r["signed_residual"] is not None for r in cmp2_rows) if cmp2_rows else False,
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
            "# P2112 S1062: strict CMP2 multibin covariance residual table",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Bins exported: `{len(cmp2_rows)}`",
            f"- Mean residual: `{payload['cmp2_multibin_residual_table']['residual_summary']['mean_signed_residual']}`",
            "",
            "This stage executes CMP2 residual rows across all available s-bins with explicit rowwise intervals.",
            "No global covariance transport theorem, full Cutkosky closure, or ToE closure is claimed.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
