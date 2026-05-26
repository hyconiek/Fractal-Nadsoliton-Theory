#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2113 = GEN / "p2113_s1063_strict_cmp2_binwise_d3_covariance_transport_table.json"
IN_2110 = GEN / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.json"
IN_2009 = GEN / "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2114_s1064_strict_cmp2_backend_binwise_covariance_and_rowmap_robustness_audit.json"
MD = GEN / "p2114_s1064_strict_cmp2_backend_binwise_covariance_and_rowmap_robustness_audit.md"

SCHEMA_VERSION = "p2114_s1064_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2113 = load(IN_2113)
    p2110 = load(IN_2110)
    p2009 = load(IN_2009)
    p2017 = load(IN_2017)

    pre_ready = (
        p2113.get("result_kind") == "PASS_STRICT_CMP2_BINWISE_D3_COVARIANCE_TRANSPORT_TABLE_WITH_TRACE"
        and p2110.get("result_kind") == "PASS_STRICT_CONTRACT_LEVEL_CMP_TABLE_WITH_CALIBRATED_ROW_MAPS_AND_D3_LINK_TRACE"
    )

    c4 = np.array((p2009.get("backend_tensor_objects", {}) or {}).get("C4_channel_covariance", []), dtype=float)
    c4_valid = c4.shape == (4, 4)

    cmp2 = next((r for r in ((p2110.get("contract_level_cmp_table", {}) or {}).get("rows") or []) if r.get("row_id") == "CMP2"), {})
    base_row_map = np.array(cmp2.get("row_map_to_channels", [0.0, 1.0, 0.0, 0.0]), dtype=float)
    if base_row_map.shape != (4,):
        base_row_map = np.array([0.0, 1.0, 0.0, 0.0], dtype=float)

    traces = ((p2017.get("quadrature_channel_covariance_candidate", {}) or {}).get("trace_profiles_by_channel", {}) or {})
    channels = ["gg", "gh", "hh", "gx"]
    profile_arrays = [np.array(traces.get(ch, []), dtype=float) for ch in channels]
    profile_lengths = [arr.size for arr in profile_arrays]
    has_backend_binwise_object = all(n > 0 for n in profile_lengths)

    rows_prev = ((p2113.get("cmp2_binwise_d3_transport_table", {}) or {}).get("rows") or [])
    n_bins = len(rows_prev)

    backend_rows = []
    if has_backend_binwise_object and c4_valid and n_bins > 0:
        m = min(profile_lengths)
        idxs = np.linspace(0, m - 1, n_bins)

        # small perturbations around calibrated row map for robustness audit
        perturbations = [
            np.array([0.0, 0.0, 0.0, 0.0]),
            np.array([+0.01, -0.01, 0.0, 0.0]),
            np.array([-0.01, +0.01, 0.0, 0.0]),
            np.array([0.0, 0.0, +0.01, -0.01]),
            np.array([0.0, 0.0, -0.01, +0.01]),
        ]

        for i, row in enumerate(rows_prev):
            j = int(round(float(idxs[i])))
            j = max(0, min(m - 1, j))
            # diagonal backend bin-wise covariance slice from exported channel traces
            diag = np.array([max(0.0, float(arr[j])) for arr in profile_arrays], dtype=float)
            c_bin = np.diag(diag)

            sigmas = []
            for p in perturbations:
                rmap = normalize(base_row_map + p)
                sigma = float(np.sqrt(max(0.0, rmap @ c_bin @ rmap)))
                sigmas.append(sigma)

            residual_center = float(row.get("signed_residual", 0.0))
            z95 = 1.96
            interval_lo = [residual_center - z95 * s for s in sigmas]
            interval_hi = [residual_center + z95 * s for s in sigmas]

            backend_rows.append(
                {
                    "s": float(row.get("s")),
                    "backend_profile_index": j,
                    "backend_binwise_covariance_diag": diag.tolist(),
                    "sigma_nominal": sigmas[0],
                    "sigma_min_under_rowmap_perturbation": float(min(sigmas)),
                    "sigma_max_under_rowmap_perturbation": float(max(sigmas)),
                    "interval_95_nominal": [interval_lo[0], interval_hi[0]],
                    "interval_95_envelope": [float(min(interval_lo)), float(max(interval_hi))],
                }
            )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2114",
        "stage_id": "S1064",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_BACKEND_BINWISE_COVARIANCE_AND_ROWMAP_ROBUSTNESS_AUDIT_WITH_TRACE"
            if pre_ready and c4_valid and has_backend_binwise_object and len(backend_rows) == n_bins and n_bins > 0
            else "OPEN_STRICT_CMP2_BACKEND_BINWISE_COVARIANCE_AND_ROWMAP_ROBUSTNESS_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2113_present": p2113.get("_missing") is None,
            "p2110_present": p2110.get("_missing") is None,
            "p2009_present": p2009.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "backend_binwise_covariance_source": {
            "source": "p2017::quadrature_channel_covariance_candidate.trace_profiles_by_channel",
            "channels": channels,
            "profile_lengths": profile_lengths,
            "backend_binwise_object_available": has_backend_binwise_object,
            "mapping_note": "s-bin rows mapped to nearest backend profile index",
        },
        "cmp2_rowmap_robustness_audit": {
            "base_row_map": base_row_map.tolist(),
            "perturbation_radius_l2_approx": 0.01414213562373095,
            "rows": backend_rows,
            "scope_limit": "operational robustness audit only; not theorem-grade global transport closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": (
                "COMPUTED_BACKEND_BINWISE_OBJECT_WITH_ROWMAP_ROBUSTNESS_AUDIT_PARTIAL"
                if has_backend_binwise_object and len(backend_rows) == n_bins and n_bins > 0
                else "OPEN"
            )
        },
        "recommended_next_honest_step": {
            "id": "P2115_candidate",
            "goal": "export explicit non-diagonal binwise channel covariance slices and rerun robustness envelope with coupled covariance terms",
        },
        "c3_gate_update": {
            "C3_cmp2_binwise_d3_transport_table": "COMPUTED" if n_bins > 0 else "OPEN",
            "C3_cmp2_rowmap_robustness_audit": "COMPUTED" if len(backend_rows) == n_bins and n_bins > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "backend_binwise_object_available": has_backend_binwise_object,
            "robustness_rows_exported": len(backend_rows) == n_bins and n_bins > 0,
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
            "# P2114 S1064: strict CMP2 backend binwise covariance + row-map robustness audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Backend binwise object available: `{has_backend_binwise_object}`",
            f"- Exported robustness rows: `{len(backend_rows)}`",
            "",
            "This stage replaces spread-based transport factor with backend-exported binwise covariance traces (when available) and audits interval stability under small row-map perturbations.",
            "No global C3/D3 transport theorem or ToE closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
