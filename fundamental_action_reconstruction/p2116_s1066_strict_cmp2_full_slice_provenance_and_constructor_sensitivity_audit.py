#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2115 = GEN / "p2115_s1065_strict_cmp2_full_binwise_channel_covariance_slices_and_coupled_robustness.json"
IN_2110 = GEN / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2116_s1066_strict_cmp2_full_slice_provenance_and_constructor_sensitivity_audit.json"
MD = GEN / "p2116_s1066_strict_cmp2_full_slice_provenance_and_constructor_sensitivity_audit.md"

SCHEMA_VERSION = "p2116_s1066_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
CHS = ["gg", "gh", "hh", "gx"]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v


def make_covariance(mats: list[np.ndarray], normalized: bool) -> np.ndarray:
    if normalized:
        mats_use = []
        for m in mats:
            fn = float(np.linalg.norm(m))
            mats_use.append(m / fn if fn > 0 else m)
    else:
        mats_use = mats

    c = np.zeros((4, 4), dtype=float)
    for a in range(4):
        for b in range(4):
            c[a, b] = float(np.sum(mats_use[a] * mats_use[b]))

    w, v = np.linalg.eigh(c)
    w_clip = np.maximum(w, 0.0)
    return (v * w_clip) @ v.T


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2115 = load(IN_2115)
    p2110 = load(IN_2110)
    p2017 = load(IN_2017)

    pre_ready = p2115.get("result_kind") == "PASS_STRICT_CMP2_FULL_BINWISE_CHANNEL_COVARIANCE_SLICES_AND_COUPLED_ROBUSTNESS_WITH_TRACE"

    cmp2 = next((r for r in ((p2110.get("contract_level_cmp_table", {}) or {}).get("rows") or []) if r.get("row_id") == "CMP2"), {})
    base_row_map = np.array(cmp2.get("row_map_to_channels", [0.0, 1.0, 0.0, 0.0]), dtype=float)
    if base_row_map.shape != (4,):
        base_row_map = np.array([0.0, 1.0, 0.0, 0.0], dtype=float)

    tensor_rows = p2017.get("tensor_candidate_table", []) or []
    prev_rows = ((p2115.get("cmp2_full_binwise_channel_covariance", {}) or {}).get("rows") or [])
    n_bins = min(len(tensor_rows), len(prev_rows))

    perturbations = [
        np.array([0.0, 0.0, 0.0, 0.0]),
        np.array([+0.01, -0.01, 0.0, 0.0]),
        np.array([-0.01, +0.01, 0.0, 0.0]),
    ]

    rows_out = []
    for i in range(n_bins):
        row_t = tensor_rows[i]
        cand = row_t.get("strict_quadrature_tensor_candidates", {}) or {}

        mats = []
        ok = True
        for ch in CHS:
            t = np.array(((cand.get(ch, {}) or {}).get("tensor_3x3", [])), dtype=float)
            if t.shape != (3, 3):
                ok = False
                break
            mats.append(t)
        if not ok:
            continue

        c_unnorm = make_covariance(mats, normalized=False)
        c_norm = make_covariance(mats, normalized=True)

        residual_center = float(prev_rows[i].get("interval_95_nominal", [0.0, 0.0])[0] + prev_rows[i].get("interval_95_nominal", [0.0, 0.0])[1]) / 2.0

        def envelope(cov: np.ndarray) -> tuple[float, float, float]:
            sigmas = []
            for p in perturbations:
                rmap = normalize(base_row_map + p)
                sigmas.append(float(np.sqrt(max(0.0, rmap @ cov @ rmap))))
            z95 = 1.96
            lo = [residual_center - z95 * s for s in sigmas]
            hi = [residual_center + z95 * s for s in sigmas]
            return float(min(lo)), float(max(hi)), sigmas[0]

        lo_u, hi_u, s_u = envelope(c_unnorm)
        lo_n, hi_n, s_n = envelope(c_norm)

        width_u = hi_u - lo_u
        width_n = hi_n - lo_n
        rel_width_delta = abs(width_u - width_n) / width_u if width_u > 0 else 0.0

        rows_out.append(
            {
                "bin_index": i,
                "s": row_t.get("s"),
                "provenance_map": {
                    "backend_tensor_candidate_row_index": i,
                    "backend_tensor_candidate_s": row_t.get("s"),
                    "channel_tuple": CHS,
                    "tensor_source_path": f"p2017::tensor_candidate_table[{i}].strict_quadrature_tensor_candidates",
                },
                "constructor_comparison": {
                    "unnormalized": {
                        "sigma_nominal": s_u,
                        "interval_95_envelope": [lo_u, hi_u],
                    },
                    "normalized": {
                        "sigma_nominal": s_n,
                        "interval_95_envelope": [lo_n, hi_n],
                    },
                    "relative_envelope_width_delta": rel_width_delta,
                },
                "constructor_stability_pass": rel_width_delta <= 0.25,
            }
        )

    pass_rate = (sum(1 for r in rows_out if r["constructor_stability_pass"]) / len(rows_out)) if rows_out else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2116",
        "stage_id": "S1066",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_FULL_SLICE_PROVENANCE_AND_CONSTRUCTOR_SENSITIVITY_AUDIT_WITH_TRACE"
            if pre_ready and len(rows_out) == n_bins and n_bins > 0
            else "OPEN_STRICT_CMP2_FULL_SLICE_PROVENANCE_AND_CONSTRUCTOR_SENSITIVITY_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2115_present": p2115.get("_missing") is None,
            "p2110_present": p2110.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp2_full_slice_provenance_and_constructor_sensitivity": {
            "rows": rows_out,
            "constructor_stability_pass_rate": pass_rate,
            "stability_threshold_relative_width_delta": 0.25,
            "scope_limit": "operational constructor-sensitivity audit only; not theorem-grade D3/C3 closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_FULL_SLICE_PROVENANCE_AND_CONSTRUCTOR_SENSITIVITY_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2117_candidate",
            "goal": "add explicit bin-alignment theorem object linking p2017 s-grid to CMP s-grid and rerun constructor sensitivity with alignment uncertainty",
        },
        "c3_gate_update": {
            "C3_cmp2_full_binwise_channel_covariance_slices": "COMPUTED",
            "C3_cmp2_full_slice_provenance_map": "COMPUTED" if len(rows_out) > 0 else "OPEN",
            "C3_cmp2_constructor_sensitivity_audit": "COMPUTED" if len(rows_out) > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "provenance_map_exported_all_rows": len(rows_out) == n_bins and n_bins > 0,
            "constructor_comparison_exported": len(rows_out) > 0,
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
            "# P2116 S1066: strict CMP2 full-slice provenance and constructor sensitivity audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(rows_out)}`",
            f"- Constructor stability pass rate: `{pass_rate}`",
            "",
            "This stage attaches explicit provenance maps for each full slice and audits sensitivity under alternative covariance constructors (normalized vs unnormalized inner-product).",
            "No global D3/C3 theorem or ToE closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
