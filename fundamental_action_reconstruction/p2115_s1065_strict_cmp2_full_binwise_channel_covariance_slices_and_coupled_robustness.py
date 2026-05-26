#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2114 = GEN / "p2114_s1064_strict_cmp2_backend_binwise_covariance_and_rowmap_robustness_audit.json"
IN_2110 = GEN / "p2110_s1060_strict_contract_level_cmp_table_d3_rowmap_calibration.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2115_s1065_strict_cmp2_full_binwise_channel_covariance_slices_and_coupled_robustness.json"
MD = GEN / "p2115_s1065_strict_cmp2_full_binwise_channel_covariance_slices_and_coupled_robustness.md"

SCHEMA_VERSION = "p2115_s1065_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
CHS = ["gg", "gh", "hh", "gx"]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2114 = load(IN_2114)
    p2110 = load(IN_2110)
    p2017 = load(IN_2017)

    pre_ready = p2114.get("result_kind") == "PASS_STRICT_CMP2_BACKEND_BINWISE_COVARIANCE_AND_ROWMAP_ROBUSTNESS_AUDIT_WITH_TRACE"

    cmp2 = next((r for r in ((p2110.get("contract_level_cmp_table", {}) or {}).get("rows") or []) if r.get("row_id") == "CMP2"), {})
    base_row_map = np.array(cmp2.get("row_map_to_channels", [0.0, 1.0, 0.0, 0.0]), dtype=float)
    if base_row_map.shape != (4,):
        base_row_map = np.array([0.0, 1.0, 0.0, 0.0], dtype=float)

    tensor_rows = p2017.get("tensor_candidate_table", []) or []
    prev_rows = ((p2114.get("cmp2_rowmap_robustness_audit", {}) or {}).get("rows") or [])

    n_bins = min(len(tensor_rows), len(prev_rows))
    perturbations = [
        np.array([0.0, 0.0, 0.0, 0.0]),
        np.array([+0.01, -0.01, 0.0, 0.0]),
        np.array([-0.01, +0.01, 0.0, 0.0]),
        np.array([0.0, 0.0, +0.01, -0.01]),
        np.array([0.0, 0.0, -0.01, +0.01]),
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

        # full coupled channel covariance slice from Frobenius inner products of channel tensors
        c_full = np.zeros((4, 4), dtype=float)
        for a in range(4):
            for b in range(4):
                c_full[a, b] = float(np.sum(mats[a] * mats[b]))

        # PSD cleanup by eigenvalue clipping (numerical hygiene only)
        w, v = np.linalg.eigh(c_full)
        w_clip = np.maximum(w, 0.0)
        c_psd = (v * w_clip) @ v.T

        sigmas = []
        for p in perturbations:
            rmap = normalize(base_row_map + p)
            sigmas.append(float(np.sqrt(max(0.0, rmap @ c_psd @ rmap))))

        residual_center = float(prev_rows[i].get("interval_95_nominal", [0.0, 0.0])[0] + prev_rows[i].get("interval_95_nominal", [0.0, 0.0])[1]) / 2.0
        z95 = 1.96
        lo = [residual_center - z95 * s for s in sigmas]
        hi = [residual_center + z95 * s for s in sigmas]

        rows_out.append(
            {
                "bin_index": i,
                "s": prev_rows[i].get("s"),
                "full_channel_covariance_slice": c_psd.tolist(),
                "offdiag_l1": float(np.sum(np.abs(c_psd - np.diag(np.diag(c_psd))))),
                "sigma_nominal": sigmas[0],
                "sigma_min_under_rowmap_perturbation": float(min(sigmas)),
                "sigma_max_under_rowmap_perturbation": float(max(sigmas)),
                "interval_95_nominal": [lo[0], hi[0]],
                "interval_95_envelope": [float(min(lo)), float(max(hi))],
            }
        )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2115",
        "stage_id": "S1065",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_FULL_BINWISE_CHANNEL_COVARIANCE_SLICES_AND_COUPLED_ROBUSTNESS_WITH_TRACE"
            if pre_ready and n_bins > 0 and len(rows_out) == n_bins
            else "OPEN_STRICT_CMP2_FULL_BINWISE_CHANNEL_COVARIANCE_SLICES_AND_COUPLED_ROBUSTNESS_BLOCKED"
        ),
        "depends_on": {
            "p2114_present": p2114.get("_missing") is None,
            "p2110_present": p2110.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "cmp2_full_binwise_channel_covariance": {
            "channels": CHS,
            "construction_rule": "Frobenius inner-product covariance from p2017 per-bin channel 3x3 tensors",
            "rows": rows_out,
            "scope_limit": "operational full-slice export and coupled robustness only; not theorem-grade D3/C3 closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_FULL_BINWISE_CHANNEL_COVARIANCE_SLICES_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2116_candidate",
            "goal": "bind full binwise covariance slices to explicit backend provenance map and perform sensitivity audit vs alternative slice constructors",
        },
        "c3_gate_update": {
            "C3_cmp2_rowmap_robustness_audit": "COMPUTED",
            "C3_cmp2_full_binwise_channel_covariance_slices": "COMPUTED" if len(rows_out) == n_bins and n_bins > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "full_binwise_slices_exported": len(rows_out) == n_bins and n_bins > 0,
            "offdiagonal_couplings_present": any(r.get("offdiag_l1", 0.0) > 0.0 for r in rows_out),
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
            "# P2115 S1065: strict CMP2 full binwise channel covariance slices + coupled robustness",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Bin rows exported: `{len(rows_out)}`",
            f"- Couplings present: `{payload['gatekeeper_checks']['offdiagonal_couplings_present']}`",
            "",
            "This stage exports full (non-diagonal) binwise channel covariance slices and recomputes robustness envelopes with inter-channel couplings.",
            "No global D3/C3 theorem or ToE closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
