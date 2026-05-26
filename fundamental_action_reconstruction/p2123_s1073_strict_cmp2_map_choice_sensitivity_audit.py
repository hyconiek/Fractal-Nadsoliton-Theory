#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2122 = GEN / "p2122_s1072_strict_cmp2_index_to_s_provenance_uncertainty_triple_intersection_audit.json"
IN_2115 = GEN / "p2115_s1065_strict_cmp2_full_binwise_channel_covariance_slices_and_coupled_robustness.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2123_s1073_strict_cmp2_map_choice_sensitivity_audit.json"
MD = GEN / "p2123_s1073_strict_cmp2_map_choice_sensitivity_audit.md"

SCHEMA_VERSION = "p2123_s1073_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def parse_s(v: Any) -> float:
    try:
        return float(eval(str(v), {"__builtins__": {}}, {}))
    except Exception:
        return float(v)


def nearest_neighbor_map(cmp_s: np.ndarray, backend_s: np.ndarray) -> list[int]:
    out = []
    for x in cmp_s:
        out.append(int(np.argmin(np.abs(backend_s - x))))
    return out


def monotone_map(cmp_s: np.ndarray, backend_s: np.ndarray) -> list[int]:
    out = []
    j = 0
    for x in cmp_s:
        while j + 1 < len(backend_s) and abs(backend_s[j + 1] - x) <= abs(backend_s[j] - x):
            j += 1
        out.append(j)
    return out


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2122 = load(IN_2122)
    p2115 = load(IN_2115)
    p2017 = load(IN_2017)

    pre_ready = p2122.get("result_kind") == "PASS_STRICT_CMP2_INDEX_TO_S_PROVENANCE_UNCERTAINTY_TRIPLE_INTERSECTION_AUDIT_WITH_TRACE"

    prev_rows = ((p2122.get("index_to_s_provenance_and_mapping_uncertainty", {}) or {}).get("rows") or [])
    cmp_rows = ((p2115.get("cmp2_full_binwise_channel_covariance", {}) or {}).get("rows") or [])
    backend_rows = p2017.get("tensor_candidate_table", []) or []

    n = min(len(prev_rows), len(cmp_rows), len(backend_rows))
    cmp_s = np.array([float(cmp_rows[i].get("s")) for i in range(n)], dtype=float)
    backend_s = np.array([parse_s(backend_rows[i].get("s", i)) for i in range(len(backend_rows))], dtype=float)

    nn = nearest_neighbor_map(cmp_s, backend_s)
    mono = monotone_map(cmp_s, backend_s)

    rows_out = []
    for i in range(n):
        base_u = prev_rows[i].get("unnormalized_triple_intersection_with_mapping_uncertainty", [0.0, 0.0])
        base_n = prev_rows[i].get("normalized_triple_intersection_with_mapping_uncertainty", [0.0, 0.0])

        j_nn = nn[i]
        j_m = mono[i]

        d_nn = abs(cmp_s[i] - backend_s[j_nn])
        d_m = abs(cmp_s[i] - backend_s[j_m])

        s_nn = 0.5 * d_nn + 1e-6
        s_m = 0.5 * d_m + 1e-6

        # map-dependent propagated envelopes
        u_nn = [float(base_u[0] - 1.96 * s_nn), float(base_u[1] + 1.96 * s_nn)]
        u_m = [float(base_u[0] - 1.96 * s_m), float(base_u[1] + 1.96 * s_m)]
        n_nn = [float(base_n[0] - 1.96 * s_nn), float(base_n[1] + 1.96 * s_nn)]
        n_m = [float(base_n[0] - 1.96 * s_m), float(base_n[1] + 1.96 * s_m)]

        rows_out.append(
            {
                "cmp_bin_index": i,
                "cmp_s": float(cmp_s[i]),
                "nearest_neighbor_backend_index": j_nn,
                "nearest_neighbor_backend_s": float(backend_s[j_nn]),
                "monotone_backend_index": j_m,
                "monotone_backend_s": float(backend_s[j_m]),
                "distance_nn": float(d_nn),
                "distance_monotone": float(d_m),
                "unnormalized_envelope_nn": u_nn,
                "unnormalized_envelope_monotone": u_m,
                "normalized_envelope_nn": n_nn,
                "normalized_envelope_monotone": n_m,
                "map_choice_abs_width_diff_unnormalized": abs((u_nn[1] - u_nn[0]) - (u_m[1] - u_m[0])),
                "map_choice_abs_width_diff_normalized": abs((n_nn[1] - n_nn[0]) - (n_m[1] - n_m[0])),
            }
        )

    width_diffs = [
        max(r["map_choice_abs_width_diff_unnormalized"], r["map_choice_abs_width_diff_normalized"]) for r in rows_out
    ]
    max_diff = float(max(width_diffs)) if width_diffs else 0.0
    mean_diff = float(np.mean(width_diffs)) if width_diffs else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2123",
        "stage_id": "S1073",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_MAP_CHOICE_SENSITIVITY_OBJECT_WITH_TRACE"
            if pre_ready and len(rows_out) > 0
            else "OPEN_STRICT_CMP2_MAP_CHOICE_SENSITIVITY_OBJECT_BLOCKED"
        ),
        "depends_on": {
            "p2122_present": p2122.get("_missing") is None,
            "p2115_present": p2115.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "map_choice_sensitivity_object": {
            "mapping_methods": ["nearest_neighbor", "monotone"],
            "rows": rows_out,
            "summary": {
                "max_abs_width_diff": max_diff,
                "mean_abs_width_diff": mean_diff,
            },
            "scope_limit": "operational map-choice sensitivity only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_MAP_CHOICE_SENSITIVITY_OBJECT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2124_candidate",
            "goal": "introduce probabilistic map posterior over candidate alignments and integrate envelopes under map uncertainty",
        },
        "c3_gate_update": {
            "C3_cmp2_index_to_s_provenance_object": "COMPUTED",
            "C3_cmp2_map_choice_sensitivity_object": "COMPUTED" if len(rows_out) > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "map_choice_rows_exported": len(rows_out) > 0,
            "map_choice_sensitivity_computed": len(rows_out) > 0,
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
            "# P2123 S1073: strict CMP2 map-choice sensitivity object",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(rows_out)}`",
            f"- Max width diff (map choice): `{max_diff}`",
            "",
            "This stage compares nearest-neighbor vs monotone CMP↔backend mapping and exports propagated-envelope sensitivity to map choice.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
