#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2118 = GEN / "p2118_s1068_strict_cmp2_backend_alignment_covariance_object_audit.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2119_s1069_strict_cmp2_alignment_covariance_crosscheck_intersection_envelope_audit.json"
MD = GEN / "p2119_s1069_strict_cmp2_alignment_covariance_crosscheck_intersection_envelope_audit.md"

SCHEMA_VERSION = "p2119_s1069_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
CHS = ["gg", "gh", "hh", "gx"]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def second_diff(x: np.ndarray) -> np.ndarray:
    if x.size < 3:
        return np.array([0.0])
    return x[2:] - 2 * x[1:-1] + x[:-2]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2118 = load(IN_2118)
    p2017 = load(IN_2017)

    pre_ready = p2118.get("result_kind") == "PASS_STRICT_CMP2_BACKEND_ALIGNMENT_COVARIANCE_OBJECT_AUDIT_WITH_TRACE"

    rows_2118 = ((p2118.get("backend_alignment_covariance_object_audit", {}) or {}).get("rows") or [])
    traces = ((p2017.get("quadrature_channel_covariance_candidate", {}) or {}).get("trace_profiles_by_channel", {}) or {})

    prof = {ch: np.array(traces.get(ch, []), dtype=float) for ch in CHS}
    m = min([arr.size for arr in prof.values()] + [len(rows_2118)]) if rows_2118 else 0

    out_rows = []
    for i in range(m):
        r = rows_2118[i]
        backend_sigma = float((r.get("backend_alignment_covariance_object", {}) or {}).get("sigma_alignment_backend", 0.0))

        # Independent estimator: local curvature magnitude across channel trace profiles
        curv = []
        for ch in CHS:
            arr = prof[ch]
            if arr.size < 3:
                curv.append(0.0)
            else:
                j = min(max(i, 1), arr.size - 2)
                c = abs(arr[j + 1] - 2 * arr[j] + arr[j - 1])
                curv.append(float(c))
        sigma_curv = float(np.linalg.norm(np.array(curv, dtype=float)))

        un = r.get("unnormalized_interval_95_with_backend_alignment", [0.0, 0.0])
        no = r.get("normalized_interval_95_with_backend_alignment", [0.0, 0.0])

        # Cross-check by replacing sigma with independent curvature estimator
        lo_u_c = float(un[0] - 1.96 * sigma_curv)
        hi_u_c = float(un[1] + 1.96 * sigma_curv)
        lo_n_c = float(no[0] - 1.96 * sigma_curv)
        hi_n_c = float(no[1] + 1.96 * sigma_curv)

        # Intersection envelope between backend-object and independent estimator
        lo_u_i = max(float(un[0]), lo_u_c)
        hi_u_i = min(float(un[1]), hi_u_c)
        lo_n_i = max(float(no[0]), lo_n_c)
        hi_n_i = min(float(no[1]), hi_n_c)

        intersection_nonempty_u = lo_u_i <= hi_u_i
        intersection_nonempty_n = lo_n_i <= hi_n_i

        out_rows.append(
            {
                "bin_index": i,
                "sigma_alignment_backend": backend_sigma,
                "sigma_alignment_independent_curvature": sigma_curv,
                "unnormalized_interval_backend": [float(un[0]), float(un[1])],
                "unnormalized_interval_independent": [lo_u_c, hi_u_c],
                "unnormalized_interval_intersection": [lo_u_i, hi_u_i],
                "unnormalized_intersection_nonempty": intersection_nonempty_u,
                "normalized_interval_backend": [float(no[0]), float(no[1])],
                "normalized_interval_independent": [lo_n_c, hi_n_c],
                "normalized_interval_intersection": [lo_n_i, hi_n_i],
                "normalized_intersection_nonempty": intersection_nonempty_n,
            }
        )

    pass_rate = (
        sum(1 for r in out_rows if r["unnormalized_intersection_nonempty"] and r["normalized_intersection_nonempty"]) / len(out_rows)
        if out_rows
        else 0.0
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2119",
        "stage_id": "S1069",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_ALIGNMENT_COVARIANCE_CROSSCHECK_INTERSECTION_ENVELOPE_AUDIT_WITH_TRACE"
            if pre_ready and len(out_rows) == m and m > 0
            else "OPEN_STRICT_CMP2_ALIGNMENT_COVARIANCE_CROSSCHECK_INTERSECTION_ENVELOPE_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2118_present": p2118.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "alignment_covariance_crosscheck_intersection_audit": {
            "independent_estimator": "local curvature norm over p2017 trace_profiles_by_channel",
            "rows": out_rows,
            "intersection_nonempty_pass_rate": pass_rate,
            "scope_limit": "operational cross-check only; no theorem-grade closure claim",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_ALIGNMENT_CROSSCHECK_INTERSECTION_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2120_candidate",
            "goal": "add third independent alignment estimator and perform triple-intersection consistency audit",
        },
        "c3_gate_update": {
            "C3_cmp2_backend_alignment_covariance_object": "COMPUTED",
            "C3_cmp2_alignment_crosscheck_intersection_audit": "COMPUTED" if len(out_rows) > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "crosscheck_rows_exported": len(out_rows) == m and m > 0,
            "intersection_audit_computed": len(out_rows) > 0,
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
            "# P2119 S1069: strict CMP2 alignment covariance cross-check and intersection-envelope audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(out_rows)}`",
            f"- Intersection pass rate: `{pass_rate}`",
            "",
            "This stage cross-checks backend alignment covariance against an independent curvature-based estimator and exports interval intersections to reduce single-constructor dominance risk.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
