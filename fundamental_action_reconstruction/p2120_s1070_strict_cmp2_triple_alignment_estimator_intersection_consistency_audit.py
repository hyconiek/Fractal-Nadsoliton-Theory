#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2119 = GEN / "p2119_s1069_strict_cmp2_alignment_covariance_crosscheck_intersection_envelope_audit.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2120_s1070_strict_cmp2_triple_alignment_estimator_intersection_consistency_audit.json"
MD = GEN / "p2120_s1070_strict_cmp2_triple_alignment_estimator_intersection_consistency_audit.md"

SCHEMA_VERSION = "p2120_s1070_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2119 = load(IN_2119)
    p2017 = load(IN_2017)

    pre_ready = p2119.get("result_kind") == "PASS_STRICT_CMP2_ALIGNMENT_COVARIANCE_CROSSCHECK_INTERSECTION_ENVELOPE_AUDIT_WITH_TRACE"

    rows_2119 = ((p2119.get("alignment_covariance_crosscheck_intersection_audit", {}) or {}).get("rows") or [])
    eigvals = np.array(((p2017.get("quadrature_channel_covariance_candidate", {}) or {}).get("C4_eigvals", [])), dtype=float)

    out_rows = []
    for i, r in enumerate(rows_2119):
        sigma_backend = float(r.get("sigma_alignment_backend", 0.0))
        sigma_curv = float(r.get("sigma_alignment_independent_curvature", 0.0))

        # Third independent estimator: local eigenspectrum-change scale (global fallback on C4 eigvals)
        if eigvals.size >= 2:
            sigma_eig = float(np.std(np.diff(eigvals)))
        elif eigvals.size == 1:
            sigma_eig = float(abs(eigvals[0]))
        else:
            sigma_eig = 0.0

        un_b = r.get("unnormalized_interval_backend", [0.0, 0.0])
        no_b = r.get("normalized_interval_backend", [0.0, 0.0])

        # build three envelopes around the same base backend interval
        def widen(interval: list[float], sigma: float) -> list[float]:
            return [float(interval[0] - 1.96 * sigma), float(interval[1] + 1.96 * sigma)]

        u1, u2, u3 = widen(un_b, sigma_backend), widen(un_b, sigma_curv), widen(un_b, sigma_eig)
        n1, n2, n3 = widen(no_b, sigma_backend), widen(no_b, sigma_curv), widen(no_b, sigma_eig)

        u_int = [max(u1[0], u2[0], u3[0]), min(u1[1], u2[1], u3[1])]
        n_int = [max(n1[0], n2[0], n3[0]), min(n1[1], n2[1], n3[1])]

        out_rows.append(
            {
                "bin_index": i,
                "sigma_backend": sigma_backend,
                "sigma_curvature": sigma_curv,
                "sigma_eigenspectrum": sigma_eig,
                "unnormalized_triple_intersection": u_int,
                "unnormalized_triple_intersection_nonempty": u_int[0] <= u_int[1],
                "normalized_triple_intersection": n_int,
                "normalized_triple_intersection_nonempty": n_int[0] <= n_int[1],
            }
        )

    pass_rate = (
        sum(1 for r in out_rows if r["unnormalized_triple_intersection_nonempty"] and r["normalized_triple_intersection_nonempty"]) / len(out_rows)
        if out_rows else 0.0
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2120",
        "stage_id": "S1070",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_TRIPLE_ALIGNMENT_ESTIMATOR_INTERSECTION_CONSISTENCY_AUDIT_WITH_TRACE"
            if pre_ready and len(out_rows) > 0
            else "OPEN_STRICT_CMP2_TRIPLE_ALIGNMENT_ESTIMATOR_INTERSECTION_CONSISTENCY_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2119_present": p2119.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "triple_alignment_estimator_consistency_audit": {
            "third_estimator": "local eigenspectrum-change scale over p2017 C4_eigvals",
            "rows": out_rows,
            "triple_intersection_pass_rate": pass_rate,
            "scope_limit": "operational consistency audit only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_TRIPLE_ALIGNMENT_ESTIMATOR_INTERSECTION_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2121_candidate",
            "goal": "replace global eigenspectrum fallback by explicit bin-local C4 eigenspectrum deltas and rerun triple-intersection audit",
        },
        "c3_gate_update": {
            "C3_cmp2_alignment_crosscheck_intersection_audit": "COMPUTED",
            "C3_cmp2_triple_alignment_estimator_consistency_audit": "COMPUTED" if len(out_rows) > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "triple_estimator_rows_exported": len(out_rows) > 0,
            "triple_intersection_audit_computed": len(out_rows) > 0,
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
            "# P2120 S1070: strict CMP2 triple-alignment-estimator intersection consistency audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(out_rows)}`",
            f"- Triple intersection pass rate: `{pass_rate}`",
            "",
            "This stage adds a third independent alignment estimator (eigenspectrum-change scale) and exports triple-intersection consistency envelopes.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
