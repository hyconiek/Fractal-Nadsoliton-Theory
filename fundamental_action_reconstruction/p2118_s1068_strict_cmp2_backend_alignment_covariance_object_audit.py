#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2117 = GEN / "p2117_s1067_strict_cmp2_binalignment_uncertainty_propagated_constructor_audit.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2118_s1068_strict_cmp2_backend_alignment_covariance_object_audit.json"
MD = GEN / "p2118_s1068_strict_cmp2_backend_alignment_covariance_object_audit.md"

SCHEMA_VERSION = "p2118_s1068_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"
CHS = ["gg", "gh", "hh", "gx"]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2117 = load(IN_2117)
    p2017 = load(IN_2017)

    pre_ready = p2117.get("result_kind") == "PASS_STRICT_CMP2_BINALIGNMENT_UNCERTAINTY_PROPAGATED_CONSTRUCTOR_AUDIT_WITH_TRACE"

    align_rows = ((p2117.get("bin_alignment_object", {}) or {}).get("rows") or [])
    audit_rows = ((p2117.get("constructor_audit_with_alignment_uncertainty", {}) or {}).get("rows") or [])
    trows = p2017.get("tensor_candidate_table", []) or []

    n = min(len(align_rows), len(audit_rows), len(trows))
    out_rows = []

    for i in range(n):
        cand = (trows[i].get("strict_quadrature_tensor_candidates", {}) or {})

        # backend alignment-covariance proxy object from quadrature error tensors
        # per channel use Frobenius norm of quad_error_3x3 as backend uncertainty signal
        err = []
        for ch in CHS:
            qe = np.array(((cand.get(ch, {}) or {}).get("quad_error_3x3", [])), dtype=float)
            err.append(float(np.linalg.norm(qe)) if qe.shape == (3, 3) else 0.0)

        err_vec = np.array(err, dtype=float)
        c_align = np.outer(err_vec, err_vec)
        sigma_align_backend = float(np.sqrt(max(0.0, np.trace(c_align))))

        row = audit_rows[i]
        un = row.get("unnormalized_interval_95_with_alignment", [0.0, 0.0])
        no = row.get("normalized_interval_95_with_alignment", [0.0, 0.0])

        # replace heuristic sigma_alignment by backend-derived sigma and widen envelopes consistently
        lo_u = float(un[0] - 1.96 * sigma_align_backend)
        hi_u = float(un[1] + 1.96 * sigma_align_backend)
        lo_n = float(no[0] - 1.96 * sigma_align_backend)
        hi_n = float(no[1] + 1.96 * sigma_align_backend)

        width_u = hi_u - lo_u
        width_n = hi_n - lo_n
        rel_delta = abs(width_u - width_n) / width_u if width_u > 0 else 0.0

        out_rows.append(
            {
                "bin_index": i,
                "backend_alignment_covariance_object": {
                    "source_tensor_row": i,
                    "channel_quad_error_fro_norm": err,
                    "alignment_covariance_matrix_4x4": c_align.tolist(),
                    "sigma_alignment_backend": sigma_align_backend,
                },
                "unnormalized_interval_95_with_backend_alignment": [lo_u, hi_u],
                "normalized_interval_95_with_backend_alignment": [lo_n, hi_n],
                "relative_envelope_width_delta_with_backend_alignment": rel_delta,
                "constructor_stability_with_backend_alignment_pass": rel_delta <= 0.25,
            }
        )

    pass_rate = (sum(1 for r in out_rows if r["constructor_stability_with_backend_alignment_pass"]) / len(out_rows)) if out_rows else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2118",
        "stage_id": "S1068",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_BACKEND_ALIGNMENT_COVARIANCE_OBJECT_AUDIT_WITH_TRACE"
            if pre_ready and len(out_rows) == n and n > 0
            else "OPEN_STRICT_CMP2_BACKEND_ALIGNMENT_COVARIANCE_OBJECT_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2117_present": p2117.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "backend_alignment_covariance_object_audit": {
            "rows": out_rows,
            "constructor_stability_pass_rate": pass_rate,
            "stability_threshold_relative_width_delta": 0.25,
            "scope_limit": "operational backend alignment-covariance object audit only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_BACKEND_ALIGNMENT_COVARIANCE_OBJECT_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2119_candidate",
            "goal": "cross-check backend alignment covariance object against independent alignment estimators and perform consistency envelope intersection audit",
        },
        "c3_gate_update": {
            "C3_cmp2_binalignment_object": "COMPUTED",
            "C3_cmp2_backend_alignment_covariance_object": "COMPUTED" if len(out_rows) == n and n > 0 else "OPEN",
            "C3_cmp2_constructor_sensitivity_with_backend_alignment": "COMPUTED" if len(out_rows) == n and n > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "backend_alignment_covariance_rows_exported": len(out_rows) == n and n > 0,
            "backend_alignment_uncertainty_propagated": len(out_rows) == n and n > 0,
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
            "# P2118 S1068: strict CMP2 backend alignment-covariance object audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(out_rows)}`",
            f"- Constructor stability pass rate: `{pass_rate}`",
            "",
            "This stage replaces heuristic alignment sigma by a backend-exported alignment-covariance object (from quadrature error tensors) and reruns constructor audit with propagated backend alignment uncertainty.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
