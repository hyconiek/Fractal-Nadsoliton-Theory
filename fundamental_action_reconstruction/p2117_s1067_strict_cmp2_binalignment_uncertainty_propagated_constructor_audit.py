#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2116 = GEN / "p2116_s1066_strict_cmp2_full_slice_provenance_and_constructor_sensitivity_audit.json"
IN_2115 = GEN / "p2115_s1065_strict_cmp2_full_binwise_channel_covariance_slices_and_coupled_robustness.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2117_s1067_strict_cmp2_binalignment_uncertainty_propagated_constructor_audit.json"
MD = GEN / "p2117_s1067_strict_cmp2_binalignment_uncertainty_propagated_constructor_audit.md"

SCHEMA_VERSION = "p2117_s1067_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2116 = load(IN_2116)
    p2115 = load(IN_2115)
    p2017 = load(IN_2017)

    pre_ready = p2116.get("result_kind") == "PASS_STRICT_CMP2_FULL_SLICE_PROVENANCE_AND_CONSTRUCTOR_SENSITIVITY_AUDIT_WITH_TRACE"

    rows_2116 = ((p2116.get("cmp2_full_slice_provenance_and_constructor_sensitivity", {}) or {}).get("rows") or [])
    rows_2115 = ((p2115.get("cmp2_full_binwise_channel_covariance", {}) or {}).get("rows") or [])
    tensor_rows = p2017.get("tensor_candidate_table", []) or []

    n = min(len(rows_2116), len(rows_2115), len(tensor_rows))
    align_rows = []
    out_rows = []

    for i in range(n):
        s_cmp = float(rows_2115[i].get("s")) if rows_2115[i].get("s") is not None else float(i)
        s_bk_raw = tensor_rows[i].get("s", i)
        try:
            s_bk = float(eval(str(s_bk_raw), {"__builtins__": {}}, {}))
        except Exception:
            s_bk = float(i)

        delta_s = abs(s_cmp - s_bk)
        sigma_align = 0.5 * delta_s + 1e-6
        align_rows.append(
            {
                "cmp_bin_index": i,
                "cmp_s": s_cmp,
                "backend_row_index": i,
                "backend_s_raw": str(s_bk_raw),
                "backend_s_numeric": s_bk,
                "abs_delta_s": delta_s,
                "sigma_alignment": sigma_align,
            }
        )

        c = rows_2116[i].get("constructor_comparison", {})
        un = c.get("unnormalized", {}).get("interval_95_envelope", [0.0, 0.0])
        no = c.get("normalized", {}).get("interval_95_envelope", [0.0, 0.0])

        # propagate bin-alignment uncertainty as additive envelope widening
        lo_u = float(un[0] - 1.96 * sigma_align)
        hi_u = float(un[1] + 1.96 * sigma_align)
        lo_n = float(no[0] - 1.96 * sigma_align)
        hi_n = float(no[1] + 1.96 * sigma_align)

        width_u = hi_u - lo_u
        width_n = hi_n - lo_n
        rel_delta = abs(width_u - width_n) / width_u if width_u > 0 else 0.0

        out_rows.append(
            {
                "bin_index": i,
                "alignment_object_ref": f"bin_alignment_rows[{i}]",
                "sigma_alignment": sigma_align,
                "unnormalized_interval_95_with_alignment": [lo_u, hi_u],
                "normalized_interval_95_with_alignment": [lo_n, hi_n],
                "relative_envelope_width_delta_with_alignment": rel_delta,
                "constructor_stability_with_alignment_pass": rel_delta <= 0.25,
            }
        )

    pass_rate = (sum(1 for r in out_rows if r["constructor_stability_with_alignment_pass"]) / len(out_rows)) if out_rows else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2117",
        "stage_id": "S1067",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_BINALIGNMENT_UNCERTAINTY_PROPAGATED_CONSTRUCTOR_AUDIT_WITH_TRACE"
            if pre_ready and len(out_rows) == n and n > 0
            else "OPEN_STRICT_CMP2_BINALIGNMENT_UNCERTAINTY_PROPAGATED_CONSTRUCTOR_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2116_present": p2116.get("_missing") is None,
            "p2115_present": p2115.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "bin_alignment_object": {
            "source_cmp": "p2115::cmp2_full_binwise_channel_covariance.rows",
            "source_backend": "p2017::tensor_candidate_table",
            "rows": align_rows,
            "scope_limit": "operational bin alignment map with uncertainty; not theorem-grade alignment proof",
        },
        "constructor_audit_with_alignment_uncertainty": {
            "rows": out_rows,
            "constructor_stability_pass_rate": pass_rate,
            "stability_threshold_relative_width_delta": 0.25,
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_BINALIGNMENT_UNCERTAINTY_PROPAGATED_CONSTRUCTOR_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2118_candidate",
            "goal": "replace heuristic sigma_alignment by backend-exported alignment covariance and rerun propagated constructor audit",
        },
        "c3_gate_update": {
            "C3_cmp2_full_slice_provenance_map": "COMPUTED",
            "C3_cmp2_binalignment_object": "COMPUTED" if len(align_rows) == n and n > 0 else "OPEN",
            "C3_cmp2_constructor_sensitivity_with_alignment_uncertainty": "COMPUTED" if len(out_rows) == n and n > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "binalignment_rows_exported": len(align_rows) == n and n > 0,
            "alignment_uncertainty_propagated": len(out_rows) == n and n > 0,
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
            "# P2117 S1067: strict CMP2 bin-alignment uncertainty propagated constructor audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Alignment rows: `{len(align_rows)}`",
            f"- Constructor stability pass rate: `{pass_rate}`",
            "",
            "This stage adds a formal bin-alignment object (CMP grid vs p2017 grid) with alignment uncertainty and propagates it into constructor-audit 95% envelopes.",
            "No global theorem-grade C3/D3 closure or ToE closure is claimed.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
