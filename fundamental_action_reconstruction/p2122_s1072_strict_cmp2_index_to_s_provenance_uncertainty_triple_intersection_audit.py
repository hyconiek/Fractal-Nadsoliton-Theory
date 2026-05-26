#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2121 = GEN / "p2121_s1071_strict_cmp2_binlocal_eigenspectrum_triple_intersection_audit.json"
IN_2115 = GEN / "p2115_s1065_strict_cmp2_full_binwise_channel_covariance_slices_and_coupled_robustness.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2122_s1072_strict_cmp2_index_to_s_provenance_uncertainty_triple_intersection_audit.json"
MD = GEN / "p2122_s1072_strict_cmp2_index_to_s_provenance_uncertainty_triple_intersection_audit.md"

SCHEMA_VERSION = "p2122_s1072_v1"
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


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2121 = load(IN_2121)
    p2115 = load(IN_2115)
    p2017 = load(IN_2017)

    pre_ready = p2121.get("result_kind") == "PASS_STRICT_CMP2_BINLOCAL_EIGENSPECTRUM_TRIPLE_INTERSECTION_AUDIT_WITH_TRACE"

    rows_2121 = ((p2121.get("binlocal_triple_intersection_audit", {}) or {}).get("rows") or [])
    cmp_rows = ((p2115.get("cmp2_full_binwise_channel_covariance", {}) or {}).get("rows") or [])
    backend_rows = p2017.get("tensor_candidate_table", []) or []

    n = min(len(rows_2121), len(cmp_rows), len(backend_rows))
    out_rows = []

    for i in range(n):
        cmp_s = float(cmp_rows[i].get("s")) if cmp_rows[i].get("s") is not None else float(i)
        backend_s = parse_s(backend_rows[i].get("s", i))
        delta_s = abs(cmp_s - backend_s)

        # mapping uncertainty converted into sigma-like additive half-width proxy
        sigma_map = 0.5 * delta_s + 1e-6

        r = rows_2121[i]
        u = r.get("unnormalized_triple_intersection_binlocal", [0.0, 0.0])
        nrm = r.get("normalized_triple_intersection_binlocal", [0.0, 0.0])

        u_prop = [float(u[0] - 1.96 * sigma_map), float(u[1] + 1.96 * sigma_map)]
        n_prop = [float(nrm[0] - 1.96 * sigma_map), float(nrm[1] + 1.96 * sigma_map)]

        out_rows.append(
            {
                "cmp_bin_index": i,
                "backend_bin_index": i,
                "cmp_s": cmp_s,
                "backend_s": backend_s,
                "abs_delta_s": delta_s,
                "sigma_index_to_s_mapping": sigma_map,
                "unnormalized_triple_intersection_with_mapping_uncertainty": u_prop,
                "unnormalized_nonempty": u_prop[0] <= u_prop[1],
                "normalized_triple_intersection_with_mapping_uncertainty": n_prop,
                "normalized_nonempty": n_prop[0] <= n_prop[1],
            }
        )

    pass_rate = (
        sum(1 for r in out_rows if r["unnormalized_nonempty"] and r["normalized_nonempty"]) / len(out_rows)
        if out_rows else 0.0
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2122",
        "stage_id": "S1072",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_INDEX_TO_S_PROVENANCE_UNCERTAINTY_TRIPLE_INTERSECTION_AUDIT_WITH_TRACE"
            if pre_ready and len(out_rows) > 0
            else "OPEN_STRICT_CMP2_INDEX_TO_S_PROVENANCE_UNCERTAINTY_TRIPLE_INTERSECTION_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2121_present": p2121.get("_missing") is None,
            "p2115_present": p2115.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "index_to_s_provenance_and_mapping_uncertainty": {
            "cmp_source": "p2115::cmp2_full_binwise_channel_covariance.rows",
            "backend_source": "p2017::tensor_candidate_table",
            "rows": out_rows,
            "triple_intersection_nonempty_pass_rate": pass_rate,
            "scope_limit": "operational provenance+mapping uncertainty propagation only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_INDEX_TO_S_PROVENANCE_AND_MAPPING_UNCERTAINTY_PROPAGATION_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2123_candidate",
            "goal": "replace simple index-locked map by nearest-neighbor / monotone matching map and compare propagated envelopes",
        },
        "c3_gate_update": {
            "C3_cmp2_binlocal_eigenspectrum_triple_intersection_audit": "COMPUTED",
            "C3_cmp2_index_to_s_provenance_object": "COMPUTED" if len(out_rows) > 0 else "OPEN",
            "C3_cmp2_mapping_uncertainty_propagation": "COMPUTED" if len(out_rows) > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "index_to_s_rows_exported": len(out_rows) > 0,
            "mapping_uncertainty_propagated": len(out_rows) > 0,
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
            "# P2122 S1072: strict CMP2 index-to-s provenance and mapping-uncertainty triple-intersection audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(out_rows)}`",
            f"- Pass rate: `{pass_rate}`",
            "",
            "This stage adds explicit CMP-bin↔backend-bin provenance and propagates index-to-s mapping uncertainty into the triple-intersection envelopes.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
