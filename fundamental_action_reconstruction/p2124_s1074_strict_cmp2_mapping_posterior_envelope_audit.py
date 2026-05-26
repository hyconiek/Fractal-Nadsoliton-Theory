#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2123 = GEN / "p2123_s1073_strict_cmp2_map_choice_sensitivity_audit.json"
OUT = GEN / "p2124_s1074_strict_cmp2_mapping_posterior_envelope_audit.json"
MD = GEN / "p2124_s1074_strict_cmp2_mapping_posterior_envelope_audit.md"

SCHEMA_VERSION = "p2124_s1074_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2123 = load(IN_2123)

    pre_ready = p2123.get("result_kind") == "PASS_STRICT_CMP2_MAP_CHOICE_SENSITIVITY_OBJECT_WITH_TRACE"
    rows = ((p2123.get("map_choice_sensitivity_object", {}) or {}).get("rows") or [])

    out_rows = []
    for row in rows:
        d_nn = float(row.get("distance_nn", 0.0))
        d_m = float(row.get("distance_monotone", 0.0))

        # posterior over map choices: higher weight for smaller mapping distance
        beta = 10.0
        w_nn = np.exp(-beta * d_nn)
        w_m = np.exp(-beta * d_m)
        z = w_nn + w_m
        p_nn = float(w_nn / z) if z > 0 else 0.5
        p_m = float(w_m / z) if z > 0 else 0.5

        u_nn = row.get("unnormalized_envelope_nn", [0.0, 0.0])
        u_m = row.get("unnormalized_envelope_monotone", [0.0, 0.0])
        n_nn = row.get("normalized_envelope_nn", [0.0, 0.0])
        n_m = row.get("normalized_envelope_monotone", [0.0, 0.0])

        # posterior-averaged envelope
        u_post = [p_nn * float(u_nn[0]) + p_m * float(u_m[0]), p_nn * float(u_nn[1]) + p_m * float(u_m[1])]
        n_post = [p_nn * float(n_nn[0]) + p_m * float(n_m[0]), p_nn * float(n_nn[1]) + p_m * float(n_m[1])]

        # conservative envelope via 95% quantile over two-map posterior mixture
        if p_nn >= 0.95:
            u_q95 = [float(u_nn[0]), float(u_nn[1])]
            n_q95 = [float(n_nn[0]), float(n_nn[1])]
        elif p_m >= 0.95:
            u_q95 = [float(u_m[0]), float(u_m[1])]
            n_q95 = [float(n_m[0]), float(n_m[1])]
        else:
            u_q95 = [min(float(u_nn[0]), float(u_m[0])), max(float(u_nn[1]), float(u_m[1]))]
            n_q95 = [min(float(n_nn[0]), float(n_m[0])), max(float(n_nn[1]), float(n_m[1]))]

        out_rows.append(
            {
                "cmp_bin_index": int(row.get("cmp_bin_index", 0)),
                "posterior_map_weights": {"nearest_neighbor": p_nn, "monotone": p_m},
                "unnormalized_posterior_envelope": u_post,
                "normalized_posterior_envelope": n_post,
                "unnormalized_conservative_q95_envelope": u_q95,
                "normalized_conservative_q95_envelope": n_q95,
            }
        )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2124",
        "stage_id": "S1074",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_MAPPING_POSTERIOR_AND_CONSERVATIVE_ENVELOPE_AUDIT_WITH_TRACE"
            if pre_ready and len(out_rows) > 0
            else "OPEN_STRICT_CMP2_MAPPING_POSTERIOR_AND_CONSERVATIVE_ENVELOPE_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2123_present": p2123.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "mapping_posterior_envelope_audit": {
            "rows": out_rows,
            "posterior_rule": "softmax(-beta*distance), beta=10",
            "conservative_rule": "95% map-quantile envelope over two-map mixture",
            "scope_limit": "operational posterior map uncertainty audit only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_MAPPING_POSTERIOR_AND_CONSERVATIVE_ENVELOPE_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2125_candidate",
            "goal": "add richer mapping family posterior (beyond two maps) and perform posterior predictive calibration audit",
        },
        "c3_gate_update": {
            "C3_cmp2_map_choice_sensitivity_object": "COMPUTED",
            "C3_cmp2_mapping_posterior_envelope_audit": "COMPUTED" if len(out_rows) > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "posterior_rows_exported": len(out_rows) > 0,
            "conservative_q95_exported": len(out_rows) > 0,
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
            "# P2124 S1074: strict CMP2 mapping posterior + conservative envelope audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(out_rows)}`",
            "",
            "This stage builds a posterior over mapping choices and exports posterior-averaged plus conservative (q95) envelopes.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
