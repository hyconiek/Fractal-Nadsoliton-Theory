#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2123 = GEN / "p2123_s1073_strict_cmp2_map_choice_sensitivity_audit.json"
IN_2112 = GEN / "p2112_s1062_strict_cmp2_multibin_covariance_residual_table.json"
OUT = GEN / "p2125_s1075_strict_cmp2_extended_mapping_posterior_predictive_calibration_audit.json"
MD = GEN / "p2125_s1075_strict_cmp2_extended_mapping_posterior_predictive_calibration_audit.md"

SCHEMA_VERSION = "p2125_s1075_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2123 = load(IN_2123)
    p2112 = load(IN_2112)

    pre_ready = p2123.get("result_kind") == "PASS_STRICT_CMP2_MAP_CHOICE_SENSITIVITY_OBJECT_WITH_TRACE"

    map_rows = ((p2123.get("map_choice_sensitivity_object", {}) or {}).get("rows") or [])
    hist_rows = ((p2112.get("cmp2_multibin_residual_table", {}) or {}).get("rows") or [])
    n = min(len(map_rows), len(hist_rows))

    out_rows = []
    hits = 0

    for i in range(n):
        r = map_rows[i]
        h = hist_rows[i]

        d_nn = float(r.get("distance_nn", 0.0))
        d_m = float(r.get("distance_monotone", 0.0))

        # Extended mapping family (4 models):
        # M1 NN, M2 Monotone, M3 NN-tie-break conservative, M4 Monotone-penalized.
        # Distances are proxy energies.
        e = np.array([
            d_nn,
            d_m,
            d_nn + 0.01,   # tie-break penalty
            d_m + 0.02,    # monotone complexity penalty
        ], dtype=float)
        beta = 10.0
        w = np.exp(-beta * e)
        p = w / np.sum(w) if np.sum(w) > 0 else np.array([0.25, 0.25, 0.25, 0.25])

        u_nn = np.array(r.get("unnormalized_envelope_nn", [0.0, 0.0]), dtype=float)
        u_m = np.array(r.get("unnormalized_envelope_monotone", [0.0, 0.0]), dtype=float)
        n_nn = np.array(r.get("normalized_envelope_nn", [0.0, 0.0]), dtype=float)
        n_m = np.array(r.get("normalized_envelope_monotone", [0.0, 0.0]), dtype=float)

        # Construct family envelopes
        fam_u = np.stack([
            u_nn,
            u_m,
            np.array([u_nn[0] - 0.005, u_nn[1] + 0.005]),
            np.array([u_m[0] - 0.010, u_m[1] + 0.010]),
        ])
        fam_n = np.stack([
            n_nn,
            n_m,
            np.array([n_nn[0] - 0.005, n_nn[1] + 0.005]),
            np.array([n_m[0] - 0.010, n_m[1] + 0.010]),
        ])

        u_post = np.sum(fam_u * p[:, None], axis=0)
        n_post = np.sum(fam_n * p[:, None], axis=0)

        # Conservative q95 over discrete posterior family
        order = np.argsort(p)[::-1]
        cum = 0.0
        chosen = []
        for idx in order:
            chosen.append(idx)
            cum += float(p[idx])
            if cum >= 0.95:
                break
        u_q95 = [float(np.min(fam_u[chosen, 0])), float(np.max(fam_u[chosen, 1]))]
        n_q95 = [float(np.min(fam_n[chosen, 0])), float(np.max(fam_n[chosen, 1]))]

        residual_hist = float(h.get("signed_residual", 0.0))
        covered = (u_q95[0] <= residual_hist <= u_q95[1]) and (n_q95[0] <= residual_hist <= n_q95[1])
        hits += int(covered)

        out_rows.append(
            {
                "cmp_bin_index": int(r.get("cmp_bin_index", i)),
                "posterior_weights_extended_family": {
                    "M1_nn": float(p[0]),
                    "M2_monotone": float(p[1]),
                    "M3_nn_tiebreak": float(p[2]),
                    "M4_monotone_penalized": float(p[3]),
                },
                "unnormalized_posterior_envelope": [float(u_post[0]), float(u_post[1])],
                "normalized_posterior_envelope": [float(n_post[0]), float(n_post[1])],
                "unnormalized_conservative_q95_envelope": u_q95,
                "normalized_conservative_q95_envelope": n_q95,
                "historical_residual": residual_hist,
                "posterior_predictive_covered": covered,
            }
        )

    coverage_rate = float(hits / n) if n > 0 else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2125",
        "stage_id": "S1075",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_EXTENDED_MAPPING_POSTERIOR_PREDICTIVE_CALIBRATION_AUDIT_WITH_TRACE"
            if pre_ready and n > 0
            else "OPEN_STRICT_CMP2_EXTENDED_MAPPING_POSTERIOR_PREDICTIVE_CALIBRATION_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2123_present": p2123.get("_missing") is None,
            "p2112_present": p2112.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "extended_mapping_posterior_predictive_calibration_audit": {
            "rows": out_rows,
            "coverage_rate": coverage_rate,
            "posterior_family": ["M1_nn", "M2_monotone", "M3_nn_tiebreak", "M4_monotone_penalized"],
            "scope_limit": "operational posterior predictive calibration audit only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_EXTENDED_MAPPING_POSTERIOR_PREDICTIVE_CALIBRATION_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2126_candidate",
            "goal": "replace synthetic penalties in extended posterior family by backend-supported map evidence scores and rerun calibration",
        },
        "c3_gate_update": {
            "C3_cmp2_mapping_posterior_envelope_audit": "COMPUTED",
            "C3_cmp2_extended_mapping_posterior_predictive_calibration": "COMPUTED" if n > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "extended_posterior_rows_exported": n > 0,
            "posterior_predictive_calibration_computed": n > 0,
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
            "# P2125 S1075: strict CMP2 extended mapping posterior predictive calibration audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{n}`",
            f"- Coverage rate: `{coverage_rate}`",
            "",
            "This stage extends mapping posterior family beyond two maps and audits posterior predictive coverage on historical residual rows.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
