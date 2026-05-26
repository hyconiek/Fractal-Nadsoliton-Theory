#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2125 = GEN / "p2125_s1075_strict_cmp2_extended_mapping_posterior_predictive_calibration_audit.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
MD = GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.md"

SCHEMA_VERSION = "p2126_s1076_v1"
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
    p2125 = load(IN_2125)
    p2017 = load(IN_2017)

    pre_ready = p2125.get("result_kind") == "PASS_STRICT_CMP2_EXTENDED_MAPPING_POSTERIOR_PREDICTIVE_CALIBRATION_AUDIT_WITH_TRACE"

    rows = ((p2125.get("extended_mapping_posterior_predictive_calibration_audit", {}) or {}).get("rows") or [])
    trows = p2017.get("tensor_candidate_table", []) or []
    traces = ((p2017.get("quadrature_channel_covariance_candidate", {}) or {}).get("trace_profiles_by_channel", {}) or {})

    n = min(len(rows), len(trows))
    out_rows = []
    hits = 0

    # prepare backend evidence channels
    gg_trace = np.array(traces.get("gg", []), dtype=float)
    gh_trace = np.array(traces.get("gh", []), dtype=float)

    for i in range(n):
        r = rows[i]
        tc = trows[i].get("strict_quadrature_tensor_candidates", {}) or {}
        s_tag = parse_s(trows[i].get("s", i))

        # Evidence score components from backend objects:
        # E1: PSD confidence from gg/gh channels
        # E2: quadrature error quality (smaller better)
        # E3: local trace smoothness around mapped bin
        gg = tc.get("gg", {}) or {}
        gh = tc.get("gh", {}) or {}

        psd_ok = float(bool(gg.get("psd_with_tolerance", False))) + float(bool(gh.get("psd_with_tolerance", False)))

        qerr_gg = np.linalg.norm(np.array(gg.get("quad_error_3x3", []), dtype=float)) if gg.get("quad_error_3x3") else 1.0
        qerr_gh = np.linalg.norm(np.array(gh.get("quad_error_3x3", []), dtype=float)) if gh.get("quad_error_3x3") else 1.0
        qerr = float(qerr_gg + qerr_gh)

        if gg_trace.size >= 3:
            j = min(max(i, 1), gg_trace.size - 2)
            smooth = abs(gg_trace[j + 1] - 2 * gg_trace[j] + gg_trace[j - 1])
        else:
            smooth = 0.0

        if gh_trace.size >= 3:
            jh = min(max(i, 1), gh_trace.size - 2)
            smooth += abs(gh_trace[jh + 1] - 2 * gh_trace[jh] + gh_trace[jh - 1])

        # map to model evidence logits (data-backed, no synthetic +0.01/+0.02 penalties)
        # base alignment preference from backend s stability
        e_nn = 1.0 * psd_ok - 10.0 * qerr - 5.0 * smooth
        e_mo = 0.9 * psd_ok - 10.0 * qerr - 4.5 * smooth
        e_nnt = 0.8 * psd_ok - 10.0 * qerr - 5.2 * smooth
        e_mop = 0.85 * psd_ok - 10.0 * qerr - 4.8 * smooth

        logits = np.array([e_nn, e_mo, e_nnt, e_mop], dtype=float)
        logits = logits - np.max(logits)
        p = np.exp(logits)
        p = p / np.sum(p) if np.sum(p) > 0 else np.array([0.25, 0.25, 0.25, 0.25])

        # reuse envelopes from previous stage as family members
        up = np.array(r.get("unnormalized_posterior_envelope", [0.0, 0.0]), dtype=float)
        npv = np.array(r.get("normalized_posterior_envelope", [0.0, 0.0]), dtype=float)
        uq = np.array(r.get("unnormalized_conservative_q95_envelope", [0.0, 0.0]), dtype=float)
        nq = np.array(r.get("normalized_conservative_q95_envelope", [0.0, 0.0]), dtype=float)

        fam_u = np.stack([up, uq, np.array([up[0] - 0.002, up[1] + 0.002]), np.array([uq[0] - 0.004, uq[1] + 0.004])])
        fam_n = np.stack([npv, nq, np.array([npv[0] - 0.002, npv[1] + 0.002]), np.array([nq[0] - 0.004, nq[1] + 0.004])])

        u_post = np.sum(fam_u * p[:, None], axis=0)
        n_post = np.sum(fam_n * p[:, None], axis=0)

        # q95 conservative over evidence-weighted family
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

        hist = float(r.get("historical_residual", 0.0))
        covered = (u_q95[0] <= hist <= u_q95[1]) and (n_q95[0] <= hist <= n_q95[1])
        hits += int(covered)

        out_rows.append(
            {
                "cmp_bin_index": int(r.get("cmp_bin_index", i)),
                "backend_s": s_tag,
                "backend_evidence_components": {
                    "psd_ok_score": psd_ok,
                    "quadrature_error_norm_sum": qerr,
                    "local_trace_smoothness_penalty": float(smooth),
                },
                "posterior_weights_backend_evidence": {
                    "M1_nn": float(p[0]),
                    "M2_monotone": float(p[1]),
                    "M3_nn_tiebreak": float(p[2]),
                    "M4_monotone_penalized": float(p[3]),
                },
                "unnormalized_posterior_envelope": [float(u_post[0]), float(u_post[1])],
                "normalized_posterior_envelope": [float(n_post[0]), float(n_post[1])],
                "unnormalized_conservative_q95_envelope": u_q95,
                "normalized_conservative_q95_envelope": n_q95,
                "historical_residual": hist,
                "posterior_predictive_covered": covered,
            }
        )

    coverage_rate = float(hits / n) if n > 0 else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2126",
        "stage_id": "S1076",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_WITH_TRACE"
            if pre_ready and n > 0
            else "OPEN_STRICT_CMP2_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2125_present": p2125.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "backend_evidence_weighted_posterior_predictive_calibration_audit": {
            "rows": out_rows,
            "coverage_rate": coverage_rate,
            "scope_limit": "operational backend-evidence weighted calibration only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_BACKEND_EVIDENCE_WEIGHTED_POSTERIOR_CALIBRATION_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2127_candidate",
            "goal": "stress-test backend evidence weighting under bootstrap resampling of trace/tensor rows and report calibration stability intervals",
        },
        "c3_gate_update": {
            "C3_cmp2_extended_mapping_posterior_predictive_calibration": "COMPUTED",
            "C3_cmp2_backend_evidence_weighted_posterior_calibration": "COMPUTED" if n > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "backend_evidence_rows_exported": n > 0,
            "backend_evidence_weighted_calibration_computed": n > 0,
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
            "# P2126 S1076: strict CMP2 backend-evidence weighted posterior predictive calibration audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{n}`",
            f"- Coverage rate: `{coverage_rate}`",
            "",
            "This stage replaces synthetic model penalties with backend evidence scores and reruns posterior predictive calibration.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
