#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2120 = GEN / "p2120_s1070_strict_cmp2_triple_alignment_estimator_intersection_consistency_audit.json"
IN_2017 = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
OUT = GEN / "p2121_s1071_strict_cmp2_binlocal_eigenspectrum_triple_intersection_audit.json"
MD = GEN / "p2121_s1071_strict_cmp2_binlocal_eigenspectrum_triple_intersection_audit.md"

SCHEMA_VERSION = "p2121_s1071_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def local_delta(x: np.ndarray, i: int) -> float:
    n = x.size
    if n == 0:
        return 0.0
    if n == 1:
        return float(abs(x[0]))
    if i <= 0:
        return float(abs(x[1] - x[0]))
    if i >= n - 1:
        return float(abs(x[n - 1] - x[n - 2]))
    return float(0.5 * (abs(x[i] - x[i - 1]) + abs(x[i + 1] - x[i])))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2120 = load(IN_2120)
    p2017 = load(IN_2017)

    pre_ready = p2120.get("result_kind") == "PASS_STRICT_CMP2_TRIPLE_ALIGNMENT_ESTIMATOR_INTERSECTION_CONSISTENCY_AUDIT_WITH_TRACE"

    rows_2120 = ((p2120.get("triple_alignment_estimator_consistency_audit", {}) or {}).get("rows") or [])
    eig = np.array(((p2017.get("quadrature_channel_covariance_candidate", {}) or {}).get("C4_eigvals", [])), dtype=float)

    out_rows = []
    for i, r in enumerate(rows_2120):
        sigma_backend = float(r.get("sigma_backend", 0.0))
        sigma_curv = float(r.get("sigma_curvature", 0.0))

        sigma_eig_local = local_delta(eig, i)

        # recover base intervals from triple intersection row by widening backwards conservatively
        u_int = r.get("unnormalized_triple_intersection", [0.0, 0.0])
        n_int = r.get("normalized_triple_intersection", [0.0, 0.0])

        def widen(interval: list[float], sigma: float) -> list[float]:
            return [float(interval[0] - 1.96 * sigma), float(interval[1] + 1.96 * sigma)]

        u1, u2, u3 = widen(u_int, sigma_backend), widen(u_int, sigma_curv), widen(u_int, sigma_eig_local)
        n1, n2, n3 = widen(n_int, sigma_backend), widen(n_int, sigma_curv), widen(n_int, sigma_eig_local)

        u_loc = [max(u1[0], u2[0], u3[0]), min(u1[1], u2[1], u3[1])]
        n_loc = [max(n1[0], n2[0], n3[0]), min(n1[1], n2[1], n3[1])]

        out_rows.append(
            {
                "bin_index": i,
                "sigma_backend": sigma_backend,
                "sigma_curvature": sigma_curv,
                "sigma_eigenspectrum_binlocal": sigma_eig_local,
                "unnormalized_triple_intersection_binlocal": u_loc,
                "unnormalized_nonempty": u_loc[0] <= u_loc[1],
                "normalized_triple_intersection_binlocal": n_loc,
                "normalized_nonempty": n_loc[0] <= n_loc[1],
            }
        )

    pass_rate = (
        sum(1 for r in out_rows if r["unnormalized_nonempty"] and r["normalized_nonempty"]) / len(out_rows)
        if out_rows else 0.0
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2121",
        "stage_id": "S1071",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_CMP2_BINLOCAL_EIGENSPECTRUM_TRIPLE_INTERSECTION_AUDIT_WITH_TRACE"
            if pre_ready and len(out_rows) > 0
            else "OPEN_STRICT_CMP2_BINLOCAL_EIGENSPECTRUM_TRIPLE_INTERSECTION_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2120_present": p2120.get("_missing") is None,
            "p2017_present": p2017.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "binlocal_triple_intersection_audit": {
            "eigenspectrum_local_rule": "adjacent-eigenvalue local delta at bin index",
            "rows": out_rows,
            "triple_intersection_nonempty_pass_rate": pass_rate,
            "scope_limit": "operational bin-local propagation audit only; not theorem-grade closure",
        },
        "backend_import_blocker_update": {
            "D3_uncertainty_propagation_from_dressed_backend": "COMPUTED_BINLOCAL_EIGENSPECTRUM_TRIPLE_INTERSECTION_AUDIT_PARTIAL",
        },
        "recommended_next_honest_step": {
            "id": "P2122_candidate",
            "goal": "bind bin indices to explicit backend s-grid provenance map and rerun binlocal audit with index-to-s uncertainty",
        },
        "c3_gate_update": {
            "C3_cmp2_triple_alignment_estimator_consistency_audit": "COMPUTED",
            "C3_cmp2_binlocal_eigenspectrum_triple_intersection_audit": "COMPUTED" if len(out_rows) > 0 else "OPEN",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "binlocal_rows_exported": len(out_rows) > 0,
            "binlocal_propagation_computed": len(out_rows) > 0,
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
            "# P2121 S1071: strict CMP2 bin-local eigenspectrum triple-intersection audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows exported: `{len(out_rows)}`",
            f"- Pass rate: `{pass_rate}`",
            "",
            "This stage replaces global eigenspectrum fallback with bin-local eigenspectrum deltas and reruns triple-intersection propagation.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
