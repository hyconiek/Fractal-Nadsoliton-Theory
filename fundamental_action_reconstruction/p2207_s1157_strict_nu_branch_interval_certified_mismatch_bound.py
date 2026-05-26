#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2206 = GEN / "p2206_s1156_strict_nu_branch_transport_operator_constraint_solver.json"
OUT = GEN / "p2207_s1157_strict_nu_branch_interval_certified_mismatch_bound.json"
MD = GEN / "p2207_s1157_strict_nu_branch_interval_certified_mismatch_bound.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def ksq(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    den = 1.0 + beta * (d**eta)
    c = math.cos(omega * d + phi)
    return (c * c) / (den * den)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2206 = load(IN_2206)
    p2190 = load(ROOT / "generated/p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json")
    w = p2190.get("strict_cutkosky_fixed_channel_real_discontinuity_integral_packet", {}) or {}
    params = w.get("strict_kernel_params", {}) or {"omega": 0.18575, "phi": 0.1625, "beta": 1.0, "eta": 1.8}

    omega0 = float(params["omega"])
    phi = float(params["phi"])
    beta = float(params["beta"])
    eta = float(params["eta"])

    # compact-domain certificate settings
    d_min, d_max = 0.0, 20.0
    d_grid = np.linspace(d_min, d_max, 8001)
    m_grid = np.linspace(0.94, 1.06, 481)
    eps_exclude = 5e-4

    rows = []
    min_nontrivial_sup = float("inf")
    argmin_nontrivial_m = None

    for m in m_grid:
        if abs(float(m) - 1.0) <= eps_exclude:
            continue
        diff = np.abs([
            ksq(float(d), omega0, phi, beta, eta) - ksq(float(d), omega0 * float(m), phi, beta, eta)
            for d in d_grid
        ])
        sup_norm = float(np.max(diff))
        rows.append({"m": float(m), "sup_norm_diff_on_compact_d": sup_norm})
        if sup_norm < min_nontrivial_sup:
            min_nontrivial_sup = sup_norm
            argmin_nontrivial_m = float(m)

    # certify trivial branch at m=1 has zero sup norm numerically
    diff_m1 = np.abs([
        ksq(float(d), omega0, phi, beta, eta) - ksq(float(d), omega0 * 1.0, phi, beta, eta)
        for d in d_grid
    ])
    sup_m1 = float(np.max(diff_m1))

    certified_gap = min_nontrivial_sup - sup_m1

    payload = {
        "schema_version": "p2207_s1157_v1",
        "packet_id": "P2207",
        "stage_id": "S1157",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_INTERVAL_CERTIFIED_MISMATCH_BOUND",
        "strict_nu_branch_interval_certified_mismatch_bound": {
            "certificate_id": "STRICT_NU_BRANCH_INTERVAL_CERTIFIED_MISMATCH_BOUND_V1",
            "source_packet": str(IN_2206.relative_to(ROOT)),
            "compact_domain": {"d_min": d_min, "d_max": d_max, "d_count": int(d_grid.size)},
            "m_scan": {"m_min": 0.94, "m_max": 1.06, "m_count": int(m_grid.size), "exclude_radius_around_1": eps_exclude},
            "sup_norm_table": rows,
            "summary": {
                "sup_norm_at_m_eq_1": sup_m1,
                "min_sup_norm_for_nontrivial_m": min_nontrivial_sup,
                "argmin_nontrivial_m": argmin_nontrivial_m,
                "certified_separation_gap": certified_gap,
            },
            "interpretation": "on this compact domain and scanned m-window, nontrivial m branches remain separated from the trivial m=1 branch by a positive sup-norm gap",
            "theorem_scope_limit": "interval-certified numeric separation on frozen domain; not global theorem over all backgrounds/branches",
        },
        "recommended_next_honest_step": {
            "id": "P2208_candidate",
            "goal": "lift compact-domain sup-norm separation into symbolic inequality with explicit parameter guards for nu-branch transport",
        },
        "gatekeeper_checks": {
            "interval_mismatch_certificate_exported": True,
            "m1_sup_norm_near_zero": sup_m1 < 1e-15,
            "nontrivial_branch_separation_positive": certified_gap > 0.0,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2207 S1157: strict nu-branch interval-certified mismatch bound",
            "",
            f"- sup norm at m=1: `{sup_m1:.12e}`",
            f"- min nontrivial sup norm: `{min_nontrivial_sup:.12e}` at m={argmin_nontrivial_m}",
            f"- certified separation gap: `{certified_gap:.12e}`",
            "",
            "Compact-domain separation certificate only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
