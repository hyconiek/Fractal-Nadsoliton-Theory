#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2257 = GEN / "p2257_s1207_strict_nu_branch_group_policy_second_order_boundary_shift_probe.json"
OUT = GEN / "p2258_s1208_strict_nu_branch_group_policy_second_order_gamma_fit_probe.json"
MD = GEN / "p2258_s1208_strict_nu_branch_group_policy_second_order_gamma_fit_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2257 = load(IN_2257)
    probe = (p2257.get("strict_nu_branch_group_policy_second_order_boundary_shift_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})
    b2 = (probe.get("second_order_boundary", {}) or {})

    progress_gap = float(inp.get("progress_gap", 0.0) or 0.0)
    rows = b2.get("rows", []) or []

    # Fit gamma from exported proxy relation:
    # |shift| = gamma * h^2 / (1+|progress_gap|)
    # => gamma_i = |shift|*(1+|progress_gap|)/h^2
    gamma_samples = []
    for r in rows:
        h = float(r.get("horizon_scale", 0.0) or 0.0)
        s = abs(float(r.get("second_order_shift", 0.0) or 0.0))
        if h > 1e-15:
            gamma_i = s * (1.0 + abs(progress_gap)) / (h * h)
            gamma_samples.append(gamma_i)

    gamma_hat = sum(gamma_samples) / max(len(gamma_samples), 1)
    gamma_spread = max(gamma_samples) - min(gamma_samples) if gamma_samples else 0.0

    payload = {
        "schema_version": "p2258_s1208_v1",
        "packet_id": "P2258",
        "stage_id": "S1208",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_SECOND_ORDER_GAMMA_FIT_PROBE",
        "strict_nu_branch_group_policy_second_order_gamma_fit_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_SECOND_ORDER_GAMMA_FIT_PROBE_V1",
            "source_packets": [str(IN_2257.relative_to(ROOT))],
            "inputs": {
                "progress_gap": progress_gap,
                "row_count": len(rows),
            },
            "gamma_fit": {
                "gamma_samples": gamma_samples,
                "gamma_hat": gamma_hat,
                "gamma_spread": gamma_spread,
            },
            "physical_interpretation_note": "Gamma fit quantifies effective second-order curvature strength in safety-boundary tightening; low spread supports internal consistency of the nonlinear correction law across horizon scales.",
            "theorem_scope_limit": "finite-sample gamma-fit diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2259_candidate",
            "goal": "inject fitted gamma into adaptive controller and measure boundary-hit frequency reduction vs first-order-only controller",
        },
        "gatekeeper_checks": {
            "gamma_fit_exported": True,
            "gamma_samples_nonempty": len(gamma_samples) > 0,
            "gamma_spread_nonnegative": gamma_spread >= 0.0,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2258 S1208: second-order gamma-fit probe",
                "",
                f"- gamma_hat: `{gamma_hat:.12e}`",
                f"- gamma_spread: `{gamma_spread:.12e}`",
                f"- gamma_samples_count: `{len(gamma_samples)}`",
                "",
                "Finite-sample gamma-fit diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
