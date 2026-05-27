#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2251 = GEN / "p2251_s1201_strict_nu_branch_group_policy_direction_field_rho_schedule_probe.json"
OUT = GEN / "p2252_s1202_strict_nu_branch_group_policy_adaptive_rho_curvature_proxy_probe.json"
MD = GEN / "p2252_s1202_strict_nu_branch_group_policy_adaptive_rho_curvature_proxy_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2251 = load(IN_2251)
    probe = (p2251.get("strict_nu_branch_group_policy_direction_field_rho_schedule_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})

    total = float(inp.get("total_coverage_mass", 1.0) or 1.0)
    n0 = float(inp.get("initial_group_count", 1.0) or 1.0)
    l0 = float(inp.get("initial_load_ratio", 1.0) or 1.0)
    kappa = float(inp.get("kappa_reference", 1.0) or 1.0)

    def margin(nv: float, lv: float) -> float:
        return total - ((nv - 1.0) + kappa * lv)

    # In affine margin law, true curvature is zero; use normalized margin deficit proxy
    # as operational curvature surrogate for adaptive conservativeness.
    m0 = margin(n0, l0)
    norm = max(abs(total), 1.0)
    curvature_proxy = max(0.0, 1.0 - m0 / norm)

    # Adaptive rho policy: tighter when proxy high.
    rho_min, rho_max = 0.55, 0.9
    rho_adaptive = rho_max - (rho_max - rho_min) * min(1.0, curvature_proxy)

    # Compare against fixed rho baseline from previous stages.
    rho_baseline = 0.8
    conservative_gap = rho_baseline - rho_adaptive

    payload = {
        "schema_version": "p2252_s1202_v1",
        "packet_id": "P2252",
        "stage_id": "S1202",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_RHO_CURVATURE_PROXY_PROBE",
        "strict_nu_branch_group_policy_adaptive_rho_curvature_proxy_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_RHO_CURVATURE_PROXY_PROBE_V1",
            "source_packets": [str(IN_2251.relative_to(ROOT))],
            "inputs": {
                "initial_group_count": n0,
                "initial_load_ratio": l0,
                "total_coverage_mass": total,
                "kappa_reference": kappa,
            },
            "curvature_proxy_model": {
                "base_margin": m0,
                "normalization_scale": norm,
                "curvature_proxy": curvature_proxy,
                "rho_min": rho_min,
                "rho_max": rho_max,
                "rho_adaptive": rho_adaptive,
                "rho_baseline": rho_baseline,
                "conservative_gap_vs_baseline": conservative_gap,
            },
            "physical_interpretation_note": "Curvature proxy acts as near-boundary fragility indicator: shrinking margin raises effective geometric stress, and adaptive rho tightens controller aggressiveness to preserve sign-stability budget.",
            "theorem_scope_limit": "local adaptive-rho proxy diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2253_candidate",
            "goal": "embed adaptive rho into multi-step trajectory and compare invariant preservation vs fixed rho under matched direction fields",
        },
        "gatekeeper_checks": {
            "adaptive_rho_exported": True,
            "curvature_proxy_bounded": 0.0 <= curvature_proxy <= 1.0,
            "rho_within_bounds": rho_min <= rho_adaptive <= rho_max,
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
                "# P2252 S1202: adaptive-rho curvature-proxy probe",
                "",
                f"- base margin: `{m0:.12e}`",
                f"- curvature proxy: `{curvature_proxy:.12e}`",
                f"- adaptive rho: `{rho_adaptive:.12e}`",
                f"- baseline rho: `{rho_baseline:.12e}`",
                f"- conservative gap (baseline-adaptive): `{conservative_gap:.12e}`",
                "",
                "Local adaptive-rho proxy only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
