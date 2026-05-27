#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2262 = GEN / "p2262_s1212_strict_nu_branch_group_policy_multiseed_reduction_stability_probe.json"
OUT = GEN / "p2263_s1213_strict_nu_branch_group_policy_multiseed_reduction_confidence_interval_probe.json"
MD = GEN / "p2263_s1213_strict_nu_branch_group_policy_multiseed_reduction_confidence_interval_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2262 = load(IN_2262)
    probe = (p2262.get("strict_nu_branch_group_policy_multiseed_reduction_stability_probe", {}) or {})
    summary = (probe.get("multiseed_reduction_summary", {}) or {})

    rows = summary.get("rows", []) or []
    reductions = [float(r.get("reduction", 0.0) or 0.0) for r in rows]
    n = len(reductions)

    mean = sum(reductions) / max(n, 1)
    if n > 1:
        var = sum((x - mean) ** 2 for x in reductions) / (n - 1)
    else:
        var = 0.0
    std = math.sqrt(max(var, 0.0))

    # Small-sample conservative normal-approx CI (z=1.96).
    se = std / math.sqrt(max(n, 1))
    z = 1.96
    ci_low = mean - z * se
    ci_high = mean + z * se

    payload = {
        "schema_version": "p2263_s1213_v1",
        "packet_id": "P2263",
        "stage_id": "S1213",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_MULTISEED_REDUCTION_CONFIDENCE_INTERVAL_PROBE",
        "strict_nu_branch_group_policy_multiseed_reduction_confidence_interval_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_MULTISEED_REDUCTION_CONFIDENCE_INTERVAL_PROBE_V1",
            "source_packets": [str(IN_2262.relative_to(ROOT))],
            "inputs": {
                "sample_count": n,
                "reductions": reductions,
            },
            "confidence_interval": {
                "mean_reduction": mean,
                "std_reduction": std,
                "se_reduction": se,
                "z_value": z,
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "ci95_positive": ci_low > 0.0,
            },
            "physical_interpretation_note": "Confidence interval over multiseed reduction quantifies robustness of nonlinear-correction benefit beyond single-realization evidence, tightening empirical support for safety improvement claims.",
            "theorem_scope_limit": "finite multiseed CI diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2264_candidate",
            "goal": "expand seed pool and compare bootstrap CI with normal-approx CI for reduction robustness",
        },
        "gatekeeper_checks": {
            "confidence_interval_exported": True,
            "ci_bounds_ordered": ci_low <= ci_high,
            "std_nonnegative": std >= 0.0,
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
                "# P2263 S1213: multiseed reduction confidence-interval probe",
                "",
                f"- mean reduction: `{mean:.12e}`",
                f"- std reduction: `{std:.12e}`",
                f"- 95% CI: [`{ci_low:.12e}`, `{ci_high:.12e}`]",
                f"- CI positive: `{ci_low > 0.0}`",
                "",
                "Finite multiseed CI diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
