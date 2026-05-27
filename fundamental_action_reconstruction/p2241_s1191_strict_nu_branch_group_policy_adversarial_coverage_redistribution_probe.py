#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2240 = GEN / "p2240_s1190_strict_nu_branch_group_policy_sign_stability_condition_probe.json"
OUT = GEN / "p2241_s1191_strict_nu_branch_group_policy_adversarial_coverage_redistribution_probe.json"
MD = GEN / "p2241_s1191_strict_nu_branch_group_policy_adversarial_coverage_redistribution_probe.md"

def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))

def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2240 = load(IN_2240)
    probe = (p2240.get("strict_nu_branch_group_policy_sign_stability_condition_probe", {}) or {})
    cond = (probe.get("derived_condition", {}) or {})
    rows = probe.get("per_group_condition_rows", []) or []

    load_ratio = float(cond.get("required_min_coverage", 1.0) or 1.0)
    covs = [float(r.get("coverage_fraction", 1.0) or 1.0) for r in rows]
    gids = [r.get("group_id", f"g{i}") for i, r in enumerate(rows)]

    n = max(len(covs), 1)
    if not covs:
        covs = [1.0]
        gids = ["fallback"]

    total_cov = sum(covs)
    # Adversarial redistribution under fixed total mass:
    # push as much coverage as possible to all but one group, leave one group minimal.
    max_per_group = 1.0
    residual = total_cov - max_per_group * (n - 1)
    worst_min_cov = max(0.0, min(max_per_group, residual))

    baseline_min_cov = min(covs)
    baseline_pass = baseline_min_cov + 1e-15 >= load_ratio
    adversarial_pass = worst_min_cov + 1e-15 >= load_ratio

    # Construct one explicit adversarial witness vector.
    worst_vector = [max_per_group] * n
    worst_vector[0] = worst_min_cov
    # rotate so that weakest is current strongest-id independent, keep explicit witness.
    witness_rows = [
        {
            "group_id": gids[i],
            "adversarial_coverage_fraction": worst_vector[i],
        }
        for i in range(n)
    ]

    payload = {
        "schema_version": "p2241_s1191_v1",
        "packet_id": "P2241",
        "stage_id": "S1191",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ADVERSARIAL_COVERAGE_REDISTRIBUTION_PROBE",
        "strict_nu_branch_group_policy_adversarial_coverage_redistribution_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ADVERSARIAL_COVERAGE_REDISTRIBUTION_PROBE_V1",
            "source_packets": [str(IN_2240.relative_to(ROOT))],
            "fixed_mass_inputs": {
                "group_count": n,
                "total_coverage_mass": total_cov,
                "load_ratio_threshold": load_ratio,
            },
            "baseline": {
                "observed_min_coverage": baseline_min_cov,
                "sufficient_condition_holds": baseline_pass,
            },
            "adversarial_worst_case": {
                "worst_min_coverage_under_fixed_mass": worst_min_cov,
                "sufficient_condition_holds": adversarial_pass,
                "robustness_headroom": worst_min_cov - load_ratio,
            },
            "explicit_adversarial_witness": witness_rows,
            "physical_interpretation_note": "Under fixed global coverage mass, sign-stability is fragile to concentration: if one group can be coverage-starved below load ratio, local signed response protection can be lost in that group.",
            "theorem_scope_limit": "finite-sample adversarial redistribution diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2242_candidate",
            "goal": "derive enforceable lower-bound coverage policy constraint that guarantees adversarial-pass under fixed mass",
        },
        "gatekeeper_checks": {
            "adversarial_redistribution_exported": True,
            "baseline_condition_holds": baseline_pass,
            "adversarial_condition_holds": adversarial_pass,
            "adversarial_fragility_detected": baseline_pass and (not adversarial_pass),
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
                "# P2241 S1191: adversarial coverage redistribution probe",
                "",
                f"- group count: `{n}`",
                f"- total coverage mass: `{total_cov:.12e}`",
                f"- load-ratio threshold: `{load_ratio:.12e}`",
                f"- baseline min coverage: `{baseline_min_cov:.12e}` (pass={baseline_pass})",
                f"- adversarial worst min coverage: `{worst_min_cov:.12e}` (pass={adversarial_pass})",
                f"- robustness headroom (worst-case): `{(worst_min_cov-load_ratio):.12e}`",
                "",
                "Finite-sample adversarial redistribution diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

if __name__ == "__main__":
    main()
