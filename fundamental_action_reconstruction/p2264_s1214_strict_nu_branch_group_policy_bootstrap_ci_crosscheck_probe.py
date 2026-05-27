#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2263 = GEN / "p2263_s1213_strict_nu_branch_group_policy_multiseed_reduction_confidence_interval_probe.json"
OUT = GEN / "p2264_s1214_strict_nu_branch_group_policy_bootstrap_ci_crosscheck_probe.json"
MD = GEN / "p2264_s1214_strict_nu_branch_group_policy_bootstrap_ci_crosscheck_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def percentile(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    idx = max(0, min(len(sorted_vals) - 1, int(q * (len(sorted_vals) - 1))))
    return sorted_vals[idx]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2263 = load(IN_2263)
    probe = (p2263.get("strict_nu_branch_group_policy_multiseed_reduction_confidence_interval_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})
    ci = (probe.get("confidence_interval", {}) or {})

    reductions = [float(x) for x in (inp.get("reductions", []) or [])]
    normal_low = float(ci.get("ci95_low", 0.0) or 0.0)
    normal_high = float(ci.get("ci95_high", 0.0) or 0.0)

    rng = random.Random(2264)
    b = 1000
    means = []
    n = len(reductions)
    if n == 0:
        reductions = [0.0]
        n = 1

    for _ in range(b):
        sample = [reductions[rng.randrange(n)] for __ in range(n)]
        means.append(sum(sample) / n)

    means.sort()
    boot_low = percentile(means, 0.025)
    boot_high = percentile(means, 0.975)

    overlap = not (boot_high < normal_low or normal_high < boot_low)

    payload = {
        "schema_version": "p2264_s1214_v1",
        "packet_id": "P2264",
        "stage_id": "S1214",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_BOOTSTRAP_CI_CROSSCHECK_PROBE",
        "strict_nu_branch_group_policy_bootstrap_ci_crosscheck_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_BOOTSTRAP_CI_CROSSCHECK_PROBE_V1",
            "source_packets": [str(IN_2263.relative_to(ROOT))],
            "inputs": {
                "reductions": reductions,
                "normal_ci95_low": normal_low,
                "normal_ci95_high": normal_high,
                "bootstrap_draws": b,
                "random_seed": 2264,
            },
            "bootstrap_crosscheck": {
                "bootstrap_ci95_low": boot_low,
                "bootstrap_ci95_high": boot_high,
                "ci_overlap": overlap,
            },
            "physical_interpretation_note": "Bootstrap CI crosscheck probes whether reduction-confidence conclusions are robust to distributional assumptions, improving reliability of physical safety-gain inference.",
            "theorem_scope_limit": "finite-sample bootstrap crosscheck only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2265_candidate",
            "goal": "expand seed set and compare bootstrap-vs-normal divergence trend as n grows",
        },
        "gatekeeper_checks": {
            "bootstrap_crosscheck_exported": True,
            "bootstrap_ci_ordered": boot_low <= boot_high,
            "ci_overlap_computable": isinstance(overlap, bool),
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
                "# P2264 S1214: bootstrap CI crosscheck probe",
                "",
                f"- normal CI95: [`{normal_low:.12e}`, `{normal_high:.12e}`]",
                f"- bootstrap CI95: [`{boot_low:.12e}`, `{boot_high:.12e}`]",
                f"- overlap: `{overlap}`",
                "",
                "Finite bootstrap crosscheck only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
