#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2242 = GEN / "p2242_s1192_strict_nu_branch_group_policy_enforceable_floor_constraint_probe.json"
OUT = GEN / "p2243_s1193_strict_nu_branch_group_policy_constraint_monte_carlo_envelope_probe.json"
MD = GEN / "p2243_s1193_strict_nu_branch_group_policy_constraint_monte_carlo_envelope_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sample_fixed_mass_vector(n: int, total_mass: float, rng: random.Random) -> list[float]:
    caps = [1.0] * n
    vec = [0.0] * n
    remaining = total_mass
    for i in range(n - 1):
        suffix_cap = sum(caps[i + 1 :])
        lo = max(0.0, remaining - suffix_cap)
        hi = min(caps[i], remaining)
        x = lo + (hi - lo) * rng.random()
        vec[i] = x
        remaining -= x
    vec[-1] = max(0.0, min(caps[-1], remaining))
    return vec


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2242 = load(IN_2242)
    probe = (p2242.get("strict_nu_branch_group_policy_enforceable_floor_constraint_probe", {}) or {})
    inputs = (probe.get("inputs", {}) or {})
    n = int(inputs.get("group_count", 1) or 1)
    total_cov = float(inputs.get("total_coverage_mass", 1.0) or 1.0)
    load_ratio = float(inputs.get("load_ratio_threshold", 1.0) or 1.0)

    draws = 2000
    rng = random.Random(2243)
    violations = 0
    min_cov_samples: list[float] = []

    for _ in range(draws):
        v = sample_fixed_mass_vector(n, total_cov, rng)
        mc_min = min(v)
        min_cov_samples.append(mc_min)
        if mc_min + 1e-15 < load_ratio:
            violations += 1

    violation_prob = violations / max(draws, 1)
    empirical_min = min(min_cov_samples) if min_cov_samples else 1.0
    empirical_max = max(min_cov_samples) if min_cov_samples else 1.0
    empirical_mean = sum(min_cov_samples) / max(len(min_cov_samples), 1)

    payload = {
        "schema_version": "p2243_s1193_v1",
        "packet_id": "P2243",
        "stage_id": "S1193",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_CONSTRAINT_MONTE_CARLO_ENVELOPE_PROBE",
        "strict_nu_branch_group_policy_constraint_monte_carlo_envelope_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_CONSTRAINT_MONTE_CARLO_ENVELOPE_PROBE_V1",
            "source_packets": [str(IN_2242.relative_to(ROOT))],
            "inputs": {
                "group_count": n,
                "total_coverage_mass": total_cov,
                "load_ratio_threshold": load_ratio,
                "draw_count": draws,
                "random_seed": 2243,
            },
            "empirical_violation_envelope": {
                "violation_probability": violation_prob,
                "empirical_min_of_min_coverages": empirical_min,
                "empirical_max_of_min_coverages": empirical_max,
                "empirical_mean_of_min_coverages": empirical_mean,
            },
            "physical_interpretation_note": "Monte-Carlo policy-mixing perturbations estimate how often grouped mixing can starve at least one lane below load-ratio safety threshold under fixed global coverage mass.",
            "theorem_scope_limit": "finite-sample empirical stress envelope only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2244_candidate",
            "goal": "derive analytic tail bound linking violation probability to distance from enforceable mass threshold",
        },
        "gatekeeper_checks": {
            "monte_carlo_envelope_exported": True,
            "draw_count_sufficient": draws >= 1000,
            "violation_probability_computable": 0.0 <= violation_prob <= 1.0,
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
                "# P2243 S1193: Monte-Carlo violation envelope probe",
                "",
                f"- draws: `{draws}` (seed=2243)",
                f"- group count: `{n}`",
                f"- total coverage mass: `{total_cov:.12e}`",
                f"- load-ratio threshold: `{load_ratio:.12e}`",
                f"- violation probability: `{violation_prob:.12e}`",
                f"- empirical min(min_coverage): `{empirical_min:.12e}`",
                f"- empirical max(min_coverage): `{empirical_max:.12e}`",
                f"- empirical mean(min_coverage): `{empirical_mean:.12e}`",
                "",
                "Finite-sample Monte-Carlo envelope only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
