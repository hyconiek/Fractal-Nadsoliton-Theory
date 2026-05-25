#!/usr/bin/env python3
"""P2046 S996: adversarial bucket-perturbation audit.

Perturbs delta_finite_part_candidate_matrix within uncertainty bounds and checks
whether norm-bucket ranking remains stable. Exports a stability margin metric.
Audit-level only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json"
MD = GEN / "p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.md"

SCHEMA_VERSION = "p2046_s996_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(v: Any) -> bool:
    return bool(v is True)


def linf(v: np.ndarray) -> float:
    return float(np.max(np.abs(v)))


def rank_desc(vals: list[float]) -> list[int]:
    idx = sorted(range(len(vals)), key=lambda i: (-vals[i], i))
    pos = [0] * len(vals)
    for r, i in enumerate(idx):
        pos[i] = r + 1
    return pos


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2045 = load("p2045_s995_strict_same_scheme_seed_norm_joint_robustness_audit.json")
    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2045 = p2045.get("gatekeeper_checks") or {}

    ready = (
        p2045.get("result_kind")
        == "PASS_SEED_NORM_JOINT_ROBUSTNESS_AUDIT_WITH_RANK_STABILITY__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2045.get("mean_spearman_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)

    ub = float(((load("p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.json").get("uncertainty_bound") or {}).get("absolute_linf_bound", 0.0)))
    propagated = float(((p2040.get("residual_audit") or {}).get("propagated_uncertainty_bound_linf", 0.0)))

    # Adversarial but bounded perturbation magnitude (small, uncertainty-scale).
    eps = max(ub * 8.0, 1.0e-7)
    base_delta_norm = float(np.max(np.abs(delta)))

    # Build base ranking from P2045 mean per bucket over seeds.
    rows = p2045.get("seed_norm_worst_case_table") or []
    bucket_vals: dict[float, list[float]] = {}
    for r in rows:
        bucket_vals.setdefault(float(r["norm_bucket"]), []).append(float(r["worst_case_residual_with_buffer_linf"]))
    buckets = sorted(bucket_vals.keys())
    base_means = [float(np.mean(bucket_vals[b])) for b in buckets]
    base_rank = rank_desc(base_means)

    rng = np.random.default_rng(2046)
    seeds = [3011, 3019, 3023]
    vectors_per_bucket = 24
    norms = buckets

    M_base = np.eye(3, dtype=float) + delta

    def worst_by_bucket(M: np.ndarray, seed: int) -> list[float]:
        rg = np.random.default_rng(seed)
        out = []
        for n in norms:
            worst = 0.0
            for _ in range(vectors_per_bucket):
                v = rg.normal(size=3)
                nv = float(np.linalg.norm(v))
                if nv == 0.0:
                    continue
                a = n * (v / nv)
                res = linf(a - M.dot(a)) + propagated
                worst = max(worst, res)
            out.append(worst)
        return out

    # adversarial search over sign-pattern perturbations clipped to eps
    adversarial_trials = []
    best_margin = math.inf
    best = None

    for t in range(120):
        sign = np.sign(rng.normal(size=delta.shape))
        sign[sign == 0] = 1.0
        pert = sign * eps
        M_adv = np.eye(3, dtype=float) + delta + pert

        per_seed_worst = [worst_by_bucket(M_adv, s) for s in seeds]
        avg_worst = list(np.mean(np.array(per_seed_worst), axis=0))
        rank = rank_desc(avg_worst)

        # stability margin: minimal adjacent gap in sorted avg_worst (descending)
        sorted_vals = sorted(avg_worst, reverse=True)
        gaps = [sorted_vals[i] - sorted_vals[i + 1] for i in range(len(sorted_vals) - 1)]
        margin = min(gaps) if gaps else math.inf

        stable = rank == base_rank
        adversarial_trials.append(
            {
                "trial": t,
                "stable_ranking": stable,
                "avg_bucket_worst": avg_worst,
                "rank": rank,
                "stability_margin": margin,
            }
        )

        if margin < best_margin:
            best_margin = margin
            best = adversarial_trials[-1]

    stable_count = sum(1 for r in adversarial_trials if r["stable_ranking"])

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2046",
        "stage_id": "S996",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_ADVERSARIAL_BUCKET_PERTURBATION_AUDIT_WITH_STABILITY_MARGIN__C3_STILL_OPEN"
            if ready and best is not None and math.isfinite(best_margin)
            else "OPEN_ADVERSARIAL_BUCKET_PERTURBATION_AUDIT_BLOCKED"
        ),
        "route": "strict_only_adversarial_bucket_perturbation_audit",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2045_present": p2045.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2045_json_sha256": file_sha256(GEN / "p2045_s995_strict_same_scheme_seed_norm_joint_robustness_audit.json"),
        },
        "audit_scope": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "norm_buckets": norms,
            "adversarial_trials": len(adversarial_trials),
            "trial_seeds": seeds,
            "vectors_per_bucket": vectors_per_bucket,
            "perturbation_linf_epsilon": eps,
            "base_delta_linf_norm": base_delta_norm,
        },
        "base_ranking": {
            "base_bucket_means": base_means,
            "base_rank": base_rank,
        },
        "adversarial_summary": {
            "stable_trial_count": stable_count,
            "total_trials": len(adversarial_trials),
            "stable_fraction": stable_count / len(adversarial_trials) if adversarial_trials else 0.0,
            "worst_case_min_stability_margin": best_margin,
            "worst_case_trial": best,
        },
        "trials": adversarial_trials,
        "c3_gate_update": {
            "C3_adversarial_bucket_perturbation_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "trials_nonempty": len(adversarial_trials) > 0,
            "margin_finite": math.isfinite(best_margin),
            "stable_fraction_in_0_1": 0.0 <= (stable_count / len(adversarial_trials)) <= 1.0 if adversarial_trials else True,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2046 S996: adversarial bucket-perturbation audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported adversarially perturbed bucket-ranking stress test within",
        "uncertainty-scale epsilon and a stability-margin summary.",
        "",
        "## Gate update",
        "",
        "- `C3`: adversarial perturbation audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
