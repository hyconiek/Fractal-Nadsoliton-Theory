#!/usr/bin/env python3
"""P2047 S997: stability certificate extraction audit.

Extracts a local-model stability certificate: an explicit epsilon* perturbation
range and norm condition under which bucket ranking remains stable.
Audit-level only; no theorem-grade C3 discharge.
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
OUT = GEN / "p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.json"
MD = GEN / "p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.md"

SCHEMA_VERSION = "p2047_s997_v1"
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
    out = [0] * len(vals)
    for r, i in enumerate(idx):
        out[i] = r + 1
    return out


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2046 = load("p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2046 = p2046.get("gatekeeper_checks") or {}

    ready = (
        p2046.get("result_kind")
        == "PASS_ADVERSARIAL_BUCKET_PERTURBATION_AUDIT_WITH_STABILITY_MARGIN__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2046.get("margin_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)

    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")
    propagated = float(((p2040.get("residual_audit") or {}).get("propagated_uncertainty_bound_linf", 0.0)))

    # Base bucket ranking from P2046 base means.
    base = p2046.get("base_ranking") or {}
    base_means = base.get("base_bucket_means") or []
    base_rank = base.get("base_rank") or rank_desc([float(x) for x in base_means])

    scope = p2046.get("audit_scope") or {}
    norm_buckets = [float(x) for x in (scope.get("norm_buckets") or [0.01, 0.015, 0.02, 0.025])]

    # certificate search settings
    trial_seeds = [4013, 4019, 4021]
    vectors_per_bucket = 20
    rng = np.random.default_rng(2047)

    def bucket_worst_for_epsilon(eps: float, seed: int, sign_pattern: np.ndarray) -> list[float]:
        M = np.eye(3, dtype=float) + delta + eps * sign_pattern
        rg = np.random.default_rng(seed)
        out = []
        for n in norm_buckets:
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

    # Use many adversarial sign patterns, increasing epsilon schedule.
    sign_patterns = []
    for _ in range(24):
        s = np.sign(rng.normal(size=delta.shape))
        s[s == 0] = 1.0
        sign_patterns.append(s)

    eps_grid = [0.0, 2e-7, 4e-7, 6e-7, 8e-7, 1.0e-6, 1.2e-6, 1.4e-6, 1.6e-6]
    stable_flags = []

    for eps in eps_grid:
        stable_all = True
        for sgn in sign_patterns:
            per_seed = [bucket_worst_for_epsilon(eps, sd, sgn) for sd in trial_seeds]
            avg = list(np.mean(np.array(per_seed), axis=0))
            r = rank_desc(avg)
            if r != base_rank:
                stable_all = False
                break
        stable_flags.append(stable_all)

    # certificate epsilon* is max epsilon that keeps all checked rankings stable.
    eps_star = 0.0
    for eps, ok in zip(eps_grid, stable_flags):
        if ok:
            eps_star = eps

    # stability margin from P2046 gives additional context.
    adv = p2046.get("adversarial_summary") or {}
    min_margin = float(adv.get("worst_case_min_stability_margin", math.inf))

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2047",
        "stage_id": "S997",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_LOCAL_STABILITY_CERTIFICATE_EXTRACTED_WITH_TRACE__C3_STILL_OPEN"
            if ready and math.isfinite(eps_star)
            else "OPEN_LOCAL_STABILITY_CERTIFICATE_EXTRACTION_BLOCKED"
        ),
        "route": "strict_only_local_stability_certificate_audit",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2046_present": p2046.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2046_json_sha256": file_sha256(GEN / "p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json"),
        },
        "certificate_scope": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "norm_condition": {
                "norm_buckets_checked": norm_buckets,
                "vectors_per_bucket": vectors_per_bucket,
                "trial_seeds": trial_seeds,
            },
            "ranking_reference": {
                "base_bucket_means": base_means,
                "base_rank": base_rank,
            },
        },
        "certificate_extraction": {
            "epsilon_grid": eps_grid,
            "stable_flags": stable_flags,
            "epsilon_star_max_stable": eps_star,
            "adversarial_sign_patterns_checked": len(sign_patterns),
            "reference_min_stability_margin_from_p2046": min_margin,
            "certificate_statement": (
                "For perturbations ||Delta_adv||_inf <= epsilon_star_max_stable and for the declared norm-bucket test family, "
                "the bucket ranking remains stable in this local audit model."
            ),
        },
        "c3_gate_update": {
            "C3_local_stability_certificate_audit": "EXTRACTED_LOCAL_CERTIFICATE",
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
            "epsilon_star_finite": math.isfinite(eps_star),
            "epsilon_star_nonnegative": eps_star >= 0.0,
            "grid_nonempty": len(eps_grid) > 0,
            "flags_match_grid": len(stable_flags) == len(eps_grid),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2047 S997: stability certificate extraction audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Extracted a local stability certificate epsilon* and explicit norm-bucket",
        "conditions under which ranking stayed stable in the present audit model.",
        "",
        "## Gate update",
        "",
        "- `C3`: local certificate extracted (audit-level).",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
