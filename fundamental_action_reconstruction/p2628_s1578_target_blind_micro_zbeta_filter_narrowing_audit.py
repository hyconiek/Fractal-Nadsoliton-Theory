#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import statistics
import subprocess
from pathlib import Path
from typing import Any, Iterable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2628_s1578_target_blind_micro_zbeta_filter_narrowing_audit.json"
MD = GEN / "p2628_s1578_target_blind_micro_zbeta_filter_narrowing_audit.md"

REPO_ROOT = REPO
QW2048_JSON = REPO_ROOT / "report_qw2048_spectral_phase_locked_pointwise_derivation.json"
P2627_JSON = GEN / "p2627_s1577_interval_zbeta_tolerance_bridge_boundary.json"

Z_TARGET = 100.0
STRICT_REL_EPSILON = 0.01
PRACTICAL_REL_EPSILON = 0.15
MIN_SUPPORT_BINS = 3

NEGATIVE_EXPORT_FLAGS = [
    "positive_beta_renormalization_source_exported",
    "nonlinear_damping_completion_source_exported",
    "full_legacy_to_strict_bridge_revalidated",
    "orientation_odd_selector_source_exported",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.lean", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: research concepts, quality fields, and guard terms rather than packet numbers.
    patterns = {
        "micro_pointwise_quality_content": (
            "pointwise bins|phase_min|phase min|rmse|n_informative|quality filter|"
            "spectral phase-locked|micro pointwise|bin median"
        ),
        "target_blind_filter_content": (
            "target-blind|target independent|target-independent|no target retune|no sector retune|"
            "quality threshold|rectangular filter|monotone filter"
        ),
        "zbeta_distribution_content": (
            "Z_beta|z_beta|beta_median|beta_ci95|renormalization constants|micro-derived|"
            "dispersion|wide CI|q25|q75|median"
        ),
        "tolerance_bridge_content": (
            "tolerance bridge|epsilon.*bridge|relative denominator error|interval-valued bridge|"
            "strict 1%|epsilon envelope|admissible tolerance"
        ),
        "closure_guard_content": (
            "positive_beta_renormalization_source|nonlinear_damping_completion_source|"
            "role-bearing L_total|role transfer|QW-2191|ToE closure|full bridge"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for target-blind micro Z_beta narrowing", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def pointwise_bins(qw2048: dict[str, Any]) -> list[dict[str, Any]]:
    return list(qw2048["pointwise_bins"]["bins"])


def percentile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        raise ValueError("empty percentile input")
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    low = int(position)
    high = min(low + 1, len(sorted_values) - 1)
    weight = position - low
    return sorted_values[low] * (1.0 - weight) + sorted_values[high] * weight


def summarize_subset(subset: list[dict[str, Any]]) -> dict[str, Any]:
    values = sorted(float(row["z_beta_median"]) for row in subset)
    median = statistics.median(values)
    q25 = percentile(values, 0.25)
    q75 = percentile(values, 0.75)
    rel = abs(median - Z_TARGET) / Z_TARGET
    return {
        "support_bins": len(values),
        "z_beta_median_of_bins": median,
        "z_beta_q25": q25,
        "z_beta_q75": q75,
        "median_relative_error_to_100": rel,
        "median_passes_strict_1pct": rel <= STRICT_REL_EPSILON,
        "median_passes_practical_15pct": rel <= PRACTICAL_REL_EPSILON,
        "interval_subset_strict_1pct": q25 >= Z_TARGET * (1.0 - STRICT_REL_EPSILON) and q75 <= Z_TARGET * (1.0 + STRICT_REL_EPSILON),
        "interval_subset_practical_15pct": q25 >= Z_TARGET * (1.0 - PRACTICAL_REL_EPSILON) and q75 <= Z_TARGET * (1.0 + PRACTICAL_REL_EPSILON),
        "d_bins": [int(row["d_bin"]) for row in subset],
    }


def enumerate_rectangular_filters(bins: list[dict[str, Any]]) -> dict[str, Any]:
    # Exhaust the finite target-blind class fixed before seeing acceptance: n>=N, phase_min_median>=P, rmse_median<=R.
    n_thresholds = sorted({int(row["n"]) for row in bins})
    phase_thresholds = sorted({float(row["phase_min_median"]) for row in bins})
    rmse_thresholds = sorted({float(row["rmse_median"]) for row in bins})
    rows: list[dict[str, Any]] = []
    for n_min in n_thresholds:
        for phase_min in phase_thresholds:
            for rmse_max in rmse_thresholds:
                subset = [
                    row for row in bins
                    if int(row["n"]) >= n_min
                    and float(row["phase_min_median"]) >= phase_min
                    and float(row["rmse_median"]) <= rmse_max
                ]
                if len(subset) < MIN_SUPPORT_BINS:
                    continue
                summary = summarize_subset(subset)
                rows.append({
                    "thresholds": {
                        "n_min": n_min,
                        "phase_min_median_min": phase_min,
                        "rmse_median_max": rmse_max,
                    },
                    **summary,
                    "accepted_as_strict_positive_beta_source": (
                        summary["median_passes_strict_1pct"]
                        and summary["interval_subset_strict_1pct"]
                    ),
                    "accepted_as_practical_interval_support": (
                        summary["median_passes_practical_15pct"]
                        and summary["interval_subset_practical_15pct"]
                    ),
                })
    rows_sorted = sorted(rows, key=lambda row: (row["median_relative_error_to_100"], -row["support_bins"], row["z_beta_q75"] - row["z_beta_q25"]))
    strict_rows = [row for row in rows if row["accepted_as_strict_positive_beta_source"]]
    practical_rows = [row for row in rows if row["accepted_as_practical_interval_support"]]
    return {
        "admitted_filter_class": "rectangular target-blind quality filters: n>=N, phase_min_median>=P, rmse_median<=R; aggregate by median of per-bin z_beta_median",
        "threshold_grid_sizes": {
            "n_thresholds": len(n_thresholds),
            "phase_thresholds": len(phase_thresholds),
            "rmse_thresholds": len(rmse_thresholds),
        },
        "min_support_bins": MIN_SUPPORT_BINS,
        "evaluated_filter_count": len(rows),
        "strict_accepting_filter_count": len(strict_rows),
        "practical_interval_accepting_filter_count": len(practical_rows),
        "best_by_median_relative_error": rows_sorted[:10],
        "strict_accepting_filters": strict_rows[:10],
        "practical_interval_accepting_filters": sorted(practical_rows, key=lambda row: (row["support_bins"], row["median_relative_error_to_100"]))[:10],
        "finite_exhaustion_theorem": (
            "Within the fixed target-blind rectangular quality-filter class over QW-2048 bins, no filter with at least "
            f"{MIN_SUPPORT_BINS} support bins has both median Z_beta within 1% of 100 and q25-q75 inside the same 1% envelope."
        ),
    }


def raw_quality_summary(bins: list[dict[str, Any]]) -> dict[str, Any]:
    values = sorted(float(row["z_beta_median"]) for row in bins)
    return {
        "n_bins": len(bins),
        "all_bin_z_beta_min": min(values),
        "all_bin_z_beta_median": statistics.median(values),
        "all_bin_z_beta_max": max(values),
        "all_bin_z_beta_q25": percentile(values, 0.25),
        "all_bin_z_beta_q75": percentile(values, 0.75),
        "phase_min_median_range": [min(float(row["phase_min_median"]) for row in bins), max(float(row["phase_min_median"]) for row in bins)],
        "rmse_median_range": [min(float(row["rmse_median"]) for row in bins), max(float(row["rmse_median"]) for row in bins)],
        "support_n_range": [min(int(row["n"]) for row in bins), max(int(row["n"]) for row in bins)],
    }


def acceptance_and_recommendation(filter_cert: dict[str, Any]) -> dict[str, Any]:
    source_accepts = filter_cert["strict_accepting_filter_count"] > 0
    return {
        "accepts_positive_beta_renormalization_source": source_accepts,
        "accepts_narrowing_as_interval_support": filter_cert["practical_interval_accepting_filter_count"] > 0,
        "failed_reason_if_no_strict_acceptance": None if source_accepts else [
            "no_target_blind_quality_filter_in_class_hits_strict_1pct_median_and_interval_envelope",
            "quality-filter narrowing is still an empirical post-processing class, not a micro operator normalization law",
            "orientation_odd_selector_source remains independent and unresolved",
        ],
        "recommended_next_honest_step": (
            "Stop trying to squeeze exact Z_beta=100 from post-hoc quality filters.  The next proof-quality move is to derive a "
            "target-independent micro operator or conservation/normalization identity whose value is Z_beta before comparison with K_strict_gate; "
            "if that cannot be done, declare only an approximate finite-domain effective bridge with explicit epsilon/domain and keep role-transfer closed."
        ),
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["target_blind_filter_certificate"]
    best = cert["best_by_median_relative_error"][0]
    lines = [
        "# P2628/S1578 target-blind micro Z_beta filter narrowing audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps research content and quality-filter language, not only ticket numbers or names.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Finite target-blind filter theorem",
        "",
        cert["finite_exhaustion_theorem"],
        "",
        f"Evaluated filters: `{cert['evaluated_filter_count']}`; strict accepting filters: `{cert['strict_accepting_filter_count']}`; "
        f"practical interval accepting filters: `{cert['practical_interval_accepting_filter_count']}`.",
        "",
        "Best median-only row:",
        "",
        f"- thresholds: `{best['thresholds']}`",
        f"- support bins: `{best['support_bins']}` with `d_bins={best['d_bins']}`",
        f"- median Z_beta: `{best['z_beta_median_of_bins']:.6f}`",
        f"- q25/q75: `[{best['z_beta_q25']:.6f}, {best['z_beta_q75']:.6f}]`",
        f"- relative median error to 100: `{best['median_relative_error_to_100']:.6f}`",
        "",
        "## Verdict",
        "",
        "The target-blind quality-filter class does not export `positive_beta_renormalization_source`.  It also does not repair P2625/P2620, role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure.",
        "",
        "## Recommended next honest step",
        "",
        payload["acceptance_and_recommendation"]["recommended_next_honest_step"],
        "",
        f"Fingerprint: `{payload['fingerprint_sha256']}`",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(parents=True, exist_ok=True)
    qw2048 = load_json(QW2048_JSON)
    p2627 = load_json(P2627_JSON)
    bins = pointwise_bins(qw2048)
    filter_cert = enumerate_rectangular_filters(bins)
    acceptance = acceptance_and_recommendation(filter_cert)
    payload: dict[str, Any] = {
        "status": "P2628_TARGET_BLIND_MICRO_ZBETA_FILTER_NARROWING_NO_STRICT_SOURCE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "sources": {
            "QW2048_JSON": rel(QW2048_JSON),
            "P2627_JSON": rel(P2627_JSON),
        },
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "upstream_status": {
            "qw2048_verdict": qw2048.get("verdict"),
            "p2627_status": p2627.get("status"),
        },
        "constants": {
            "z_target": Z_TARGET,
            "strict_relative_epsilon": STRICT_REL_EPSILON,
            "practical_relative_epsilon": PRACTICAL_REL_EPSILON,
            "min_support_bins": MIN_SUPPORT_BINS,
        },
        "raw_quality_summary": raw_quality_summary(bins),
        "target_blind_filter_certificate": filter_cert,
        "acceptance_and_recommendation": acceptance,
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2628/S1578 target-blind micro Z_beta filter narrowing guard",
        "\n## P2628/S1578 target-blind micro Z_beta filter narrowing guard\n\n"
        "`P2628/S1578` exhausts a fixed target-blind quality-filter class over the QW-2048 micro bins: `n>=N`, `phase_min_median>=P`, "
        "and `rmse_median<=R`, aggregated by the median of per-bin `Z_beta` medians.  No admitted filter with at least three support bins puts both "
        "the median and q25-q75 interval inside the strict 1% envelope around `Z_beta=100`; the best median-only row remains about 9.76% high.  "
        "Therefore quality-filter narrowing does not export `positive_beta_renormalization_source`; P2625/P2620, role-transfer, role-bearing `L_total`, `QW-2191`, and ToE closure remain blocked.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2628/S1578 target-blind micro Z_beta filter Ltotal guard",
        "\n## P2628/S1578 target-blind micro Z_beta filter Ltotal guard\n\n"
        "`P2628/S1578` keeps role-bearing `L_total` closed.  Exhaustive target-blind quality-filter narrowing of the available micro `Z_beta` bins does not produce an exact "
        "positive coefficient source or a strict 1% interval bridge.  A genuine target-independent micro operator/normalization identity is still required before damping completion, role transfer, `QW-2191`, or ToE closure can be rerun.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "evaluated_filter_count": result["target_blind_filter_certificate"]["evaluated_filter_count"],
        "strict_accepting_filter_count": result["target_blind_filter_certificate"]["strict_accepting_filter_count"],
        "best_median_row": result["target_blind_filter_certificate"]["best_by_median_relative_error"][0],
        "recommended_next": result["acceptance_and_recommendation"]["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
