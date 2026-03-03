#!/usr/bin/env python3
"""
QW-1947: Coupling-completeness audit + mass-operator frontier under frozen kernel.

Strict rules:
- kernel parameters (omega, phi, beta, eta) frozen from QW-1932,
- no fitting to masses,
- deterministic gamma extraction and deterministic operator maps only,
- complexity-aware comparison of operator classes.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1947_coupling_completeness_and_mass_operator_frontier.json"
OUT_MD = ROOT / "RAPORT_QW1947_COUPLING_COMPLETENESS_AND_MASS_OPERATOR_FRONTIER.md"


PARTICLES = [
    ("Top", 173_000.0),
    ("Bottom", 4_180.0),
    ("Tau", 1_776.9),
    ("Charm", 1_270.0),
    ("Muon", 105.7),
    ("Electron", 0.511),
]

CANONICAL_COUPLINGS = [
    "amplitude_decay",
    "sign_oscillation",
    "phase_offset",
    "parity_structure",
    "local_gradient",
    "local_curvature",
    "memory_cumulative",
    "nonlocal_scale_mix",
]

THRESHOLDS = {
    "mass_mean_rel_pct_max": 15.0,
    "mass_max_rel_pct_max": 75.0,
}


@dataclass(frozen=True)
class OperatorClass:
    op_id: str
    s_map: str
    gamma_key: str
    complexity_k: int
    features: Tuple[str, ...]


OPERATOR_CLASSES: List[OperatorClass] = [
    OperatorClass(
        op_id="O1_hard_linear_local14",
        s_map="linear",
        gamma_key="local_1_to_4",
        complexity_k=1,
        features=("amplitude_decay",),
    ),
    OperatorClass(
        op_id="O2_hard_linear_wls16",
        s_map="linear",
        gamma_key="wls_1_to_6",
        complexity_k=2,
        features=("amplitude_decay", "local_gradient"),
    ),
    OperatorClass(
        op_id="O3_memory_cumulative_wls16",
        s_map="memory",
        gamma_key="wls_1_to_6",
        complexity_k=2,
        features=("amplitude_decay", "memory_cumulative"),
    ),
    OperatorClass(
        op_id="O4_weighted_nonlocal_shortslope",
        s_map="weighted_memory",
        gamma_key="short_slope_mean",
        complexity_k=3,
        features=("amplitude_decay", "local_gradient", "nonlocal_scale_mix"),
    ),
    OperatorClass(
        op_id="O5_oscillatory_phase_wls16",
        s_map="oscillatory",
        gamma_key="wls_1_to_6",
        complexity_k=3,
        features=("amplitude_decay", "sign_oscillation", "phase_offset"),
    ),
    OperatorClass(
        op_id="O6_parity_augmented_wls16",
        s_map="parity",
        gamma_key="wls_1_to_6",
        complexity_k=3,
        features=("amplitude_decay", "parity_structure", "sign_oscillation"),
    ),
    OperatorClass(
        op_id="O7_curvature_augmented_wls16",
        s_map="curvature",
        gamma_key="wls_1_to_6",
        complexity_k=3,
        features=("amplitude_decay", "local_curvature"),
    ),
    OperatorClass(
        op_id="O8_hybrid_full_curvmod",
        s_map="hybrid",
        gamma_key="curvature_modulated",
        complexity_k=6,
        features=(
            "amplitude_decay",
            "sign_oscillation",
            "phase_offset",
            "parity_structure",
            "local_curvature",
            "memory_cumulative",
            "nonlocal_scale_mix",
        ),
    ),
]


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def gamma_from_two_points(k_a: float, k_b: float, d_a: float, d_b: float) -> float:
    slope = math.log(max(k_b, 1e-15) / max(k_a, 1e-15)) / max(d_b - d_a, 1e-15)
    return float(-4.0 * slope / math.log(4.0))


def gamma_wls_short_range(d: np.ndarray, k_abs: np.ndarray, d_max: float = 6.0) -> Dict[str, float]:
    m = d <= d_max
    x = d[m]
    y = np.log(np.clip(k_abs[m], 1e-15, None))
    w = np.clip(k_abs[m], 1e-12, None)

    x_bar = float(np.sum(w * x) / np.sum(w))
    y_bar = float(np.sum(w * y) / np.sum(w))
    num = float(np.sum(w * (x - x_bar) * (y - y_bar)))
    den = float(np.sum(w * (x - x_bar) ** 2))
    slope = num / max(den, 1e-15)
    intercept = y_bar - slope * x_bar

    y_hat = intercept + slope * x
    ss_res = float(np.sum(w * (y - y_hat) ** 2))
    ss_tot = float(np.sum(w * (y - y_bar) ** 2))
    r2 = float(1.0 - ss_res / max(ss_tot, 1e-15))
    gamma = float(-4.0 * slope / math.log(4.0))
    return {"gamma": gamma, "slope": slope, "intercept": intercept, "r2_weighted": r2}


def build_kernel_characteristics(kernel: Dict[str, float]) -> Dict[str, object]:
    d = np.arange(1.0, 25.0, dtype=float)
    k_signed = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    k_abs = np.abs(k_signed)

    log_abs = np.log(np.clip(k_abs, 1e-15, None))
    slope = -(log_abs[1:] - log_abs[:-1])
    curvature = np.abs(log_abs[2:] - 2.0 * log_abs[1:-1] + log_abs[:-2])

    sign_bits = np.sign(k_signed)
    sign_bits[sign_bits == 0.0] = 1.0
    flips = (sign_bits[1:] != sign_bits[:-1]).astype(float)
    flip_count_cum = np.concatenate([[0.0], np.cumsum(flips)])

    c_abs = np.cumsum(k_abs) / max(np.sum(k_abs), 1e-15)
    c_weighted = np.cumsum(d * k_abs) / max(np.sum(d * k_abs), 1e-15)

    c_curv = np.zeros_like(d)
    c_curv[2:] = np.cumsum(curvature)
    c_curv = c_curv / max(c_curv[-1], 1e-15)

    even_sum = float(np.sum(k_abs[(d.astype(int) % 2) == 0]))
    odd_sum = float(np.sum(k_abs[(d.astype(int) % 2) == 1]))
    parity_imbalance = float(abs(even_sum - odd_sum) / max(even_sum + odd_sum, 1e-15))

    short_slope_mean = float(np.mean(slope[:6]))
    short_curv_mean = float(np.mean(curvature[:6]))
    chi_osc = float(np.clip(flip_count_cum[-1] / max(len(d) - 1, 1), 0.0, 1.0))
    chi_curv = float(np.clip(short_curv_mean / max(abs(short_slope_mean), 1e-15), 0.0, 1.0))
    chi_mem = float(np.clip(c_abs[5], 0.0, 1.0))
    chi_parity = float(np.clip(parity_imbalance, 0.0, 1.0))

    g12 = gamma_from_two_points(float(k_abs[0]), float(k_abs[1]), 1.0, 2.0)
    g14 = gamma_from_two_points(float(k_abs[0]), float(k_abs[3]), 1.0, 4.0)
    gwls = gamma_wls_short_range(d, k_abs, d_max=6.0)
    g_short = float(-4.0 * short_slope_mean / math.log(4.0))
    g_curv = float(gwls["gamma"] * (1.0 + 0.20 * chi_curv))

    gammas = {
        "local_1_to_2": float(g12),
        "local_1_to_4": float(g14),
        "wls_1_to_6": float(gwls["gamma"]),
        "short_slope_mean": float(g_short),
        "curvature_modulated": float(g_curv),
    }

    q_axis = np.arange(0, 25, dtype=float)
    flip_rate = np.zeros_like(q_axis)
    curv_rate = np.zeros_like(q_axis)
    for q in range(1, 25):
        flip_rate[q] = flip_count_cum[q - 1] / max(q - 1, 1)
    for q in range(2, 25):
        curv_rate[q] = c_curv[q - 1]

    s_maps = {
        "linear": np.maximum.accumulate(q_axis.copy()),
        "memory": np.maximum.accumulate(np.concatenate([[0.0], 24.0 * c_abs])),
        "weighted_memory": np.maximum.accumulate(np.concatenate([[0.0], 24.0 * c_weighted])),
        "oscillatory": np.maximum.accumulate(q_axis * (1.0 + chi_osc * flip_rate)),
        "parity": np.maximum.accumulate(q_axis * (1.0 + chi_parity * np.abs(np.sin(np.pi * q_axis / 2.0)))),
        "curvature": np.maximum.accumulate(q_axis * (1.0 + chi_curv * curv_rate)),
    }
    s_maps["hybrid"] = np.maximum.accumulate(
        0.35 * s_maps["memory"]
        + 0.25 * s_maps["weighted_memory"]
        + 0.20 * s_maps["oscillatory"]
        + 0.10 * s_maps["parity"]
        + 0.10 * s_maps["curvature"]
    )
    s_maps["hybrid"][0] = 0.0

    return {
        "d_1_to_24": [float(x) for x in d],
        "k_abs_1_to_24": [float(x) for x in k_abs],
        "k_signed_1_to_24": [float(x) for x in k_signed],
        "diagnostics": {
            "short_slope_mean": short_slope_mean,
            "short_curvature_mean": short_curv_mean,
            "sign_flip_rate": chi_osc,
            "short_memory_fraction_d1_to_d6": chi_mem,
            "parity_imbalance": chi_parity,
            "chi_curvature": chi_curv,
        },
        "gamma_recipes": gammas,
        "gamma_wls_fit_quality": {
            "slope": float(gwls["slope"]),
            "intercept": float(gwls["intercept"]),
            "r2_weighted": float(gwls["r2_weighted"]),
        },
        "s_maps": {k: [float(v) for v in arr] for k, arr in s_maps.items()},
    }


def make_mass_scenarios(d1943: Dict[str, object]) -> Dict[str, Dict[str, int]]:
    baseline = d1943["baseline_row"]["mass_q"]
    best_mass = d1943["best_mass_row"]["mass_q"]
    best_joint = d1943["best_joint_row"]["mass_q"]
    return {
        "baseline_assignment": {k: int(v) for k, v in baseline.items()},
        "audit_best_mass_assignment": {k: int(v) for k, v in best_mass.items()},
        "audit_best_joint_assignment": {k: int(v) for k, v in best_joint.items()},
    }


def evaluate_mass_operator(
    op: OperatorClass,
    gamma: float,
    s_map: np.ndarray,
    mass_q: Dict[str, int],
) -> Dict[str, object]:
    rows = []
    errs = []
    log_res = []
    for p_name, exp_mev in PARTICLES:
        q = int(mass_q[p_name])
        s_q = float(s_map[min(max(q, 0), len(s_map) - 1)])
        pred = 173_000.0 * (4.0 ** (-(gamma * s_q / 4.0)))
        err = float(abs(pred - exp_mev) / max(exp_mev, 1e-15) * 100.0)
        errs.append(err)
        log_res.append(float(math.log(max(pred, 1e-15)) - math.log(max(exp_mev, 1e-15))))
        rows.append(
            {
                "particle": p_name,
                "q": int(q),
                "s_q": float(s_q),
                "exp_mev": float(exp_mev),
                "pred_mev": float(pred),
                "rel_err_pct": err,
            }
        )

    n = len(log_res)
    rss = float(np.sum(np.square(np.array(log_res, dtype=float))))
    bic = float(n * math.log(max(rss / max(n, 1), 1e-15)) + op.complexity_k * math.log(max(n, 2)))
    aic = float(n * math.log(max(rss / max(n, 1), 1e-15)) + 2.0 * op.complexity_k)
    mean_rel = float(np.mean(errs))
    max_rel = float(np.max(errs))

    norm_loss = float(mean_rel / THRESHOLDS["mass_mean_rel_pct_max"] + max_rel / THRESHOLDS["mass_max_rel_pct_max"])
    strict_pass = bool(
        mean_rel <= THRESHOLDS["mass_mean_rel_pct_max"] and max_rel <= THRESHOLDS["mass_max_rel_pct_max"]
    )

    return {
        "rows": rows,
        "mean_rel_err_pct": mean_rel,
        "max_rel_err_pct": max_rel,
        "normalized_mass_loss": norm_loss,
        "bic_like": bic,
        "aic_like": aic,
        "strict_mass_pass": strict_pass,
    }


def pareto_front(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    sorted_rows = sorted(rows, key=lambda r: (int(r["strict"]["complexity_k"]), float(r["strict"]["mean_rel_err_pct"])))
    out = []
    best_err = float("inf")
    for r in sorted_rows:
        err = float(r["strict"]["mean_rel_err_pct"])
        if err < best_err - 1e-12:
            out.append(r)
            best_err = err
    return out


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    d1943 = json.loads((ROOT / "report_qw1943_topological_q_assignment_audit.json").read_text(encoding="utf-8"))

    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }

    kchar = build_kernel_characteristics(kernel)
    gamma_recipes = {k: float(v) for k, v in kchar["gamma_recipes"].items()}
    s_maps = {k: np.array(v, dtype=float) for k, v in kchar["s_maps"].items()}
    scenarios = make_mass_scenarios(d1943)

    rows = []
    tested_feature_union = set()
    for op in OPERATOR_CLASSES:
        tested_feature_union.update(op.features)

        per_scenario = []
        for scen_name, mass_q in scenarios.items():
            ev = evaluate_mass_operator(
                op=op,
                gamma=gamma_recipes[op.gamma_key],
                s_map=s_maps[op.s_map],
                mass_q=mass_q,
            )
            per_scenario.append({"scenario": scen_name, "metrics": ev})

        strict_ev = next(x["metrics"] for x in per_scenario if x["scenario"] == "baseline_assignment")
        best_expl = min(per_scenario, key=lambda x: x["metrics"]["normalized_mass_loss"])
        coverage_ratio = float(len(set(op.features)) / len(CANONICAL_COUPLINGS))

        rows.append(
            {
                "operator": {
                    "op_id": op.op_id,
                    "s_map": op.s_map,
                    "gamma_key": op.gamma_key,
                    "gamma_value": gamma_recipes[op.gamma_key],
                    "complexity_k": int(op.complexity_k),
                    "features": list(op.features),
                    "feature_coverage_ratio": coverage_ratio,
                },
                "strict": {
                    "scenario": "baseline_assignment",
                    "mean_rel_err_pct": strict_ev["mean_rel_err_pct"],
                    "max_rel_err_pct": strict_ev["max_rel_err_pct"],
                    "normalized_mass_loss": strict_ev["normalized_mass_loss"],
                    "bic_like": strict_ev["bic_like"],
                    "aic_like": strict_ev["aic_like"],
                    "strict_mass_pass": strict_ev["strict_mass_pass"],
                    "complexity_k": int(op.complexity_k),
                },
                "exploratory_best_scenario": {
                    "scenario": best_expl["scenario"],
                    "mean_rel_err_pct": best_expl["metrics"]["mean_rel_err_pct"],
                    "max_rel_err_pct": best_expl["metrics"]["max_rel_err_pct"],
                    "normalized_mass_loss": best_expl["metrics"]["normalized_mass_loss"],
                    "strict_mass_pass": best_expl["metrics"]["strict_mass_pass"],
                },
                "per_scenario": per_scenario,
            }
        )

    best_strict = min(rows, key=lambda r: (r["strict"]["normalized_mass_loss"], r["strict"]["complexity_k"]))
    best_expl = min(rows, key=lambda r: (r["exploratory_best_scenario"]["normalized_mass_loss"], r["operator"]["complexity_k"]))
    strict_pass_exists = any(bool(r["strict"]["strict_mass_pass"]) for r in rows)
    expl_pass_exists = any(bool(r["exploratory_best_scenario"]["strict_mass_pass"]) for r in rows)

    frontier = pareto_front(rows)
    frontier_short = [
        {
            "op_id": r["operator"]["op_id"],
            "strict_mean_rel_err_pct": r["strict"]["mean_rel_err_pct"],
            "strict_max_rel_err_pct": r["strict"]["max_rel_err_pct"],
            "complexity_k": r["strict"]["complexity_k"],
        }
        for r in frontier
    ]

    tested_feature_union = sorted(tested_feature_union)
    never_tested = sorted(set(CANONICAL_COUPLINGS) - set(tested_feature_union))
    best_missing = sorted(set(CANONICAL_COUPLINGS) - set(best_strict["operator"]["features"]))
    legacy_q1939_features = ["amplitude_decay"]
    legacy_missing = sorted(set(CANONICAL_COUPLINGS) - set(legacy_q1939_features))

    if strict_pass_exists:
        verdict = "COUPLING_AUDIT_AND_MASS_OPERATOR_FRONTIER_STRICT_PASS"
        required_next = "LOCK_BEST_STRICT_OPERATOR_AND_RUN_FLAVOR_LINK_TEST"
    else:
        verdict = "COUPLING_AUDIT_AND_MASS_OPERATOR_FRONTIER_STRICT_FAIL"
        required_next = "MASS_LINK_REWORK_REQUIRED_BEYOND_TESTED_DETERMINISTIC_OPERATOR_CLASSES"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "thresholds": THRESHOLDS,
        "canonical_couplings": CANONICAL_COUPLINGS,
        "kernel_characteristics": {
            "diagnostics": kchar["diagnostics"],
            "gamma_recipes": gamma_recipes,
            "gamma_wls_fit_quality": kchar["gamma_wls_fit_quality"],
        },
        "mass_scenarios": scenarios,
        "legacy_mass_operator_qw1939": {
            "features": legacy_q1939_features,
            "missing_canonical_features": legacy_missing,
        },
        "operator_scan": rows,
        "pareto_frontier_strict": frontier_short,
        "coupling_completeness_audit": {
            "tested_feature_union": tested_feature_union,
            "never_tested_canonical_features": never_tested,
            "best_strict_operator_features": best_strict["operator"]["features"],
            "best_strict_missing_canonical_features": best_missing,
        },
        "summary": {
            "strict_pass_exists": strict_pass_exists,
            "exploratory_pass_exists": expl_pass_exists,
            "best_strict_operator": {
                "op_id": best_strict["operator"]["op_id"],
                "mean_rel_err_pct": best_strict["strict"]["mean_rel_err_pct"],
                "max_rel_err_pct": best_strict["strict"]["max_rel_err_pct"],
                "complexity_k": best_strict["strict"]["complexity_k"],
                "feature_coverage_ratio": best_strict["operator"]["feature_coverage_ratio"],
            },
            "best_exploratory_operator": {
                "op_id": best_expl["operator"]["op_id"],
                "scenario": best_expl["exploratory_best_scenario"]["scenario"],
                "mean_rel_err_pct": best_expl["exploratory_best_scenario"]["mean_rel_err_pct"],
                "max_rel_err_pct": best_expl["exploratory_best_scenario"]["max_rel_err_pct"],
            },
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1947: COUPLING COMPLETENESS + MASS OPERATOR FRONTIER",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Strict (Baseline Assignment)",
        f"- strict_pass_exists: {strict_pass_exists}",
        (
            "- best strict operator: "
            f"{best_strict['operator']['op_id']} | "
            f"mean/max rel%={best_strict['strict']['mean_rel_err_pct']:.3f}/{best_strict['strict']['max_rel_err_pct']:.3f} | "
            f"complexity_k={best_strict['strict']['complexity_k']}"
        ),
        "",
        "## Exploratory (Assignment Pool from QW-1943)",
        f"- exploratory_pass_exists: {expl_pass_exists}",
        (
            "- best exploratory operator: "
            f"{best_expl['operator']['op_id']} @ {best_expl['exploratory_best_scenario']['scenario']} | "
            f"mean/max rel%={best_expl['exploratory_best_scenario']['mean_rel_err_pct']:.3f}/"
            f"{best_expl['exploratory_best_scenario']['max_rel_err_pct']:.3f}"
        ),
        "",
        "## Coupling Completeness Audit",
        f"- legacy QW-1939 feature coverage: 1/{len(CANONICAL_COUPLINGS)} (amplitude only)",
        f"- tested union coverage: {len(tested_feature_union)}/{len(CANONICAL_COUPLINGS)}",
        f"- best strict operator feature coverage: {len(best_strict['operator']['features'])}/{len(CANONICAL_COUPLINGS)}",
        f"- never tested canonical features: {never_tested if never_tested else 'NONE'}",
        f"- best strict missing canonical features: {best_missing if best_missing else 'NONE'}",
        "",
        "## Pareto Frontier (Strict mean error vs complexity)",
    ]
    for row in frontier_short:
        lines.append(
            f"- {row['op_id']}: mean/max rel%={row['strict_mean_rel_err_pct']:.3f}/{row['strict_max_rel_err_pct']:.3f}, "
            f"k={row['complexity_k']}"
        )
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1947] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1947] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1947] verdict={verdict}")


if __name__ == "__main__":
    main()
