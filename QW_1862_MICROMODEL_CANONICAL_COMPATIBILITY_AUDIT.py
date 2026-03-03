#!/usr/bin/env python3
"""
QW-1862: Micromodel compatibility audit versus reconstructed canonical kernel.

Checks whether micromodel estimates (QW-1739..1744) are numerically compatible
with the canonical tuple reconstructed in QW-1861.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1862_micromodel_canonical_compatibility_audit.json"
OUT_MD = ROOT / "RAPORT_QW1862_MICROMODEL_CANONICAL_COMPATIBILITY_AUDIT.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def circular_diff(a: float, b: float) -> float:
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(d)


def safe_scale(v: float, floor: float) -> float:
    if v is None or not math.isfinite(v) or v <= 0:
        return floor
    return max(v, floor)


def ci_to_sigma(ci: List[float], floor: float) -> float:
    if not isinstance(ci, list) or len(ci) != 2:
        return floor
    try:
        lo = float(ci[0])
        hi = float(ci[1])
    except Exception:
        return floor
    width = max(0.0, hi - lo)
    return safe_scale(width / 3.92, floor)


def param_score(delta: float, sigma: float) -> Dict[str, float]:
    z = delta / sigma if sigma > 0 else float("inf")
    z_cap = min(z, 20.0)
    score = math.exp(-0.5 * z_cap * z_cap)
    return {"delta": delta, "sigma": sigma, "z": z, "score": score}


def tuple_from_1861(d1861: Dict) -> Dict[str, float]:
    w = d1861.get("canonical_tuple", {})

    phi_map = {
        "pi_over_6": math.pi / 6.0,
        "zero": 0.0,
    }
    beta_map = {
        "0.01": 0.01,
        "0.05": 0.05,
    }
    omega_map = {
        "pi_over_4": math.pi / 4.0,
    }

    phi = phi_map.get(w.get("phi"), math.pi / 6.0)
    beta = beta_map.get(w.get("beta_tors"), 0.01)
    omega = omega_map.get(w.get("omega"), math.pi / 4.0)

    return {
        "omega": omega,
        "phi": phi,
        "beta": beta,
        "source": {
            "omega": "QW-1861 winner" if w.get("omega") in omega_map else "fallback_pi_over_4",
            "phi": "QW-1861 winner" if w.get("phi") in phi_map else "fallback_pi_over_6",
            "beta": "QW-1861 winner" if w.get("beta_tors") in beta_map else "fallback_0.01",
        },
    }


def source_row(name: str, est: Dict[str, float], sig: Dict[str, float], target: Dict[str, float], fit_quality: float, nonboundary: bool) -> Dict:
    s_omega = param_score(abs(est["omega"] - target["omega"]), sig["omega"])
    s_phi = param_score(circular_diff(est["phi"], target["phi"]), sig["phi"])
    s_beta = param_score(abs(est["beta"] - target["beta"]), sig["beta"])

    raw_score = (s_omega["score"] * s_phi["score"] * s_beta["score"]) ** (1.0 / 3.0)
    fit_factor = max(0.0, min(1.0, fit_quality))
    boundary_factor = 1.0 if nonboundary else 0.60
    final_score = raw_score * (0.70 + 0.30 * fit_factor) * boundary_factor

    return {
        "source": name,
        "estimate": est,
        "sigma_used": sig,
        "delta_scores": {
            "omega": s_omega,
            "phi": s_phi,
            "beta": s_beta,
        },
        "fit_quality": fit_quality,
        "nonboundary": nonboundary,
        "raw_score": raw_score,
        "final_score": final_score,
    }


def main() -> None:
    d1861 = read_json("report_qw1861_canonical_kernel_reconstruction_700_1600.json")
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1741 = read_json("report_qw1741_constrained_global_derivation.json")
    d1742 = read_json("report_qw1742_profile_likelihood_identifiability.json")
    d1743 = read_json("report_qw1743_oscillatory_cohort_derivation.json")
    d1744 = read_json("report_qw1744_oscillatory_cohort_identifiability.json")

    target = tuple_from_1861(d1861)

    s1739 = d1739.get("summary", {})
    est_1739 = {
        "omega": float(s1739.get("omega_median", math.nan)),
        "phi": float(s1739.get("phi_circular_mean", math.nan)),
        "beta": float(s1739.get("beta_median", math.nan)),
    }
    sig_1739 = {
        "omega": safe_scale(float(s1739.get("omega_iqr", 0.0)) / 1.349, 0.05),
        "phi": safe_scale(float(s1739.get("phi_circular_std", 0.0)), 0.08),
        "beta": safe_scale(float(s1739.get("beta_iqr", 0.0)) / 1.349, 0.01),
    }

    b1741 = d1741.get("bootstrap_summary", {})
    o1741 = d1741.get("optimum", {})
    est_1741 = {
        "omega": float(o1741.get("omega", math.nan)),
        "phi": float(o1741.get("phi", math.nan)),
        "beta": float(o1741.get("beta", math.nan)),
    }
    sig_1741 = {
        "omega": ci_to_sigma(b1741.get("omega_ci95", []), 0.03),
        "phi": ci_to_sigma(b1741.get("phi_quantile_ci95_raw", []), 0.08),
        "beta": ci_to_sigma(b1741.get("beta_ci95", []), 0.01),
    }

    b1743 = d1743.get("bootstrap_summary", {})
    o1743 = d1743.get("optimum", {})
    est_1743 = {
        "omega": float(o1743.get("omega", math.nan)),
        "phi": float(o1743.get("phi", math.nan)),
        "beta": float(o1743.get("beta", math.nan)),
    }
    sig_1743 = {
        "omega": ci_to_sigma(b1743.get("omega_ci95", []), 0.03),
        "phi": safe_scale(float(b1743.get("phi_circular_std", 0.0)), 0.08),
        "beta": ci_to_sigma(b1743.get("beta_ci95", []), 0.01),
    }

    row_1739 = source_row(
        "QW-1739",
        est_1739,
        sig_1739,
        target,
        fit_quality=float(s1739.get("median_r2", 0.0)),
        nonboundary=bool(d1739.get("pass_flags", {}).get("stability", False)),
    )
    row_1741 = source_row(
        "QW-1741",
        est_1741,
        sig_1741,
        target,
        fit_quality=float(d1741.get("run_metrics", {}).get("median_r2", 0.0)),
        nonboundary=bool(d1741.get("pass_flags", {}).get("nonboundary_solution", False)),
    )
    row_1743 = source_row(
        "QW-1743",
        est_1743,
        sig_1743,
        target,
        fit_quality=float(d1743.get("cohort_metrics", {}).get("median_r2", 0.0)),
        nonboundary=bool(d1743.get("pass_flags", {}).get("nonboundary_solution", False)),
    )

    rows = [row_1739, row_1741, row_1743]

    ident1742_flags = d1742.get("pass_flags", {})
    ident1744_flags = d1744.get("pass_flags", {})

    ident_penalty = 1.0
    if not all(bool(v) for v in ident1742_flags.values()):
        ident_penalty *= 0.65
    if not all(bool(v) for v in ident1744_flags.values()):
        ident_penalty *= 0.75

    cond1742 = float(d1742.get("hessian_local", {}).get("condition_number", float("inf")))
    cond1744 = float(d1744.get("hessian", {}).get("condition_number", float("inf")))

    if cond1742 > 1e8:
        ident_penalty *= 0.85
    if cond1744 > 1e8:
        ident_penalty *= 0.85

    mean_base_score = sum(r["final_score"] for r in rows) / len(rows)
    compatibility_index = mean_base_score * ident_penalty

    strongest_row = max(rows, key=lambda x: x["final_score"])

    if compatibility_index >= 0.70 and all(r["nonboundary"] for r in rows):
        verdict = "MICROMODEL_CANONICAL_COMPATIBILITY_PASS"
    elif compatibility_index >= 0.40:
        verdict = "MICROMODEL_CANONICAL_COMPATIBILITY_PARTIAL"
    else:
        verdict = "MICROMODEL_CANONICAL_COMPATIBILITY_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_tuple": target,
        "sources": rows,
        "identifiability": {
            "qw1742_pass_flags": ident1742_flags,
            "qw1744_pass_flags": ident1744_flags,
            "qw1742_condition_number": cond1742,
            "qw1744_condition_number": cond1744,
            "ident_penalty": ident_penalty,
        },
        "summary": {
            "mean_base_score": mean_base_score,
            "compatibility_index": compatibility_index,
            "strongest_source": strongest_row["source"],
            "strongest_source_score": strongest_row["final_score"],
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1862: MICROMODEL CANONICAL COMPATIBILITY AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Target Tuple",
        f"- omega target: {target['omega']:.6f} ({target['source']['omega']})",
        f"- phi target: {target['phi']:.6f} ({target['source']['phi']})",
        f"- beta target: {target['beta']:.6f} ({target['source']['beta']})",
        "",
        "## Source Compatibility",
    ]

    for r in rows:
        lines.append(
            f"- {r['source']}: final_score={r['final_score']:.4f}, raw_score={r['raw_score']:.4f}, "
            f"fit={r['fit_quality']:.3f}, nonboundary={r['nonboundary']}"
        )
        lines.append(
            f"  omega: delta={r['delta_scores']['omega']['delta']:.4f}, z={r['delta_scores']['omega']['z']:.2f}; "
            f"phi: delta={r['delta_scores']['phi']['delta']:.4f}, z={r['delta_scores']['phi']['z']:.2f}; "
            f"beta: delta={r['delta_scores']['beta']['delta']:.4f}, z={r['delta_scores']['beta']['z']:.2f}"
        )

    lines += [
        "",
        "## Identifiability Penalty",
        f"- QW-1742 condition number: {cond1742:.3e}",
        f"- QW-1744 condition number: {cond1744:.3e}",
        f"- ident_penalty: {ident_penalty:.3f}",
        "",
        "## Summary",
        f"- mean_base_score: {mean_base_score:.4f}",
        f"- compatibility_index: {compatibility_index:.4f}",
        f"- strongest_source: {strongest_row['source']} ({strongest_row['final_score']:.4f})",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1862] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1862] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
