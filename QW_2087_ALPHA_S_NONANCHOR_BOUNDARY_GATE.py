#!/usr/bin/env python3
"""
QW-2087: alpha_s non-anchor boundary gate.

Purpose:
- attempt strict non-anchor derivation of alpha_s(MZ) from an independent boundary
  and validation points,
- enforce provenance checks to block anchor leakage,
- produce deterministic update for package integration.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "t1_nonanchor_alpha_s_input_qw2087.json"
TEMPLATE_IN = ROOT / "t1_nonanchor_alpha_s_input_qw2087.template.json"
OUT_JSON = ROOT / "report_qw2087_alpha_s_nonanchor_boundary_gate.json"
OUT_MD = ROOT / "RAPORT_QW2087_ALPHA_S_NONANCHOR_BOUNDARY_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def get_registry_item(groups: Dict, pid: str) -> Dict:
    for _, items in groups.items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry parameter: {pid}")


def rel_err_pct(pred: float, ref: float) -> float:
    return abs(pred - ref) / max(abs(ref), 1e-15) * 100.0


def active_nf_qcd(mu_gev: float) -> int:
    quark_masses = [0.00216, 0.00467, 0.093, 1.27, 4.18, 173.0]
    nf = sum(1 for m in quark_masses if mu_gev >= m)
    return int(max(3, min(6, nf)))


def run_alpha_s_one_loop(mu0: float, alpha0: float, mu_target: float) -> float:
    if mu0 <= 0.0 or mu_target <= 0.0 or alpha0 <= 0.0:
        raise ValueError("Scales and alpha0 must be > 0.")
    if abs(mu_target - mu0) < 1e-15:
        return float(alpha0)

    thresholds = sorted([1.27, 4.18, 173.0])
    forward = mu_target > mu0
    cuts = [x for x in thresholds if mu0 < x < mu_target] if forward else [x for x in thresholds if mu_target < x < mu0]
    cuts = sorted(cuts)
    if not forward:
        cuts = list(reversed(cuts))

    boundaries = [mu0] + cuts + [mu_target]
    inv_alpha = 1.0 / alpha0

    for i in range(len(boundaries) - 1):
        a = boundaries[i]
        b = boundaries[i + 1]
        mu_mid = math.sqrt(a * b)
        nf = active_nf_qcd(mu_mid)
        beta0 = 11.0 - (2.0 / 3.0) * nf
        inv_alpha = inv_alpha + (beta0 / (2.0 * math.pi)) * math.log(b / a)

    if inv_alpha <= 0.0:
        raise ValueError("Unphysical QCD running: inverse alpha_s <= 0.")
    return float(1.0 / inv_alpha)


def parse_validation_points(raw_points: List[Dict]) -> List[Dict]:
    out = []
    for r in raw_points:
        if not isinstance(r, dict):
            continue
        mu = r.get("mu_gev")
        a = r.get("alpha_s_obs")
        s = r.get("sigma_total")
        if not (isinstance(mu, (int, float)) and isinstance(a, (int, float)) and isinstance(s, (int, float))):
            continue
        if float(mu) <= 0.0 or float(a) <= 0.0 or float(s) <= 0.0:
            continue
        out.append(
            {
                "mu_gev": float(mu),
                "alpha_s_obs": float(a),
                "sigma_total": float(s),
                "origin": str(r.get("origin", "unknown")),
                "source": str(r.get("source", "")),
            }
        )
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2087 alpha_s non-anchor boundary gate")
    p.add_argument(
        "--input",
        default=str(DEFAULT_IN),
        help="JSON with alpha_s boundary and validation points.",
    )
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    qcd_baseline = load_json(ROOT / "report_qw2070_full_radiative_program_baseline.json")["qcd_one_loop_baseline"]
    groups = reg["groups"]

    alpha_s_ref = float(get_registry_item(groups, "alpha_s_mz")["value"])
    alpha_s_tol = float(get_registry_item(groups, "alpha_s_mz")["tolerance_rel_pct"])
    mz_ref = float(get_registry_item(groups, "m_z")["value"])

    in_path = Path(args.input).resolve()
    if in_path.exists():
        obs = load_json(in_path)
    elif TEMPLATE_IN.exists():
        obs = load_json(TEMPLATE_IN)
        in_path = TEMPLATE_IN
    else:
        obs = {}

    boundary = obs.get("alpha_s_boundary", {})
    mu0 = boundary.get("mu0_gev")
    alpha0 = boundary.get("alpha_s_mu0")
    nf0 = boundary.get("n_f_active_at_mu0")
    b_origin = str(boundary.get("origin", "unknown"))
    b_source = str(boundary.get("source", "")).strip()

    boundary_available = bool(
        isinstance(mu0, (int, float))
        and isinstance(alpha0, (int, float))
        and isinstance(nf0, (int, float))
        and float(mu0) > 0.0
        and float(alpha0) > 0.0
        and int(nf0) in {3, 4, 5, 6}
        and bool(b_source)
    )
    boundary_kernel_derived = bool(b_origin.lower() == "kernel_derived")

    raw_points = obs.get("validation_points", [])
    val_points = parse_validation_points(raw_points if isinstance(raw_points, list) else [])
    validation_points_available = len(val_points) >= 3
    validation_kernel_derived = bool(
        validation_points_available
        and all(p["origin"].lower() == "kernel_derived" and bool(p["source"]) for p in val_points)
    )

    fallback_alpha_s_mz = float(qcd_baseline["alpha_s_mu0"])
    fallback_mu0 = float(qcd_baseline["mu0_gev"])
    fallback_anchor_leak = bool(abs(fallback_mu0 - mz_ref) < 1e-12 and abs(fallback_alpha_s_mz - alpha_s_ref) < 1e-15)

    residuals_z: List[float] = []
    if boundary_available:
        alpha_s_mz_candidate = run_alpha_s_one_loop(float(mu0), float(alpha0), float(mz_ref))
        for vp in val_points:
            pred = run_alpha_s_one_loop(float(mu0), float(alpha0), vp["mu_gev"])
            z = abs(pred - vp["alpha_s_obs"]) / max(vp["sigma_total"], 1e-15)
            residuals_z.append(float(z))
        candidate_method = "qw2087_nonanchor_boundary_running"
        anchor_leak_detected = False
    else:
        alpha_s_mz_candidate = fallback_alpha_s_mz
        candidate_method = "qw2087_fallback_from_qw2070_anchored_boundary"
        anchor_leak_detected = fallback_anchor_leak

    if residuals_z:
        z_mean = float(sum(residuals_z) / len(residuals_z))
        z_max = float(max(residuals_z))
        validation_consistency_pass = bool(z_mean <= 2.0 and z_max <= 3.0)
    else:
        z_mean = None
        z_max = None
        validation_consistency_pass = False

    rel = rel_err_pct(alpha_s_mz_candidate, alpha_s_ref)
    within_tol = rel <= alpha_s_tol

    strict_nonanchor_pass = bool(
        boundary_available
        and boundary_kernel_derived
        and validation_points_available
        and validation_kernel_derived
        and validation_consistency_pass
        and within_tol
        and (not anchor_leak_detected)
    )

    update = {
        "id": "alpha_s_mz",
        "predicted_value": alpha_s_mz_candidate,
        "method": candidate_method,
        "status": "derived" if strict_nonanchor_pass else "derived_nofit_anchor_dependent",
        "strict_level": "strict_internal_gate" if strict_nonanchor_pass else "physical_relation_anchor_dependent",
        "notes": (
            "Strict non-anchor pass from kernel-derived alpha_s boundary and validation points."
            if strict_nonanchor_pass
            else "Strict non-anchor failed: boundary/provenance/validation constraints not satisfied."
        ),
    }

    flags = {
        "deterministic_no_retune_no_scan": True,
        "boundary_available": boundary_available,
        "boundary_kernel_derived": boundary_kernel_derived,
        "validation_points_available": validation_points_available,
        "validation_kernel_derived": validation_kernel_derived,
        "validation_consistency_pass": validation_consistency_pass,
        "within_tolerance": within_tol,
        "strict_nonanchor_pass": strict_nonanchor_pass,
        "anchor_leak_detected": anchor_leak_detected,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    verdict = (
        "ALPHA_S_NONANCHOR_BOUNDARY_GATE_PASS"
        if strict_nonanchor_pass
        else "ALPHA_S_NONANCHOR_BOUNDARY_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_observables": str(in_path),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "qcd_baseline": "report_qw2070_full_radiative_program_baseline.json",
            "template": TEMPLATE_IN.name,
        },
        "checks": {
            "alpha_s_mz_candidate": alpha_s_mz_candidate,
            "alpha_s_mz_reference": alpha_s_ref,
            "rel_err_pct": rel,
            "tolerance_rel_pct": alpha_s_tol,
            "boundary_mu0_gev": mu0,
            "boundary_alpha_s_mu0": alpha0,
            "boundary_n_f_active": nf0,
            "boundary_origin": b_origin,
            "boundary_source": b_source,
            "n_validation_points": len(val_points),
            "z_mean": z_mean,
            "z_max": z_max,
        },
        "update": update,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROVIDE_KERNEL_DERIVED_ALPHA_S_BOUNDARY_AND_VALIDATION_POINTS"
            if verdict != "ALPHA_S_NONANCHOR_BOUNDARY_GATE_PASS"
            else "PIPE_UPDATE_TO_QW2069_AND_RERUN_QW2071"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2087: ALPHA_S NONANCHOR BOUNDARY GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Pass count: `{pass_count}/{total_flags}`",
        "",
        "## Checks",
        f"- alpha_s_mz_candidate: `{alpha_s_mz_candidate:.9f}`",
        f"- alpha_s_mz_reference: `{alpha_s_ref:.9f}`",
        f"- rel_err_pct: `{rel:.6f}`",
        f"- tolerance_rel_pct: `{alpha_s_tol:.6f}`",
        f"- n_validation_points: `{len(val_points)}`",
        f"- z_mean: `{z_mean}`",
        f"- z_max: `{z_max}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2087] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2087] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2087] verdict={verdict} pass_count={pass_count}/{total_flags} "
        f"strict_nonanchor_pass={strict_nonanchor_pass}"
    )


if __name__ == "__main__":
    main()

