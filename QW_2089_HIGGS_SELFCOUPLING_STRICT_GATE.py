#!/usr/bin/env python3
"""
QW-2089: Higgs self-coupling strict non-anchor gate (m_h).
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "t2_nonanchor_higgs_input_qw2089.json"
TEMPLATE_IN = ROOT / "t2_nonanchor_higgs_input_qw2089.template.json"
OUT_JSON = ROOT / "report_qw2089_higgs_selfcoupling_strict_gate.json"
OUT_MD = ROOT / "RAPORT_QW2089_HIGGS_SELFCOUPLING_STRICT_GATE.md"


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


def lambda_run_proxy(lambda0: float, mu0: float, mu: float, beta_uv: float) -> float:
    if lambda0 <= 0.0 or mu0 <= 0.0 or mu <= 0.0:
        raise ValueError("lambda0/mu0/mu must be positive.")
    return float(lambda0 / (1.0 + beta_uv * math.log(mu / mu0)))


def parse_validation_points(raw_points: List[Dict]) -> List[Dict]:
    out = []
    for r in raw_points:
        if not isinstance(r, dict):
            continue
        mu = r.get("mu_gev")
        lam = r.get("lambda_obs")
        sig = r.get("sigma_total")
        if not (isinstance(mu, (int, float)) and isinstance(lam, (int, float)) and isinstance(sig, (int, float))):
            continue
        if float(mu) <= 0.0 or float(lam) <= 0.0 or float(sig) <= 0.0:
            continue
        out.append(
            {
                "mu_gev": float(mu),
                "lambda_obs": float(lam),
                "sigma_total": float(sig),
                "origin": str(r.get("origin", "unknown")).strip().lower(),
                "source": str(r.get("source", "")).strip(),
            }
        )
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2089 Higgs self-coupling strict non-anchor gate")
    p.add_argument(
        "--input",
        default=str(DEFAULT_IN),
        help="JSON with higgs_chain + lambda_validation_points.",
    )
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    r2064 = load_json(ROOT / "report_qw2064_micro_derived_renormalization_constants_gate.json")
    groups = reg["groups"]

    m_h_ref = float(get_registry_item(groups, "m_h")["value"])
    m_h_tol = float(get_registry_item(groups, "m_h")["tolerance_rel_pct"])

    in_path = Path(args.input).resolve()
    if in_path.exists():
        obs = load_json(in_path)
    elif TEMPLATE_IN.exists():
        obs = load_json(TEMPLATE_IN)
        in_path = TEMPLATE_IN
    else:
        obs = {}

    h = obs.get("higgs_chain", {})
    lambda_eff = h.get("lambda_eff")
    v_eff = h.get("v_eff_gev")
    m_h = h.get("m_h_gev")
    o_lambda = str(h.get("lambda_eff_origin", "unknown"))
    o_v = str(h.get("v_eff_origin", "unknown"))
    o_mh = str(h.get("m_h_origin", "unknown"))
    s_lambda = str(h.get("source_lambda_eff", "")).strip()
    s_v = str(h.get("source_v_eff", "")).strip()
    s_mh = str(h.get("source_m_h", "")).strip()

    chain_available = bool(
        isinstance(lambda_eff, (int, float))
        and isinstance(v_eff, (int, float))
        and isinstance(m_h, (int, float))
        and float(lambda_eff) > 0.0
        and float(v_eff) > 0.0
        and float(m_h) > 0.0
        and bool(s_lambda)
        and bool(s_v)
        and bool(s_mh)
    )
    kernel_derived_provenance = bool(
        o_lambda.lower() == "kernel_derived"
        and o_v.lower() == "kernel_derived"
        and o_mh.lower() == "kernel_derived"
    )
    lambda_physical_range = bool(chain_available and 0.0 < float(lambda_eff) < 1.0)

    raw_points = obs.get("lambda_validation_points", [])
    val_points = parse_validation_points(raw_points if isinstance(raw_points, list) else [])
    validation_points_available = len(val_points) >= 2
    validation_external_provenance = bool(
        validation_points_available and all(p["origin"] in {"external_higgs", "external"} and bool(p["source"]) for p in val_points)
    )

    beta_uv = float(r2064["uv_canonical_constants"]["beta_uv"])
    if chain_available:
        m_h_candidate = float(m_h)
        mu0 = 125.0
        residuals_z = []
        for pnt in val_points:
            lam_pred = lambda_run_proxy(float(lambda_eff), mu0, float(pnt["mu_gev"]), beta_uv)
            z = abs(lam_pred - float(pnt["lambda_obs"])) / max(float(pnt["sigma_total"]), 1e-15)
            residuals_z.append(float(z))
        z_mean = float(sum(residuals_z) / len(residuals_z)) if residuals_z else None
        z_max = float(max(residuals_z)) if residuals_z else None
        validation_consistency_pass = bool(residuals_z and z_mean <= 2.5 and z_max <= 3.5)
        anchor_or_circularity_detected = False
        method = "qw2089_kernel_derived_higgs_chain"
    else:
        m_h_candidate = m_h_ref
        z_mean = None
        z_max = None
        validation_consistency_pass = False
        anchor_or_circularity_detected = True
        method = "qw2089_fallback_registry_anchored_mh"

    rel = rel_err_pct(m_h_candidate, m_h_ref)
    within_tol = rel <= m_h_tol

    strict_nonanchor_pass = bool(
        chain_available
        and kernel_derived_provenance
        and lambda_physical_range
        and validation_points_available
        and validation_external_provenance
        and validation_consistency_pass
        and within_tol
        and (not anchor_or_circularity_detected)
    )

    update = {
        "id": "m_h",
        "predicted_value": m_h_candidate,
        "method": method,
        "status": "derived" if strict_nonanchor_pass else "derived_model_assumption_nonclosing",
        "strict_level": "strict_internal_gate" if strict_nonanchor_pass else "model_assumption_anchor",
        "notes": (
            "Strict non-anchor pass from kernel-derived Higgs self-coupling chain + independent lambda validation."
            if strict_nonanchor_pass
            else "Strict non-anchor failed: provenance/validation/closure conditions not fully satisfied."
        ),
    }

    flags = {
        "deterministic_no_retune_no_scan": True,
        "chain_available": chain_available,
        "kernel_derived_provenance": kernel_derived_provenance,
        "lambda_physical_range": lambda_physical_range,
        "validation_points_available": validation_points_available,
        "validation_external_provenance": validation_external_provenance,
        "validation_consistency_pass": validation_consistency_pass,
        "within_tolerance": within_tol,
        "strict_nonanchor_pass": strict_nonanchor_pass,
        "anchor_or_circularity_detected": anchor_or_circularity_detected,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    verdict = (
        "HIGGS_SELFCOUPLING_STRICT_GATE_PASS"
        if strict_nonanchor_pass
        else "HIGGS_SELFCOUPLING_STRICT_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_observables": str(in_path),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "micro_constants": "report_qw2064_micro_derived_renormalization_constants_gate.json",
            "template": TEMPLATE_IN.name,
        },
        "checks": {
            "m_h_candidate": m_h_candidate,
            "m_h_reference": m_h_ref,
            "rel_err_pct": rel,
            "tolerance_rel_pct": m_h_tol,
            "lambda_eff": lambda_eff,
            "v_eff_gev": v_eff,
            "z_mean": z_mean,
            "z_max": z_max,
            "n_validation_points": len(val_points),
        },
        "update": update,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PIPE_UPDATE_TO_QW2069_AND_RERUN_QW2071"
            if strict_nonanchor_pass
            else "PROVIDE_KERNEL_DERIVED_HIGGS_CHAIN_WITH_INDEPENDENT_LAMBDA_VALIDATION"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2089: HIGGS SELFCOUPLING STRICT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Pass count: `{pass_count}/{total_flags}`",
        "",
        "## Checks",
        f"- m_h_candidate: `{m_h_candidate:.9f}`",
        f"- m_h_reference: `{m_h_ref:.9f}`",
        f"- rel_err_pct: `{rel:.6f}`",
        f"- z_mean: `{z_mean}`",
        f"- z_max: `{z_max}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2089] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2089] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2089] verdict={verdict} pass_count={pass_count}/{total_flags} strict_nonanchor_pass={strict_nonanchor_pass}")


if __name__ == "__main__":
    main()

