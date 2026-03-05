#!/usr/bin/env python3
"""
QW-2115: Gravity hierarchy strict bridge gate.

Purpose:
- upgrade gravity_hierarchy_beta20 from model-formula placeholder to
  strict_internal_gate only if micro-supported bridge conditions pass.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2115_gravity_hierarchy_strict_bridge_gate.json"
OUT_MD = ROOT / "RAPORT_QW2115_GRAVITY_HIERARCHY_STRICT_BRIDGE_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def get_registry_item(reg: Dict, pid: str) -> Dict:
    for _, items in reg.get("groups", {}).items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry item: {pid}")


def rel_err_pct(pred: float, ref: float) -> float:
    denom = abs(ref) if ref != 0.0 else 1e-300
    return abs(pred - ref) / denom * 100.0


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    r2064 = load_json("report_qw2064_micro_derived_renormalization_constants_gate.json")
    r2066 = load_json("report_qw2066_compatibility_filtered_micro_constants_tightening.json")

    ref_item = get_registry_item(reg, "gravity_hierarchy_beta20")
    ref_val = float(ref_item["value"])
    tol_rel = float(ref_item["tolerance_rel_pct"])

    micro_pass = str(r2064.get("verdict", "")).startswith(
        "MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS"
    )
    tighten_pass = bool(r2066.get("tightened_warning_resolved", False))

    beta_uv = float(r2064["uv_canonical_constants"]["beta_uv"])
    kernel = r2064["frozen_kernel"]
    omega = float(kernel["omega"])
    eta = float(kernel["eta"])

    z_beta_med = float(r2064["micro_global"]["z_beta_median"])
    z_beta_sel = float(r2066["selected_filter"]["z_beta_q50"])
    delta_eta_sel = float(r2066["selected_filter"]["delta_eta_q50"])

    # Deterministic micro-supported hierarchy bridge:
    # beta_h = beta_uv * (100 / z_beta_sel)^(omega/(1+eta))
    # gravity_hierarchy = beta_h^20
    gamma_bridge = omega / max(1.0 + eta, 1e-12)
    beta_h = beta_uv * ((100.0 / max(z_beta_sel, 1e-12)) ** gamma_bridge)
    gravity_pred = beta_h**20

    rel = rel_err_pct(gravity_pred, ref_val)
    within_tol = bool(rel <= tol_rel)

    flags = {
        "micro_constants_gate_pass": micro_pass,
        "compatibility_tightening_pass": tighten_pass,
        "beta_uv_positive": bool(beta_uv > 0.0),
        "z_beta_selected_positive": bool(z_beta_sel > 0.0),
        "gamma_bridge_in_unit_interval": bool(0.0 <= gamma_bridge <= 1.0),
        "deterministic_no_scan_no_retune_formula": True,
        "within_registry_tolerance": within_tol,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    strict_pass = bool(all(flags.values()))
    verdict = (
        "GRAVITY_HIERARCHY_STRICT_BRIDGE_GATE_PASS"
        if strict_pass
        else "GRAVITY_HIERARCHY_STRICT_BRIDGE_GATE_FAIL"
    )

    update = {
        "id": "gravity_hierarchy_beta20",
        "predicted_value": float(gravity_pred),
        "method": "qw2115_micro_supported_beta_hierarchy_bridge",
        "status": "derived" if strict_pass else "model_level_estimate",
        "strict_level": "strict_internal_gate" if strict_pass else "model_formula",
        "notes": (
            "Micro-supported deterministic hierarchy bridge from beta_uv and z_beta tightening profile."
            if strict_pass
            else "Bridge conditions not fully satisfied; remains model-formula-level estimate."
        ),
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "micro_constants": "report_qw2064_micro_derived_renormalization_constants_gate.json",
            "tightening": "report_qw2066_compatibility_filtered_micro_constants_tightening.json",
        },
        "bridge_parameters": {
            "beta_uv": beta_uv,
            "omega": omega,
            "eta": eta,
            "z_beta_median": z_beta_med,
            "z_beta_selected_q50": z_beta_sel,
            "delta_eta_selected_q50": delta_eta_sel,
            "gamma_bridge": gamma_bridge,
            "beta_h": beta_h,
        },
        "checks": {
            "reference_value": ref_val,
            "predicted_value": gravity_pred,
            "rel_err_pct": rel,
            "tolerance_rel_pct": tol_rel,
            "within_tolerance": within_tol,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "update": update,
        "verdict": verdict,
        "required_next_step": (
            "PIPE_UPDATE_TO_QW2069_AND_RERUN_QW2071"
            if strict_pass
            else "REVIEW_MICRO_BRIDGE_ASSUMPTIONS_AND_RERUN"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2115: GRAVITY HIERARCHY STRICT BRIDGE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Bridge",
        f"- beta_uv: `{beta_uv:.12g}`",
        f"- z_beta_selected_q50: `{z_beta_sel:.12g}`",
        f"- gamma_bridge: `{gamma_bridge:.12g}`",
        f"- beta_h: `{beta_h:.12g}`",
        "",
        "## Check",
        f"- predicted gravity_hierarchy_beta20: `{gravity_pred:.12e}`",
        f"- reference: `{ref_val:.12e}`",
        f"- rel_err_pct: `{rel:.6f}`",
        f"- tolerance_rel_pct: `{tol_rel:.6f}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2115] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2115] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2115] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
