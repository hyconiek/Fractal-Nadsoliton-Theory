#!/usr/bin/env python3
"""
QW-2098: EW secondary non-anchor closure gate.

Purpose:
- propagate strict non-anchor EW chain (QW-2085/QW-2086) into secondary
  parameters:
  1) v_higgs,
  2) sin2_theta_w_mz,
  3) m_w,
  4) alpha_em_inv_mz,
- explicitly mark target-miss parameters without promoting false closure.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "t1_nonanchor_observables_input_qw2085_2086.json"
OUT_JSON = ROOT / "report_qw2098_ew_secondary_nonanchor_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def get_registry_item(groups: Dict, pid: str) -> Dict:
    for _, items in groups.items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry parameter: {pid}")


def rel_err_pct(pred: float, ref: float) -> float:
    return abs(pred - ref) / max(abs(ref), 1e-30) * 100.0


def make_update(
    pid: str,
    value: float,
    method: str,
    strict_chain_ready: bool,
    within_tol: bool,
    note_pass: str,
    note_fail: str,
) -> Dict:
    if strict_chain_ready:
        if within_tol:
            status = "derived"
        else:
            status = "derived_strict_target_miss"
        strict_level = "strict_internal_gate"
        notes = note_pass if within_tol else note_fail
    else:
        status = "derived_nofit_anchor_dependent"
        strict_level = "physical_relation_anchor_dependent"
        notes = "Strict non-anchor EW chain unavailable; entry not closure-ready."

    return {
        "id": pid,
        "predicted_value": float(value),
        "method": method,
        "status": status,
        "strict_level": strict_level,
        "notes": notes,
    }


def main() -> None:
    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    r2085 = load_json(ROOT / "report_qw2085_gf_nonanchor_lifetime_gate.json")
    r2086 = load_json(ROOT / "report_qw2086_mz_nonanchor_ew_pole_gate.json")

    input_path = DEFAULT_INPUT.resolve()
    obs = load_json(input_path) if input_path.exists() else {}
    mz_chain = obs.get("mz_pole_chain", {})

    groups = reg["groups"]
    g_f_item = get_registry_item(groups, "g_f")
    v_item = get_registry_item(groups, "v_higgs")
    mw_item = get_registry_item(groups, "m_w")
    sin2_item = get_registry_item(groups, "sin2_theta_w_mz")
    alpha_inv_item = get_registry_item(groups, "alpha_em_inv_mz")

    # Upstream strict gates must pass first.
    strict_2085 = bool(r2085.get("flags", {}).get("strict_nonanchor_pass", False))
    strict_2086 = bool(r2086.get("flags", {}).get("strict_nonanchor_pass", False))

    # Pull deterministic EW observables from non-anchor input chain.
    mw_in = mz_chain.get("mw_pole_gev")
    sin2_in = mz_chain.get("sin2_theta_w_eff")
    delta_r = mz_chain.get("delta_r_full")
    mw_origin = str(mz_chain.get("mw_pole_origin", "")).strip().lower()
    sin2_origin = str(mz_chain.get("sin2_theta_w_eff_origin", "")).strip().lower()
    dr_origin = str(mz_chain.get("delta_r_full_origin", "")).strip().lower()
    src_mw = str(mz_chain.get("source_mw_pole", "")).strip()
    src_sin2 = str(mz_chain.get("source_sin2_theta_w_eff", "")).strip()
    src_dr = str(mz_chain.get("source_delta_r_full", "")).strip()

    chain_inputs_available = bool(
        isinstance(mw_in, (int, float))
        and isinstance(sin2_in, (int, float))
        and isinstance(delta_r, (int, float))
        and float(mw_in) > 0.0
        and 0.0 < float(sin2_in) < 1.0
        and abs(float(delta_r)) < 0.5
        and bool(src_mw)
        and bool(src_sin2)
        and bool(src_dr)
    )
    chain_provenance_kernel = bool(
        mw_origin == "kernel_derived"
        and sin2_origin == "kernel_derived"
        and dr_origin == "kernel_derived"
    )

    g_f_pred = r2085.get("update", {}).get("predicted_value")
    g_f_available = bool(isinstance(g_f_pred, (int, float)) and float(g_f_pred) > 0.0)

    strict_chain_ready = bool(
        strict_2085
        and strict_2086
        and g_f_available
        and chain_inputs_available
        and chain_provenance_kernel
    )

    # Deterministic secondary predictions.
    if g_f_available:
        v_higgs_pred = 1.0 / math.sqrt(math.sqrt(2.0) * float(g_f_pred))
    else:
        v_higgs_pred = float(v_item["value"])

    if chain_inputs_available:
        m_w_pred = float(mw_in)
        sin2_pred = float(sin2_in)
        # On-shell one-loop relation rearranged for alpha(MZ):
        # G_F = [pi * alpha] / [sqrt(2) m_W^2 s_W^2 (1 - Delta r)]
        # => alpha = sqrt(2) G_F m_W^2 s_W^2 (1 - Delta r) / pi
        if g_f_available:
            alpha_mz = (
                math.sqrt(2.0)
                * float(g_f_pred)
                * m_w_pred
                * m_w_pred
                * sin2_pred
                * (1.0 - float(delta_r))
                / math.pi
            )
            alpha_inv_pred = 1.0 / max(alpha_mz, 1e-30)
        else:
            alpha_inv_pred = float(alpha_inv_item["value"])
    else:
        m_w_pred = float(mw_item["value"])
        sin2_pred = float(sin2_item["value"])
        alpha_inv_pred = float(alpha_inv_item["value"])

    # Tolerance checks.
    v_rel = rel_err_pct(v_higgs_pred, float(v_item["value"]))
    mw_rel = rel_err_pct(m_w_pred, float(mw_item["value"]))
    sin2_rel = rel_err_pct(sin2_pred, float(sin2_item["value"]))
    alpha_rel = rel_err_pct(alpha_inv_pred, float(alpha_inv_item["value"]))

    v_ok = bool(v_rel <= float(v_item["tolerance_rel_pct"]))
    mw_ok = bool(mw_rel <= float(mw_item["tolerance_rel_pct"]))
    sin2_ok = bool(sin2_rel <= float(sin2_item["tolerance_rel_pct"]))
    alpha_ok = bool(alpha_rel <= float(alpha_inv_item["tolerance_rel_pct"]))

    updates: List[Dict] = [
        make_update(
            "v_higgs",
            v_higgs_pred,
            "qw2098_v_from_gf_nonanchor_chain",
            strict_chain_ready,
            v_ok,
            "Derived from strict non-anchor G_F chain.",
            "Strict non-anchor chain available, but tolerance miss.",
        ),
        make_update(
            "m_w",
            m_w_pred,
            "qw2098_mw_from_nonanchor_ew_pole_chain",
            strict_chain_ready,
            mw_ok,
            "Derived from strict non-anchor EW-pole chain input.",
            "Strict non-anchor chain available, but tolerance miss.",
        ),
        make_update(
            "sin2_theta_w_mz",
            sin2_pred,
            "qw2098_sin2_from_nonanchor_ew_pole_chain",
            strict_chain_ready,
            sin2_ok,
            "Derived from strict non-anchor EW-pole chain input.",
            "Strict non-anchor chain available, but tolerance miss.",
        ),
        make_update(
            "alpha_em_inv_mz",
            alpha_inv_pred,
            "qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r",
            strict_chain_ready,
            alpha_ok,
            "Derived from strict non-anchor on-shell EW relation.",
            "Strict non-anchor on-shell relation yields target miss at current precision.",
        ),
    ]

    flags = {
        "upstream_qw2085_strict_pass": strict_2085,
        "upstream_qw2086_strict_pass": strict_2086,
        "g_f_available": g_f_available,
        "chain_inputs_available": chain_inputs_available,
        "chain_provenance_kernel": chain_provenance_kernel,
        "strict_chain_ready": strict_chain_ready,
        "v_higgs_within_tolerance": v_ok,
        "m_w_within_tolerance": mw_ok,
        "sin2_theta_w_mz_within_tolerance": sin2_ok,
        "alpha_em_inv_mz_within_tolerance": alpha_ok,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    if strict_chain_ready and all([v_ok, mw_ok, sin2_ok, alpha_ok]):
        verdict = "EW_SECONDARY_NONANCHOR_CLOSURE_GATE_PASS_STRICT"
    elif strict_chain_ready:
        verdict = "EW_SECONDARY_NONANCHOR_CLOSURE_GATE_TARGET_MISS"
    else:
        verdict = "EW_SECONDARY_NONANCHOR_CLOSURE_GATE_FAIL_CHAIN_NOT_READY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "qw2085": "report_qw2085_gf_nonanchor_lifetime_gate.json",
            "qw2086": "report_qw2086_mz_nonanchor_ew_pole_gate.json",
            "input_observables": str(input_path) if input_path.exists() else None,
        },
        "derived": {
            "g_f_pred": float(g_f_pred) if g_f_available else None,
            "v_higgs_pred": float(v_higgs_pred),
            "m_w_pred": float(m_w_pred),
            "sin2_theta_w_mz_pred": float(sin2_pred),
            "alpha_em_inv_mz_pred": float(alpha_inv_pred),
            "delta_r_used": float(delta_r) if chain_inputs_available else None,
        },
        "checks": {
            "v_higgs_rel_err_pct": float(v_rel),
            "m_w_rel_err_pct": float(mw_rel),
            "sin2_theta_w_mz_rel_err_pct": float(sin2_rel),
            "alpha_em_inv_mz_rel_err_pct": float(alpha_rel),
            "v_higgs_tolerance_rel_pct": float(v_item["tolerance_rel_pct"]),
            "m_w_tolerance_rel_pct": float(mw_item["tolerance_rel_pct"]),
            "sin2_theta_w_mz_tolerance_rel_pct": float(sin2_item["tolerance_rel_pct"]),
            "alpha_em_inv_mz_tolerance_rel_pct": float(alpha_inv_item["tolerance_rel_pct"]),
        },
        "updates": updates,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PIPE_QW2098_UPDATES_TO_QW2069_AND_RERUN_QW2071"
            if strict_chain_ready
            else "MAKE_QW2085_AND_QW2086_STRICT_NONANCHOR_PASS_FIRST"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2098: EW SECONDARY NONANCHOR CLOSURE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Derived",
        f"- v_higgs: {v_higgs_pred:.9f} GeV (rel_err={v_rel:.6f}%)",
        f"- m_w: {m_w_pred:.9f} GeV (rel_err={mw_rel:.6f}%)",
        f"- sin2_theta_w_mz: {sin2_pred:.12f} (rel_err={sin2_rel:.6f}%)",
        f"- alpha_em_inv_mz: {alpha_inv_pred:.9f} (rel_err={alpha_rel:.6f}%)",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2098] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2098] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2098] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
