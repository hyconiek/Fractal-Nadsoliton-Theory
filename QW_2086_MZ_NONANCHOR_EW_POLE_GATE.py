#!/usr/bin/env python3
"""
QW-2086: M_Z non-anchor electroweak pole gate.

Purpose:
- attempt non-anchor strict derivation of M_Z from independent EW-pole observables,
- enforce explicit provenance checks (kernel-derived vs anchor-derived),
- block false strict closure when current chain injects M_Z anchor.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "t1_nonanchor_observables_input_qw2085_2086.json"
TEMPLATE_IN = ROOT / "t1_nonanchor_observables_input_qw2085_2086.template.json"
OUT_JSON = ROOT / "report_qw2086_mz_nonanchor_ew_pole_gate.json"
OUT_MD = ROOT / "RAPORT_QW2086_MZ_NONANCHOR_EW_POLE_GATE.md"


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


def mz_from_mw_sin2(mw_gev: float, sin2_theta_w: float) -> float:
    if mw_gev <= 0.0:
        raise ValueError("mw_gev must be positive.")
    if not (0.0 < sin2_theta_w < 1.0):
        raise ValueError("sin2_theta_w must be in (0,1).")
    return mw_gev / math.sqrt(max(1.0 - sin2_theta_w, 1e-15))


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2086 M_Z non-anchor EW pole gate")
    p.add_argument(
        "--input",
        default=str(DEFAULT_IN),
        help="JSON with mz_pole_chain block (default: t1_nonanchor_observables_input_qw2085_2086.json).",
    )
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    ew = load_json(ROOT / "report_qw2072_ew_yukawa_flavor_radiative_baselines.json")["electroweak_baseline"]
    groups = reg["groups"]

    mz_ref = float(get_registry_item(groups, "m_z")["value"])
    mz_tol = float(get_registry_item(groups, "m_z")["tolerance_rel_pct"])

    in_path = Path(args.input).resolve()
    if in_path.exists():
        obs = load_json(in_path)
    elif TEMPLATE_IN.exists():
        obs = load_json(TEMPLATE_IN)
        in_path = TEMPLATE_IN
    else:
        obs = {}

    mz_chain = obs.get("mz_pole_chain", {})
    mw_pole_gev = mz_chain.get("mw_pole_gev")
    sin2_eff = mz_chain.get("sin2_theta_w_eff")
    delta_r_full = mz_chain.get("delta_r_full")
    mw_origin = str(mz_chain.get("mw_pole_origin", "unknown"))
    sin2_origin = str(mz_chain.get("sin2_theta_w_eff_origin", "unknown"))
    delta_r_origin = str(mz_chain.get("delta_r_full_origin", "unknown"))
    src_mw = str(mz_chain.get("source_mw_pole", "")).strip()
    src_sin2 = str(mz_chain.get("source_sin2_theta_w_eff", "")).strip()
    src_deltar = str(mz_chain.get("source_delta_r_full", "")).strip()

    has_mw = isinstance(mw_pole_gev, (int, float)) and float(mw_pole_gev) > 0.0
    has_sin2 = isinstance(sin2_eff, (int, float)) and 0.0 < float(sin2_eff) < 1.0
    has_deltar = isinstance(delta_r_full, (int, float))
    has_sources = bool(src_mw) and bool(src_sin2) and bool(src_deltar)

    ew_pole_chain_available = bool(has_mw and has_sin2 and has_deltar and has_sources)
    kernel_derived_provenance = bool(
        mw_origin.lower() == "kernel_derived"
        and sin2_origin.lower() == "kernel_derived"
        and delta_r_origin.lower() == "kernel_derived"
    )

    fallback_mz = float(
        ew["m_w_delta_r0_from_alpha_gf_mz"]
        / math.sqrt(max(1.0 - float(ew["sin2_theta_w_input"]), 1e-15))
    )

    if ew_pole_chain_available:
        mz_candidate = float(mz_from_mw_sin2(float(mw_pole_gev), float(sin2_eff)))
        candidate_method = "qw2086_mz_from_independent_mw_and_sin2"
        anchor_leak_detected = False
    else:
        mz_candidate = fallback_mz
        candidate_method = "qw2086_fallback_from_qw2072_anchored_ew_baseline"
        anchor_leak_detected = True

    mz_rel = rel_err_pct(mz_candidate, mz_ref)
    within_tol = mz_rel <= mz_tol

    strict_nonanchor_pass = bool(
        ew_pole_chain_available and kernel_derived_provenance and within_tol and (not anchor_leak_detected)
    )

    update = {
        "id": "m_z",
        "predicted_value": mz_candidate,
        "method": candidate_method,
        "status": "derived" if strict_nonanchor_pass else "derived_nofit_anchor_dependent",
        "strict_level": "strict_internal_gate" if strict_nonanchor_pass else "physical_relation_anchor_dependent",
        "notes": (
            "Strict non-anchor pass from kernel-derived EW-pole chain."
            if strict_nonanchor_pass
            else (
                "Strict non-anchor failed: missing kernel-derived EW-pole inputs."
                if not ew_pole_chain_available
                else "Strict non-anchor failed: provenance is not kernel_derived."
            )
        ),
    }

    flags = {
        "deterministic_no_retune_no_scan": True,
        "ew_pole_chain_available": ew_pole_chain_available,
        "kernel_derived_provenance": kernel_derived_provenance,
        "within_tolerance": within_tol,
        "strict_nonanchor_pass": strict_nonanchor_pass,
        "anchor_leak_detected": anchor_leak_detected,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    verdict = (
        "MZ_NONANCHOR_EW_POLE_GATE_PASS"
        if strict_nonanchor_pass
        else "MZ_NONANCHOR_EW_POLE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_observables": str(in_path),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "ew_baseline": "report_qw2072_ew_yukawa_flavor_radiative_baselines.json",
            "template": TEMPLATE_IN.name,
        },
        "checks": {
            "m_z_candidate": mz_candidate,
            "m_z_reference": mz_ref,
            "rel_err_pct": mz_rel,
            "tolerance_rel_pct": mz_tol,
            "mw_pole_gev": mw_pole_gev,
            "sin2_theta_w_eff": sin2_eff,
            "delta_r_full": delta_r_full,
            "mw_pole_origin": mw_origin,
            "sin2_theta_w_eff_origin": sin2_origin,
            "delta_r_full_origin": delta_r_origin,
            "source_mw_pole": src_mw,
            "source_sin2_theta_w_eff": src_sin2,
            "source_delta_r_full": src_deltar,
        },
        "update": update,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROVIDE_KERNEL_DERIVED_EW_POLE_CHAIN_AND_REPEAT"
            if verdict != "MZ_NONANCHOR_EW_POLE_GATE_PASS"
            else "PIPE_UPDATE_TO_QW2069_AND_RERUN_QW2071"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2086: M_Z NONANCHOR EW POLE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Pass count: `{pass_count}/{total_flags}`",
        "",
        "## Checks",
        f"- m_z_candidate: `{mz_candidate:.9f}`",
        f"- m_z_reference: `{mz_ref:.9f}`",
        f"- rel_err_pct: `{mz_rel:.6f}`",
        f"- tolerance_rel_pct: `{mz_tol:.6f}`",
        f"- ew_pole_chain_available: `{ew_pole_chain_available}`",
        f"- kernel_derived_provenance: `{kernel_derived_provenance}`",
        f"- anchor_leak_detected: `{anchor_leak_detected}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2086] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2086] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2086] verdict={verdict} pass_count={pass_count}/{total_flags} "
        f"strict_nonanchor_pass={strict_nonanchor_pass}"
    )


if __name__ == "__main__":
    main()

