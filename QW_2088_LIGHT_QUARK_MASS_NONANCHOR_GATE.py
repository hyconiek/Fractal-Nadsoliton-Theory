#!/usr/bin/env python3
"""
QW-2088: Light-quark mass non-anchor strict gate (m_up, m_down, m_strange).
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "t2_nonanchor_light_quark_input_qw2088.json"
TEMPLATE_IN = ROOT / "t2_nonanchor_light_quark_input_qw2088.template.json"
OUT_JSON = ROOT / "report_qw2088_light_quark_mass_nonanchor_gate.json"
OUT_MD = ROOT / "RAPORT_QW2088_LIGHT_QUARK_MASS_NONANCHOR_GATE.md"


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


def _is_pos_num(x) -> bool:
    return isinstance(x, (int, float)) and float(x) > 0.0


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2088 light-quark mass non-anchor strict gate")
    p.add_argument(
        "--input",
        default=str(DEFAULT_IN),
        help="JSON with light_quark_chain + hadronic_crosschecks.",
    )
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    groups = reg["groups"]

    refs = {}
    for pid in ["m_up", "m_down", "m_strange"]:
        refs[pid] = get_registry_item(groups, pid)

    in_path = Path(args.input).resolve()
    if in_path.exists():
        obs = load_json(in_path)
    elif TEMPLATE_IN.exists():
        obs = load_json(TEMPLATE_IN)
        in_path = TEMPLATE_IN
    else:
        obs = {}

    chain = obs.get("light_quark_chain", {})
    mu = chain.get("m_up_gev")
    md = chain.get("m_down_gev")
    ms = chain.get("m_strange_gev")
    o_mu = str(chain.get("m_up_origin", "unknown"))
    o_md = str(chain.get("m_down_origin", "unknown"))
    o_ms = str(chain.get("m_strange_origin", "unknown"))
    s_mu = str(chain.get("source_m_up", "")).strip()
    s_md = str(chain.get("source_m_down", "")).strip()
    s_ms = str(chain.get("source_m_strange", "")).strip()

    chain_available = bool(
        _is_pos_num(mu)
        and _is_pos_num(md)
        and _is_pos_num(ms)
        and bool(s_mu)
        and bool(s_md)
        and bool(s_ms)
    )
    kernel_derived_provenance = bool(
        o_mu.lower() == "kernel_derived"
        and o_md.lower() == "kernel_derived"
        and o_ms.lower() == "kernel_derived"
    )

    had = obs.get("hadronic_crosschecks", {})
    r_md_mu_obs = had.get("ratio_md_mu_obs")
    r_md_mu_sig = had.get("ratio_md_mu_sigma")
    r_ms_md_obs = had.get("ratio_ms_md_obs")
    r_ms_md_sig = had.get("ratio_ms_md_sigma")
    had_origin = str(had.get("origin", "unknown")).strip().lower()
    had_source = str(had.get("source", "")).strip()

    hadronic_crosschecks_available = bool(
        _is_pos_num(r_md_mu_obs)
        and _is_pos_num(r_md_mu_sig)
        and _is_pos_num(r_ms_md_obs)
        and _is_pos_num(r_ms_md_sig)
        and bool(had_source)
    )
    hadronic_external_provenance = bool(had_origin in {"external_hadronic", "external", "lattice_qcd"})

    if chain_available:
        mu_c = float(mu)
        md_c = float(md)
        ms_c = float(ms)
        anchor_or_circularity_detected = False
        method = "qw2088_kernel_derived_light_quark_chain"
    else:
        mu_c = float(refs["m_up"]["value"])
        md_c = float(refs["m_down"]["value"])
        ms_c = float(refs["m_strange"]["value"])
        anchor_or_circularity_detected = True
        method = "qw2088_fallback_registry_anchored_values"

    rel_mu = rel_err_pct(mu_c, float(refs["m_up"]["value"]))
    rel_md = rel_err_pct(md_c, float(refs["m_down"]["value"]))
    rel_ms = rel_err_pct(ms_c, float(refs["m_strange"]["value"]))
    within_tol_mu = rel_mu <= float(refs["m_up"]["tolerance_rel_pct"])
    within_tol_md = rel_md <= float(refs["m_down"]["tolerance_rel_pct"])
    within_tol_ms = rel_ms <= float(refs["m_strange"]["tolerance_rel_pct"])
    within_tol_all = bool(within_tol_mu and within_tol_md and within_tol_ms)

    if chain_available and hadronic_crosschecks_available:
        r_md_mu_pred = md_c / max(mu_c, 1e-15)
        r_ms_md_pred = ms_c / max(md_c, 1e-15)
        z_md_mu = abs(r_md_mu_pred - float(r_md_mu_obs)) / max(float(r_md_mu_sig), 1e-15)
        z_ms_md = abs(r_ms_md_pred - float(r_ms_md_obs)) / max(float(r_ms_md_sig), 1e-15)
        z_mean = 0.5 * (z_md_mu + z_ms_md)
        z_max = max(z_md_mu, z_ms_md)
        hadronic_consistency_pass = bool(z_mean <= 2.0 and z_max <= 3.0)
    else:
        r_md_mu_pred = None
        r_ms_md_pred = None
        z_md_mu = None
        z_ms_md = None
        z_mean = None
        z_max = None
        hadronic_consistency_pass = False

    strict_nonanchor_pass = bool(
        chain_available
        and kernel_derived_provenance
        and hadronic_crosschecks_available
        and hadronic_external_provenance
        and hadronic_consistency_pass
        and within_tol_all
        and (not anchor_or_circularity_detected)
    )

    strict_level = "strict_internal_gate" if strict_nonanchor_pass else "physical_relation_anchor_dependent"
    status = "derived" if strict_nonanchor_pass else "derived_nofit_anchor_dependent"
    notes = (
        "Strict non-anchor pass from kernel-derived light-quark chain + independent hadronic checks."
        if strict_nonanchor_pass
        else "Strict non-anchor failed: provenance/hadronic consistency/availability constraints not satisfied."
    )

    updates: List[Dict] = [
        {"id": "m_up", "predicted_value": mu_c, "method": method, "status": status, "strict_level": strict_level, "notes": notes},
        {"id": "m_down", "predicted_value": md_c, "method": method, "status": status, "strict_level": strict_level, "notes": notes},
        {"id": "m_strange", "predicted_value": ms_c, "method": method, "status": status, "strict_level": strict_level, "notes": notes},
    ]

    flags = {
        "deterministic_no_retune_no_scan": True,
        "chain_available": chain_available,
        "kernel_derived_provenance": kernel_derived_provenance,
        "hadronic_crosschecks_available": hadronic_crosschecks_available,
        "hadronic_external_provenance": hadronic_external_provenance,
        "hadronic_consistency_pass": hadronic_consistency_pass,
        "within_tolerance_all": within_tol_all,
        "strict_nonanchor_pass": strict_nonanchor_pass,
        "anchor_or_circularity_detected": anchor_or_circularity_detected,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    verdict = (
        "LIGHT_QUARK_MASS_NONANCHOR_GATE_PASS"
        if strict_nonanchor_pass
        else "LIGHT_QUARK_MASS_NONANCHOR_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_observables": str(in_path),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "template": TEMPLATE_IN.name,
        },
        "checks": {
            "m_up": {"candidate": mu_c, "reference": float(refs["m_up"]["value"]), "rel_err_pct": rel_mu, "tolerance_rel_pct": float(refs["m_up"]["tolerance_rel_pct"])},
            "m_down": {"candidate": md_c, "reference": float(refs["m_down"]["value"]), "rel_err_pct": rel_md, "tolerance_rel_pct": float(refs["m_down"]["tolerance_rel_pct"])},
            "m_strange": {"candidate": ms_c, "reference": float(refs["m_strange"]["value"]), "rel_err_pct": rel_ms, "tolerance_rel_pct": float(refs["m_strange"]["tolerance_rel_pct"])},
            "hadronic": {
                "ratio_md_mu_pred": r_md_mu_pred,
                "ratio_md_mu_obs": r_md_mu_obs,
                "ratio_ms_md_pred": r_ms_md_pred,
                "ratio_ms_md_obs": r_ms_md_obs,
                "z_md_mu": z_md_mu,
                "z_ms_md": z_ms_md,
                "z_mean": z_mean,
                "z_max": z_max,
                "origin": had_origin,
                "source": had_source,
            },
        },
        "updates": updates,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PIPE_UPDATES_TO_QW2069_AND_RERUN_QW2071"
            if strict_nonanchor_pass
            else "PROVIDE_KERNEL_DERIVED_LIGHT_QUARK_CHAIN_WITH_INDEPENDENT_HADRONIC_CHECKS"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2088: LIGHT QUARK MASS NONANCHOR GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Pass count: `{pass_count}/{total_flags}`",
        "",
        "## Checks",
        f"- m_up rel_err_pct: `{rel_mu:.6f}`",
        f"- m_down rel_err_pct: `{rel_md:.6f}`",
        f"- m_strange rel_err_pct: `{rel_ms:.6f}`",
        f"- hadronic z_mean: `{z_mean}`",
        f"- hadronic z_max: `{z_max}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2088] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2088] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2088] verdict={verdict} pass_count={pass_count}/{total_flags} strict_nonanchor_pass={strict_nonanchor_pass}")


if __name__ == "__main__":
    main()

