#!/usr/bin/env python3
"""
QW-2085: G_F non-anchor lifetime gate.

Purpose:
- attempt non-anchor strict derivation of G_F from an independent muon-lifetime chain,
- enforce explicit provenance checks (kernel-derived vs anchor-derived),
- prevent false strict promotion when circular/anchored inputs are used.
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
OUT_JSON = ROOT / "report_qw2085_gf_nonanchor_lifetime_gate.json"
OUT_MD = ROOT / "RAPORT_QW2085_GF_NONANCHOR_LIFETIME_GATE.md"


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


def gf_from_tau_mu(tau_mu_s: float, m_mu_gev: float, delta_q: float) -> float:
    # Natural-units conversion: tau[GeV^-1] = tau[s] / hbar[GeV*s]
    hbar_gev_s = 6.582119569e-25
    tau_nat = tau_mu_s / hbar_gev_s
    pref = 192.0 * (math.pi**3)
    denom = tau_nat * (m_mu_gev**5) * (1.0 + delta_q)
    if denom <= 0.0:
        raise ValueError("Invalid denominator in G_F lifetime derivation.")
    return math.sqrt(pref / denom)


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2085 G_F non-anchor lifetime gate")
    p.add_argument(
        "--input",
        default=str(DEFAULT_IN),
        help="JSON with gf_lifetime_chain block (default: t1_nonanchor_observables_input_qw2085_2086.json).",
    )
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    ew = load_json(ROOT / "report_qw2072_ew_yukawa_flavor_radiative_baselines.json")["electroweak_baseline"]
    groups = reg["groups"]

    gf_ref = float(get_registry_item(groups, "g_f")["value"])
    gf_tol = float(get_registry_item(groups, "g_f")["tolerance_rel_pct"])
    m_mu_ref = float(get_registry_item(groups, "m_muon")["value"])

    in_path = Path(args.input).resolve()
    if in_path.exists():
        obs = load_json(in_path)
    elif TEMPLATE_IN.exists():
        obs = load_json(TEMPLATE_IN)
        in_path = TEMPLATE_IN
    else:
        obs = {}

    lf = obs.get("gf_lifetime_chain", {})
    m_mu_in = lf.get("m_mu_gev")
    tau_mu_s = lf.get("tau_mu_s")
    delta_q = lf.get("delta_q")
    m_mu_origin = str(lf.get("m_mu_origin", "unknown"))
    tau_origin = str(lf.get("tau_mu_origin", "unknown"))
    delta_origin = str(lf.get("delta_q_origin", "unknown"))
    src_m_mu = str(lf.get("source_m_mu", "")).strip()
    src_tau = str(lf.get("source_tau_mu", "")).strip()
    src_delta = str(lf.get("source_delta_q", "")).strip()

    has_m_mu = isinstance(m_mu_in, (int, float)) and float(m_mu_in) > 0.0
    has_tau = isinstance(tau_mu_s, (int, float)) and float(tau_mu_s) > 0.0
    has_delta = isinstance(delta_q, (int, float)) and -0.5 < float(delta_q) < 0.5
    has_sources = bool(src_m_mu) and bool(src_tau) and bool(src_delta)

    lifetime_chain_available = bool(has_m_mu and has_tau and has_delta and has_sources)
    kernel_derived_provenance = bool(
        m_mu_origin.lower() == "kernel_derived"
        and tau_origin.lower() == "kernel_derived"
        and delta_origin.lower() == "kernel_derived"
    )

    fallback_gf = float(1.0 / (math.sqrt(2.0) * float(ew["v_from_gf"]) ** 2))

    if lifetime_chain_available:
        gf_candidate = float(gf_from_tau_mu(float(tau_mu_s), float(m_mu_in), float(delta_q)))
        candidate_method = "qw2085_muon_lifetime_chain"
        circularity_detected = False
    else:
        gf_candidate = fallback_gf
        candidate_method = "qw2085_fallback_identity_over_v_from_gf"
        circularity_detected = True

    gf_rel = rel_err_pct(gf_candidate, gf_ref)
    within_tol = gf_rel <= gf_tol

    strict_nonanchor_pass = bool(
        lifetime_chain_available and kernel_derived_provenance and within_tol and (not circularity_detected)
    )

    update = {
        "id": "g_f",
        "predicted_value": gf_candidate,
        "method": candidate_method,
        "status": "derived" if strict_nonanchor_pass else "derived_nofit_anchor_dependent",
        "strict_level": "strict_internal_gate" if strict_nonanchor_pass else "physical_relation_anchor_dependent",
        "notes": (
            "Strict non-anchor pass from kernel-derived muon-lifetime chain."
            if strict_nonanchor_pass
            else (
                "Strict non-anchor failed: missing kernel-derived lifetime inputs."
                if not lifetime_chain_available
                else "Strict non-anchor failed: provenance is not kernel_derived."
            )
        ),
    }

    flags = {
        "deterministic_no_retune_no_scan": True,
        "lifetime_chain_available": lifetime_chain_available,
        "kernel_derived_provenance": kernel_derived_provenance,
        "within_tolerance": within_tol,
        "strict_nonanchor_pass": strict_nonanchor_pass,
        "circularity_detected": circularity_detected,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    verdict = (
        "GF_NONANCHOR_LIFETIME_GATE_PASS"
        if strict_nonanchor_pass
        else "GF_NONANCHOR_LIFETIME_GATE_FAIL"
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
            "g_f_candidate": gf_candidate,
            "g_f_reference": gf_ref,
            "rel_err_pct": gf_rel,
            "tolerance_rel_pct": gf_tol,
            "tau_mu_s": tau_mu_s,
            "delta_q": delta_q,
            "m_mu_gev": m_mu_in,
            "m_mu_reference_gev": m_mu_ref,
            "m_mu_origin": m_mu_origin,
            "tau_mu_origin": tau_origin,
            "delta_q_origin": delta_origin,
            "source_m_mu": src_m_mu,
            "source_tau_mu": src_tau,
            "source_delta_q": src_delta,
        },
        "update": update,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROVIDE_KERNEL_DERIVED_MUON_LIFETIME_CHAIN_AND_REPEAT"
            if verdict != "GF_NONANCHOR_LIFETIME_GATE_PASS"
            else "PIPE_UPDATE_TO_QW2069_AND_RERUN_QW2071"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2085: G_F NONANCHOR LIFETIME GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Pass count: `{pass_count}/{total_flags}`",
        "",
        "## Checks",
        f"- g_f_candidate: `{gf_candidate:.12e}`",
        f"- g_f_reference: `{gf_ref:.12e}`",
        f"- rel_err_pct: `{gf_rel:.6f}`",
        f"- tolerance_rel_pct: `{gf_tol:.6f}`",
        f"- lifetime_chain_available: `{lifetime_chain_available}`",
        f"- kernel_derived_provenance: `{kernel_derived_provenance}`",
        f"- circularity_detected: `{circularity_detected}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2085] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2085] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2085] verdict={verdict} pass_count={pass_count}/{total_flags} "
        f"strict_nonanchor_pass={strict_nonanchor_pass}"
    )


if __name__ == "__main__":
    main()
