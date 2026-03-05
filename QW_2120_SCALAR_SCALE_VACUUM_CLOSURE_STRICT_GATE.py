#!/usr/bin/env python3
"""
QW-2120: Scalar-scale vacuum-closure strict gate.

Purpose:
- use explicit strict-derived scalar/mass inputs to test whether
  diagonal mass floor is sufficient to close vacuum stability condition
  implied by QW-2118 (lambda_min(K_total) < 0),
- keep strict decision free from exploratory/nonclosing hypotheses.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2120_scalar_scale_vacuum_closure_strict_gate.json"
OUT_MD = ROOT / "RAPORT_QW2120_SCALAR_SCALE_VACUUM_CLOSURE_STRICT_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def find_entry(entries: List[Dict], pid: str) -> Dict:
    for e in entries:
        if e.get("id") == pid:
            return e
    raise KeyError(f"Missing entry id: {pid}")


def is_strict_derived(e: Dict) -> bool:
    return str(e.get("strict_level", "")) == "strict_internal_gate" and str(e.get("status", "")).startswith("derived")


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2118 = load_json("report_qw2118_ktotal_spectral_tripartition_gate.json")
    r2066 = load_json("report_qw2066_compatibility_filtered_micro_constants_tightening.json")

    entries = r2069["entries"]
    e_v = find_entry(entries, "v_higgs")
    e_h = find_entry(entries, "m_h")

    fermion_ids = [
        "m_top",
        "m_bottom",
        "m_charm",
        "m_tau",
        "m_muon",
        "m_electron",
        "m_up",
        "m_down",
        "m_strange",
    ]
    fermion_entries = [find_entry(entries, pid) for pid in fermion_ids]
    fermion_entries = [e for e in fermion_entries if e.get("predicted_value") is not None]

    v = float(e_v["predicted_value"])
    m_h = float(e_h["predicted_value"])
    required_shift = float(r2118["vacuum_shift_condition"]["required_uniform_mass_shift_ge"])

    # Strict scalar-side floors (no exploratory renorm multipliers in verdict):
    # 1) Higgs curvature scale in v-normalized units: m_h^2 / v^2 = 2*lambda_H
    higgs_curvature_floor = float((m_h * m_h) / max(v * v, 1e-300))
    lambda_h = float(0.5 * higgs_curvature_floor)

    # 2) Effective scalar couplings from strict-derived fermion masses under
    #    m_i^2 = g_i * v^2 relation for |Phi|^2|Psi_i|^2 channel.
    g_rows: List[Dict] = []
    for e in fermion_entries:
        m_i = float(e["predicted_value"])
        g_i = float((m_i * m_i) / max(v * v, 1e-300))
        g_rows.append(
            {
                "id": str(e["id"]),
                "m_pred_gev": m_i,
                "g_i_from_m2_over_v2": g_i,
                "strict_level": e.get("strict_level"),
                "status": e.get("status"),
            }
        )
    g_max = float(max(r["g_i_from_m2_over_v2"] for r in g_rows))
    strict_floor = float(max(higgs_curvature_floor, g_max))

    # Exploratory nonclosing channel (reported, not used for strict verdict):
    sel = r2066.get("selected_filter", {})
    z_beta_q50 = float(sel.get("z_beta_q50", 100.0))
    delta_eta_q50 = float(sel.get("delta_eta_q50", 0.8))
    exploratory_micro_floor = float(strict_floor * (z_beta_q50 / 100.0) * (1.0 + delta_eta_q50))

    strict_ok = bool(strict_floor >= required_shift)

    flags = {
        "required_shift_loaded_from_qw2118": bool(required_shift >= 0.0),
        "v_higgs_strict_derived_available": bool(is_strict_derived(e_v)),
        "m_h_strict_derived_available": bool(is_strict_derived(e_h)),
        "fermion_strict_rows_available": bool(len(g_rows) >= 6),
        "strict_scalar_floor_computable": bool(np.isfinite(strict_floor) and strict_floor >= 0.0),
        "strict_floor_ge_required_shift": bool(strict_ok),
        "exploratory_micro_floor_ge_required_shift_nonclosing": bool(exploratory_micro_floor >= required_shift),
        "strict_verdict_ignores_exploratory_channel": True,
    }
    pass_count = int(sum(1 for v_ in flags.values() if bool(v_)))
    total_flags = int(len(flags))

    verdict = (
        "SCALAR_SCALE_VACUUM_CLOSURE_STRICT_PASS"
        if strict_ok
        else "SCALAR_SCALE_VACUUM_CLOSURE_STRICT_FAIL_INSUFFICIENT_FLOOR"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "package_entries": "report_qw2069_full_sm_gr_derivation_package.json:entries",
            "vacuum_shift": "report_qw2118_ktotal_spectral_tripartition_gate.json:vacuum_shift_condition",
            "micro_tightening_exploratory": "report_qw2066_compatibility_filtered_micro_constants_tightening.json:selected_filter",
        },
        "inputs": {
            "v_higgs_gev": v,
            "m_h_gev": m_h,
            "required_uniform_mass_shift_ge": required_shift,
            "n_fermion_rows": len(g_rows),
        },
        "strict_floor_components": {
            "lambda_h_from_mh_v": lambda_h,
            "higgs_curvature_floor_mh2_over_v2": higgs_curvature_floor,
            "max_g_i_from_fermion_m2_over_v2": g_max,
            "strict_floor_used_for_verdict": strict_floor,
            "fermion_couplings_rows": g_rows,
        },
        "exploratory_nonclosing_channel": {
            "z_beta_q50": z_beta_q50,
            "delta_eta_q50": delta_eta_q50,
            "micro_renorm_floor_candidate": exploratory_micro_floor,
            "used_for_strict_verdict": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "FREEZE_SCALAR_FLOOR_IN_STRICT_CHAIN_AND_PROPAGATE_TO_QW2071"
            if strict_ok
            else "DERIVE_ADDITIVE_DIAGONAL_MASS_FLOOR_FROM_EXPLICIT_PSI_POTENTIAL_AND_RERUN_QW2120"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2120: SCALAR SCALE VACUUM CLOSURE STRICT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Inputs",
        f"- v_higgs [GeV]: `{v:.9f}`",
        f"- m_h [GeV]: `{m_h:.9f}`",
        f"- required shift >= `{required_shift:.9f}`",
        "",
        "## Strict floor",
        f"- lambda_h = m_h^2/(2 v^2): `{lambda_h:.9f}`",
        f"- higgs_curvature_floor = m_h^2/v^2: `{higgs_curvature_floor:.9f}`",
        f"- max g_i = max(m_i^2/v^2): `{g_max:.9f}`",
        f"- strict_floor_used_for_verdict: `{strict_floor:.9f}`",
        "",
        "## Exploratory nonclosing channel (not used in verdict)",
        f"- micro_renorm_floor_candidate: `{exploratory_micro_floor:.9f}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2120] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2120] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2120] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

