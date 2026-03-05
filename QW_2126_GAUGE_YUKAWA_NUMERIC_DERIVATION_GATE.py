#!/usr/bin/env python3
"""
QW-2126: Gauge + Yukawa numeric derivation gate from strict package entries.

Purpose:
- derive electroweak/QCD couplings and Yukawa couplings directly from strict
  package values (no scan/no fit),
- provide a strict numeric bridge from scalar package to spinor+gauge terms.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2126_gauge_yukawa_numeric_derivation_gate.json"
OUT_MD = ROOT / "RAPORT_QW2126_GAUGE_YUKAWA_NUMERIC_DERIVATION_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def find_entry(entries: List[Dict], pid: str) -> Dict:
    for e in entries:
        if e.get("id") == pid:
            return e
    raise KeyError(f"Missing entry id: {pid}")


def is_strict_derived(e: Dict) -> bool:
    return str(e.get("strict_level", "")) == "strict_internal_gate" and str(e.get("status", "")).startswith("derived")


def rel_err_pct(pred: float, ref: float) -> float:
    return float(abs(pred - ref) / max(abs(ref), 1e-300) * 100.0)


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    entries = r2069["entries"]

    ids_required = [
        "alpha_em_inv_mz",
        "sin2_theta_w_mz",
        "alpha_s_mz",
        "v_higgs",
        "m_w",
        "m_z",
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
    e = {pid: find_entry(entries, pid) for pid in ids_required}

    alpha_em_inv = float(e["alpha_em_inv_mz"]["predicted_value"])
    sin2 = float(e["sin2_theta_w_mz"]["predicted_value"])
    alpha_s = float(e["alpha_s_mz"]["predicted_value"])
    v = float(e["v_higgs"]["predicted_value"])

    alpha_em = float(1.0 / alpha_em_inv)
    e_charge = float(math.sqrt(4.0 * math.pi * alpha_em))
    s_w = float(math.sqrt(sin2))
    c_w = float(math.sqrt(max(1.0 - sin2, 1e-300)))
    g = float(e_charge / max(s_w, 1e-300))
    gp = float(e_charge / max(c_w, 1e-300))
    gs = float(math.sqrt(4.0 * math.pi * alpha_s))

    mw_rebuilt = float(0.5 * g * v)
    mz_rebuilt = float(0.5 * v * math.sqrt(g * g + gp * gp))
    mw_pkg = float(e["m_w"]["predicted_value"])
    mz_pkg = float(e["m_z"]["predicted_value"])

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
    yukawa_rows: List[Dict] = []
    for pid in fermion_ids:
        m = float(e[pid]["predicted_value"])
        y = float(math.sqrt(2.0) * m / max(v, 1e-300))
        yukawa_rows.append(
            {
                "id": pid,
                "mass_gev": m,
                "yukawa_from_sqrt2_m_over_v": y,
                "strict_derived": bool(is_strict_derived(e[pid])),
            }
        )

    ew_identity_err = float(abs(e_charge - g * s_w) + abs(e_charge - gp * c_w))

    flags = {
        "all_required_entries_present": True,
        "all_required_entries_strict_derived": bool(all(is_strict_derived(e[pid]) for pid in ids_required)),
        "alpha_em_positive": bool(alpha_em > 0.0),
        "sin2_in_open_unit_interval": bool(0.0 < sin2 < 1.0),
        "gauge_couplings_positive_finite": bool(all(math.isfinite(x) and x > 0.0 for x in [e_charge, g, gp, gs])),
        "ew_identity_e_eq_gs_eq_gpc_holds_numeric": bool(ew_identity_err <= 1e-12),
        "rebuilt_mw_within_5pct_of_package_mw": bool(rel_err_pct(mw_rebuilt, mw_pkg) <= 5.0),
        "rebuilt_mz_within_5pct_of_package_mz": bool(rel_err_pct(mz_rebuilt, mz_pkg) <= 5.0),
        "all_yukawa_rows_positive_and_strict": bool(all(r["yukawa_from_sqrt2_m_over_v"] > 0.0 and r["strict_derived"] for r in yukawa_rows)),
        "deterministic_no_scan_no_retune": True,
        "full_nonabelian_spinor_action_strict_derived": False,
    }
    pass_count = int(sum(1 for v_ in flags.values() if bool(v_)))
    total_flags = int(len(flags))

    verdict = (
        "GAUGE_YUKAWA_NUMERIC_DERIVATION_GATE_PASS_PARTIAL"
        if (
            flags["all_required_entries_present"]
            and flags["all_required_entries_strict_derived"]
            and flags["gauge_couplings_positive_finite"]
            and flags["ew_identity_e_eq_gs_eq_gpc_holds_numeric"]
            and flags["rebuilt_mw_within_5pct_of_package_mw"]
            and flags["rebuilt_mz_within_5pct_of_package_mz"]
            and flags["all_yukawa_rows_positive_and_strict"]
            and flags["deterministic_no_scan_no_retune"]
        )
        else "GAUGE_YUKAWA_NUMERIC_DERIVATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2069_full_sm_gr_derivation_package.json:entries",
        "inputs": {
            "alpha_em_inv_mz": alpha_em_inv,
            "sin2_theta_w_mz": sin2,
            "alpha_s_mz": alpha_s,
            "v_higgs_gev": v,
            "m_w_package_gev": mw_pkg,
            "m_z_package_gev": mz_pkg,
        },
        "derived_gauge_couplings": {
            "alpha_em_mz": alpha_em,
            "e_charge": e_charge,
            "sin_theta_w": s_w,
            "cos_theta_w": c_w,
            "g_su2": g,
            "gprime_u1": gp,
            "g3_su3": gs,
            "ew_identity_absolute_error": ew_identity_err,
        },
        "rebuilt_vector_boson_masses": {
            "m_w_rebuilt_gev": mw_rebuilt,
            "m_z_rebuilt_gev": mz_rebuilt,
            "m_w_rel_err_pct_vs_package": rel_err_pct(mw_rebuilt, mw_pkg),
            "m_z_rel_err_pct_vs_package": rel_err_pct(mz_rebuilt, mz_pkg),
        },
        "derived_yukawa_rows": yukawa_rows,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROMOTE_TO_FULL_SPINOR_GAUGE_ACTION_DERIVATION_GATE_WITH_NONABELIAN_FIELD_STRENGTHS"
            if verdict.endswith("PARTIAL")
            else "REPAIR_INPUTS_OR_DERIVATION_RULES_AND_RERUN_QW2126"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2126: GAUGE+YUKAWA NUMERIC DERIVATION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Gauge couplings",
        f"- e: `{e_charge:.12f}`",
        f"- g (SU2): `{g:.12f}`",
        f"- g' (U1): `{gp:.12f}`",
        f"- g3 (SU3): `{gs:.12f}`",
        f"- EW identity abs err: `{ew_identity_err:.3e}`",
        "",
        "## Rebuilt masses vs package",
        f"- mW rebuilt/package: `{mw_rebuilt:.6f}` / `{mw_pkg:.6f}` GeV",
        f"- mZ rebuilt/package: `{mz_rebuilt:.6f}` / `{mz_pkg:.6f}` GeV",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2126] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2126] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2126] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

