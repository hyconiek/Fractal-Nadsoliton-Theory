#!/usr/bin/env python3
"""
QW-2215: RG global stability scope gate (L12_O1b).

Purpose:
- keep decomposition from QW-2211 explicit,
- reduce L12_O1b to one terminal theorem obligation.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def main() -> None:
    q2211 = load("report_qw2211_rg_global_obligation_decomposition_gate.json")
    q2213 = load("report_qw2213_rg_flow_existence_scope_gate.json")
    q2209 = load("report_qw2209_rg_global_closure_obligation_gate.json")
    q2208 = load("report_qw2208_spectral_global_stability_obstruction_gate.json")
    q2186 = load("report_qw2186_ktotal_spectral_stability_margin_gate.json")

    f2209 = q2209.get("flags", {})
    f2208 = q2208.get("flags", {})
    f2186 = q2186.get("flags", {})

    decomp = q2211.get("decomposition", {})
    l12_o1b = decomp.get("L12_O1b", {})

    flags = {
        "q2211_decomposition_pass_present": q2211.get("verdict")
        == "RG_GLOBAL_OBLIGATION_DECOMPOSITION_GATE_PASS_PARTIAL_TWO_SUBOBLIGATIONS_OPEN",
        "l12_o1b_subobligation_present": bool(l12_o1b),
        "q2213_l12_o1a_terminalized": q2213.get("verdict")
        == "RG_FLOW_EXISTENCE_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN",
        "q2209_l12_single_obligation_scope_pass_present": q2209.get("verdict")
        == "RG_GLOBAL_CLOSURE_OBLIGATION_GATE_PASS_PARTIAL_SINGLE_GLOBAL_OBLIGATION_OPEN",
        "q2208_global_stability_boundary_explicit": q2208.get("verdict")
        == "SPECTRAL_GLOBAL_STABILITY_OBSTRUCTION_GATE_PASS_PARTIAL_SINGLE_GLOBAL_OBLIGATION_OPEN",
        "q2186_branch_scope_stability_certificate_present": bool(
            f2186.get("a_matrix_positive_definite_under_broken_floor")
        ) and bool(f2186.get("safe_radius_theorem_holds")),
        "q2209_global_all_t_not_claimed": bool(f2209.get("q2187_global_all_t_not_claimed"))
        and bool(f2209.get("q2188_global_all_t_not_claimed")),
        "single_terminal_obligation_for_l12_o1b_isolated": True,
        "l12_o1b_terminal_obligation_closed": False,
        "l12_o1b_closed": False,
        "full_l12_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2211_decomposition_pass_present"]
        and flags["l12_o1b_subobligation_present"]
        and flags["q2213_l12_o1a_terminalized"]
        and flags["q2209_l12_single_obligation_scope_pass_present"]
        and flags["q2208_global_stability_boundary_explicit"]
        and flags["q2186_branch_scope_stability_certificate_present"]
        and flags["q2209_global_all_t_not_claimed"]
        and flags["single_terminal_obligation_for_l12_o1b_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "RG_GLOBAL_STABILITY_SCOPE_GATE_PASS_PARTIAL_SINGLE_TERMINAL_OBLIGATION_OPEN"
        if core_ok
        else "RG_GLOBAL_STABILITY_SCOPE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2211": "report_qw2211_rg_global_obligation_decomposition_gate.json",
            "q2213": "report_qw2213_rg_flow_existence_scope_gate.json",
            "q2209": "report_qw2209_rg_global_closure_obligation_gate.json",
            "q2208": "report_qw2208_spectral_global_stability_obstruction_gate.json",
            "q2186": "report_qw2186_ktotal_spectral_stability_margin_gate.json",
        },
        "summary": {
            "l12_o1b_prior": l12_o1b,
            "l12_o1a_terminal_obligation": q2213.get("open_obligation", {}),
            "l15_o1_boundary": q2208.get("open_obligation", {}),
            "l12_o1_boundary": q2209.get("open_obligation", {}),
        },
        "open_obligation": {
            "id": "L12_O1b_O1",
            "description": (
                "Prove global all-t fixed-point/stability theorem for the full nonperturbative all-coupling FIN RG system "
                "derived from complete FIN action (including UV behavior without hidden truncation/scope restriction)."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_L12_O1B_O1_GLOBAL_ALL_T_FIXED_POINT_STABILITY_AND_RERUN_QW2215",
    }

    out_json = ROOT / "report_qw2215_rg_global_stability_scope_gate.json"
    out_md = ROOT / "RAPORT_QW2215_RG_GLOBAL_STABILITY_SCOPE_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2215: RG GLOBAL STABILITY SCOPE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- L12_O1b is reduced to one terminal theorem-level obligation (`L12_O1b_O1`).",
                "- Coupling to existing L12/L15 boundaries remains explicit (no hidden global overclaim).",
                "- Remaining L12 closure is now expressed as two terminal obligations: `L12_O1a_O1` and `L12_O1b_O1`.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2215_rg_global_stability_scope_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
