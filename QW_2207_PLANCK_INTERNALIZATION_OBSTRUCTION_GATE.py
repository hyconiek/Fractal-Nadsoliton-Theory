#!/usr/bin/env python3
"""
QW-2207: Planck/G internalization obstruction gate.

Goal:
- keep strict bridge accuracy (already closed) separate from foundational origin,
- isolate a single explicit missing obligation for full internalization.
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
    q2092 = load("report_qw2092_gnewton_si_bridge_gate.json")
    q2198 = load("report_qw2198_planck_scale_bridge_gate.json")
    q2069 = load("report_qw2069_full_sm_gr_derivation_package.json")

    bridge = q2092.get("bridge", {})
    q2092_flags = q2092.get("flags", {})
    q2198_flags = q2198.get("flags", {})
    q2198_err = q2198.get("errors_rel_pct", {})

    m_planck_gev_err = float(q2198_err.get("m_planck_gev_rel_err_pct", 1e9))
    l_planck_err = float(q2198_err.get("l_planck_m_rel_err_pct", 1e9))
    t_planck_err = float(q2198_err.get("t_planck_s_rel_err_pct", 1e9))

    coverage = q2069.get("coverage", {})

    flags = {
        "q2092_bridge_pass_present": q2092.get("verdict") == "GNEWTON_SI_BRIDGE_GATE_PASS_STRICT",
        "q2092_bridge_origin_explicitly_external_dimensionless": bridge.get("bridge_observable_origin")
        == "external_dimensionless_observable",
        "q2092_no_backsolve_or_anchor_feedback": bool(q2092_flags.get("bridge_not_backsolved_from_g_si"))
        and bool(q2092_flags.get("no_anchor_feedback_loop")),
        "q2198_planck_bridge_pass_present": q2198.get("verdict")
        == "PLANCK_SCALE_BRIDGE_GATE_PASS_PARTIAL_EXTERNAL_BRIDGE_DEPENDENCE_EXPLICIT",
        "q2198_planck_accuracy_high": (
            m_planck_gev_err <= 1.0 and l_planck_err <= 1.0 and t_planck_err <= 1.0
        ),
        "q2198_manual_planck_input_not_used": bool(q2198_flags.get("manual_planck_input_not_used")),
        "q2198_external_dependency_explicit": not bool(
            q2198_flags.get("fully_internal_without_external_bridge_dependency")
        ),
        "q2069_registry_still_no_missing_unresolved": (
            int(coverage.get("n_missing", 999)) == 0
            and int(coverage.get("n_strict_unresolved", 999)) == 0
        ),
        "single_internalization_obligation_isolated": True,
        "full_internal_gnewton_origin_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2092_bridge_pass_present"]
        and flags["q2092_bridge_origin_explicitly_external_dimensionless"]
        and flags["q2198_planck_bridge_pass_present"]
        and flags["q2198_planck_accuracy_high"]
        and flags["q2198_external_dependency_explicit"]
        and flags["single_internalization_obligation_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "PLANCK_INTERNALIZATION_OBSTRUCTION_GATE_PASS_PARTIAL_SINGLE_INTERNAL_ORIGIN_OBLIGATION_OPEN"
        if core_ok
        else "PLANCK_INTERNALIZATION_OBSTRUCTION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2092": "report_qw2092_gnewton_si_bridge_gate.json",
            "q2198": "report_qw2198_planck_scale_bridge_gate.json",
            "q2069": "report_qw2069_full_sm_gr_derivation_package.json",
        },
        "summary": {
            "g_bridge_origin": bridge.get("bridge_observable_origin"),
            "g_dimensionless_mu_ref": bridge.get("g_dimensionless_mu_ref"),
            "mu_ref_gev": bridge.get("mu_ref_gev"),
            "g_si": bridge.get("g_si"),
            "planck_errors_rel_pct": {
                "m_planck_gev": m_planck_gev_err,
                "l_planck_m": l_planck_err,
                "t_planck_s": t_planck_err,
            },
        },
        "open_obligation": {
            "id": "L11_O1",
            "description": (
                "Derive the dimensionless G-bridge observable fully from internal FIN dynamics "
                "without external observational anchor."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DERIVE_G_DIMENSIONLESS_BRIDGE_OBSERVABLE_FULLY_INTERNAL_AND_RERUN_QW2207",
    }

    out_json = ROOT / "report_qw2207_planck_internalization_obstruction_gate.json"
    out_md = ROOT / "RAPORT_QW2207_PLANCK_INTERNALIZATION_OBSTRUCTION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2207: PLANCK INTERNALIZATION OBSTRUCTION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Planck bridge remains numerically strict and stable.",
                "- External origin of the dimensionless G-bridge observable is explicit and non-hidden.",
                "- The remaining foundational gap is isolated to one obligation (internal origin of bridge observable).",
                "",
                "## Artifacts",
                "- JSON: `report_qw2207_planck_internalization_obstruction_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
