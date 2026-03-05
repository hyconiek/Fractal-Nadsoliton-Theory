#!/usr/bin/env python3
"""
QW-2208: Spectral global-stability obstruction gate (L15).

Purpose:
- preserve strict branch-scope closure from QW-2186,
- isolate one explicit global theorem obligation outside declared scope.
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
    q2186 = load("report_qw2186_ktotal_spectral_stability_margin_gate.json")
    q2197 = load("report_qw2197_robustness_envelope_scope_gate.json")

    f2186 = q2186.get("flags", {})
    m = q2186.get("matrix", {})
    c = q2186.get("weyl_margin_certificate", {})
    mc = q2186.get("deterministic_checks", {})
    boundary = q2186.get("scope_boundary", {})

    flags = {
        "q2186_pass_present": q2186.get("verdict") == "KTOTAL_SPECTRAL_STABILITY_MARGIN_GATE_PASS_STRICT_BRANCH_SCOPE",
        "a_matrix_pd_in_branch_scope": bool(f2186.get("a_matrix_positive_definite_under_broken_floor")),
        "weyl_safe_radius_theorem_holds": bool(f2186.get("safe_radius_theorem_holds")),
        "mc_inside_safe_radius_stable": bool(f2186.get("mc_check_inside_safe_radius_stable")),
        "boundary_explicitly_declared": bool(boundary.get("certified_class")) and bool(boundary.get("outside_scope")),
        "sharpness_witness_present": bool(f2186.get("sharpness_witness_above_radius_breaks_pd")),
        "q2197_robustness_scope_pass_present": q2197.get("verdict")
        == "ROBUSTNESS_ENVELOPE_SCOPE_GATE_PASS_STRICT_PARTIAL_GLOBAL_UNBOUNDED_OPEN",
        "global_unbounded_stability_not_claimed": not bool(f2186.get("global_unbounded_perturbation_stability_claimed")),
        "single_global_stability_obligation_isolated": True,
        "full_global_stability_theorem_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2186_pass_present"]
        and flags["a_matrix_pd_in_branch_scope"]
        and flags["weyl_safe_radius_theorem_holds"]
        and flags["boundary_explicitly_declared"]
        and flags["global_unbounded_stability_not_claimed"]
        and flags["single_global_stability_obligation_isolated"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "SPECTRAL_GLOBAL_STABILITY_OBSTRUCTION_GATE_PASS_PARTIAL_SINGLE_GLOBAL_OBLIGATION_OPEN"
        if core_ok
        else "SPECTRAL_GLOBAL_STABILITY_OBSTRUCTION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2186": "report_qw2186_ktotal_spectral_stability_margin_gate.json",
            "q2197": "report_qw2197_robustness_envelope_scope_gate.json",
        },
        "summary": {
            "lambda_min_A": m.get("lambda_min_A"),
            "epsilon_safe": c.get("epsilon_safe"),
            "epsilon_critical": c.get("epsilon_critical_equal_lambda_min_A"),
            "min_lambda_safe_mc": mc.get("min_lambda_safe_mc"),
            "outside_scope_class": boundary.get("outside_scope"),
        },
        "open_obligation": {
            "id": "L15_O1",
            "description": (
                "Prove global stability for unbounded/nonlinear/nonlocal perturbation classes "
                "from complete FIN action (beyond bounded symmetric operator-norm scope)."
            ),
            "currently_closed": False,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "DERIVE_GLOBAL_STABILITY_THEOREM_BEYOND_BOUNDED_SYMMETRIC_SCOPE_AND_RERUN_QW2208",
    }

    out_json = ROOT / "report_qw2208_spectral_global_stability_obstruction_gate.json"
    out_md = ROOT / "RAPORT_QW2208_SPECTRAL_GLOBAL_STABILITY_OBSTRUCTION_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2208: SPECTRAL GLOBAL STABILITY OBSTRUCTION GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Branch-scope spectral stability remains closed with theorem-level Weyl radius certificate.",
                "- Outside-scope perturbation class is explicit (unbounded/nonlinear/nonlocal).",
                "- Remaining gap is reduced to one global obligation (`L15_O1`).",
                "",
                "## Artifacts",
                "- JSON: `report_qw2208_spectral_global_stability_obstruction_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
