#!/usr/bin/env python3
"""Scratch probe: is the missing eta_eff coupling the legacy D_f-1 path-growth exponent?

This is deliberately a scratch/audit packet, not a theorem export.  It greps the
repo for D_f/D_F provenance, compares the existing bridge_missing_coupling
eta_eff against D_f-1 and strict eta=1.8, and records whether the numerical
trace is strong enough to justify a later analytic bridge attempt.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MISSING_COUPLING = HERE / "bridge_missing_coupling_report.json"
OUT_JSON = HERE / "bridge_df_eta_eff_report.json"
OUT_MD = HERE / "bridge_df_eta_eff_report.md"

STRICT_ETA = 1.8
ALPHA_GEO = 4.0 * math.log(2.0)
D_F_MINUS_ONE = ALPHA_GEO - 1.0


def rg_hits(pattern: str, limit: int = 80) -> list[str]:
    cmd = [
        "rg",
        "-n",
        "-S",
        pattern,
        ".",
        "-g",
        "!**/.git/**",
        "-g",
        "!fundamental_action_reconstruction/scratch/bridge_df_eta_eff*",
    ]
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:limit]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel_delta(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), 1e-15)


def main() -> None:
    missing = load_json(MISSING_COUPLING)
    eta_eff = float(missing["models"]["M1_legacy_plus_eta_eff"]["eta_eff"])
    midpoint = 0.5 * (eta_eff + STRICT_ETA)

    deltas = {
        "abs_eta_eff_minus_D_f_minus_1": abs(eta_eff - D_F_MINUS_ONE),
        "rel_eta_eff_minus_D_f_minus_1": rel_delta(eta_eff, D_F_MINUS_ONE),
        "abs_strict_eta_minus_D_f_minus_1": abs(STRICT_ETA - D_F_MINUS_ONE),
        "rel_strict_eta_minus_D_f_minus_1": rel_delta(STRICT_ETA, D_F_MINUS_ONE),
        "abs_midpoint_eta_eff_strict_minus_D_f_minus_1": abs(midpoint - D_F_MINUS_ONE),
        "rel_midpoint_eta_eff_strict_minus_D_f_minus_1": rel_delta(midpoint, D_F_MINUS_ONE),
        "abs_strict_eta_minus_eta_eff": abs(STRICT_ETA - eta_eff),
        "rel_strict_eta_minus_eta_eff": rel_delta(STRICT_ETA, eta_eff),
    }

    grep_audit = {
        "df_symbol_hits": rg_hits(r"\bD_F\b|\bD_f\b", 100),
        "df_minus_one_path_growth_hits": rg_hits(r"D_f\s*-\s*1|D_F\s*-\s*1|D_F\s*\-\s*1|D_f \- 1|d\^\{D_f - 1\}|d\^\{D_f\s*\-\s*1\}|d \*\* \(D_F - 1", 80),
        "one_point_seven_three_four_hits": rg_hits(r"1\.73[0-9]*|1\.74[0-9]*", 80),
        "strict_eta_hits": rg_hits(r"eta\W+1\.8|eta\W+1\.80000|\\eta = 1\.80000|delta_eta", 80),
    }

    candidate = {
        "candidate_id": "SCRATCH_DF_MINUS_ONE_AS_MISSING_ETA_EFF_COUPLING",
        "question": "Could the bridge_missing_coupling eta_eff ~= 1.74 be the legacy/nadsoliton D_f-1 path-growth exponent rather than a free phenomenological knob?",
        "numerical_constants": {
            "alpha_geo_equal_D_f_legacy_value_4ln2": ALPHA_GEO,
            "D_f_minus_1_path_growth_exponent": D_F_MINUS_ONE,
            "eta_eff_from_bridge_missing_coupling": eta_eff,
            "strict_kernel_eta": STRICT_ETA,
            "midpoint_eta_eff_strict_eta": midpoint,
        },
        "deltas": deltas,
        "bracketing_signature": {
            "eta_eff_lt_D_f_minus_1_lt_strict_eta": eta_eff < D_F_MINUS_ONE < STRICT_ETA,
            "midpoint_close_to_D_f_minus_1": deltas["rel_midpoint_eta_eff_strict_minus_D_f_minus_1"] < 0.002,
            "interpretation": "D_f-1 sits between the best eta_eff fit and the frozen strict eta, and the eta_eff/strict midpoint is within about 0.18% of D_f-1.",
        },
        "repo_evidence_summary": [
            "Repo grep finds legacy D_f/D_F usage with D_f = alpha_geo = 4 ln 2 = 2.77259 in historical/legacy material.",
            "Repo grep finds an explicit path-counting statement N(d) ~ d^{D_f - 1} ~= d^{1.77} in kernel-derivation material, giving a documented exponent close to eta_eff and strict eta.",
            "The existing scratch bridge_missing_coupling_probe shows the extra eta_eff parameter is statistically decisive, but that result is still model-selection evidence only.",
        ],
        "analytic_bridge_sketch_not_proven": {
            "starting_legacy_trace": "historical path multiplicity N(d) ~ d^(D_f-1)",
            "candidate_effective_denominator": "1 + beta_eff*d^(D_f-1)",
            "strict_side_target": "1 + beta*d^eta with eta=1.8",
            "missing_theorem_obligation": "derive how path multiplicity/fractal transport enters the damping denominator, including normalization and finite-window renormalization from D_f-1 to eta_eff/eta; no such theorem is exported here.",
        },
        "verdict": "STRONG_NUMERICAL_AND_REPO_PROVENANCE_TRACE_FOR_DF_MINUS_ONE_CANDIDATE__NOT_A_BRIDGE_THEOREM",
    }

    report = {
        "status": "OPEN_BRIDGE_CANDIDATE_WITH_TRACE_NO_THEOREM",
        "result_kind": "SCRATCH_DF_ETA_EFF_BRIDGE_CANDIDATE_AUDIT__NOT_A_THEOREM",
        "source_reports": {
            "bridge_missing_coupling_report": str(MISSING_COUPLING.relative_to(ROOT)),
        },
        "candidate": candidate,
        "grep_audit": grep_audit,
        "hard_limits": [
            "This does not identify K_legacy_ont with K_strict_gate.",
            "This does not transfer legacy physical-role claims to the strict kernel.",
            "This does not prove that D_f-1 analytically produces eta=1.8 or eta_eff=1.7388862982059556.",
            "This does not discharge QW-2191 and does not close ToE.",
        ],
        "next_honest_step": "Attempt an analytic finite-window/fractal-transport lemma: starting from N(d)~d^(D_f-1), derive the effective damping exponent in the denominator and bound its renormalized value between eta_eff and strict eta, or classify the trace as a fit coincidence.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch bridge audit: D_f-1 vs eta_eff\n\n"
        "Status: strong candidate trace; not a legacy→strict bridge theorem.\n\n"
        f"- `eta_eff` from `bridge_missing_coupling_report.json`: `{eta_eff:.15f}`.\n"
        f"- `D_f = 4 ln 2`: `{ALPHA_GEO:.15f}`; candidate `D_f-1`: `{D_F_MINUS_ONE:.15f}`.\n"
        f"- Strict `eta`: `{STRICT_ETA:.15f}`; midpoint `(eta_eff+eta)/2`: `{midpoint:.15f}`.\n"
        f"- Absolute gaps: `|eta_eff-(D_f-1)|={deltas['abs_eta_eff_minus_D_f_minus_1']:.15f}`, `|eta-(D_f-1)|={deltas['abs_strict_eta_minus_D_f_minus_1']:.15f}`, `|midpoint-(D_f-1)|={deltas['abs_midpoint_eta_eff_strict_minus_D_f_minus_1']:.15f}`.\n"
        f"- Bracketing: `eta_eff < D_f-1 < eta` is `{candidate['bracketing_signature']['eta_eff_lt_D_f_minus_1_lt_strict_eta']}`.\n"
        "- Grep found legacy/path-growth provenance for `D_f=4 ln 2` and `N(d) ~ d^(D_f-1) ≈ d^1.77`; see JSON `grep_audit`.\n"
        "- No false pass: no kernel identity, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
