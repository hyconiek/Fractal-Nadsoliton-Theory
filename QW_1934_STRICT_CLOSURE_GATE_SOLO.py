#!/usr/bin/env python3
"""
QW-1934: Strict solo closure gate after eta reparameterization branch.

Purpose:
- integrate key gates after QW-1932/1933,
- determine whether the "maximum closure without external team" is achieved,
- separate closure status from later SM+GR head-to-head prediction stage.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1934_strict_closure_gate_solo.json"
OUT_MD = ROOT / "RAPORT_QW1934_STRICT_CLOSURE_GATE_SOLO.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1927 = load("report_qw1927_true_external_beta_channel_readiness_gate.json")
    d1922 = load("report_qw1922_beta_observable_blind_external_intervention.json")
    d1918 = load("report_qw1918_triad_blind_external_validation.json")
    d1931 = load("report_qw1931_strict_triad_feasibility_frontier.json")
    d1932 = load("report_qw1932_physical_reparameterization_eta_scan.json")
    d1933 = load("report_qw1933_reparam_hash_split_robustness.json")

    external_ready = d1927.get("readiness") == "READY"
    intervention_pass = d1922.get("verdict") == "BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS"
    blind_external_pass = str(d1918.get("verdict", "")).startswith("TRIAD_BLIND_EXTERNAL_VALIDATION_PASS")

    sel1932 = d1932.get("selected", {})
    comp1932 = sel1932.get("comparison", {})
    strict1932_flags = comp1932.get("strict_flags", {})
    reparam_strict_pass = (
        d1932.get("verdict") == "PHYSICAL_REPARAMETERIZATION_STRICT_PASS"
        and int(d1932.get("strict_pass_count", 0)) >= 1
        and bool(comp1932.get("strict_all_pass", False))
        and bool(all(bool(v) for v in strict1932_flags.values()))
    )

    flags1933 = d1933.get("flags", {})
    robust1933_pass = d1933.get("verdict") == "REPARAM_HASH_SPLIT_ROBUSTNESS_PASS" and bool(
        all(bool(v) for v in flags1933.values())
    )

    gap1931_resolved = (
        d1931.get("verdict") == "STRICT_TRIAD_FEASIBILITY_FAIL"
        and d1932.get("verdict") == "PHYSICAL_REPARAMETERIZATION_STRICT_PASS"
    )

    primary_gain_median = float(d1933["datasets"]["primary"]["summary"]["gain_median"])
    primary_corr_median = float(d1933["datasets"]["primary"]["summary"]["corr_median"])
    stress_gain_median = float(d1933["datasets"]["stress"]["summary"]["gain_median"])
    stress_corr_median = float(d1933["datasets"]["stress"]["summary"]["corr_median"])

    signal_floor_flags = {
        "primary_corr_median_ge_0p05": bool(primary_corr_median >= 0.05),
        "primary_gain_median_ge_0p001": bool(primary_gain_median >= 0.001),
        "stress_corr_median_ge_0p25": bool(stress_corr_median >= 0.25),
        "stress_gain_median_ge_0p03": bool(stress_gain_median >= 0.03),
    }

    solo_hard_flags = {
        "external_ready": bool(external_ready),
        "intervention_pass": bool(intervention_pass),
        "blind_external_pass": bool(blind_external_pass),
        "gap1931_resolved_by_reparam": bool(gap1931_resolved),
        "reparam_strict_pass": bool(reparam_strict_pass),
        "reparam_multisplit_robust_pass": bool(robust1933_pass),
        **signal_floor_flags,
    }

    if all(solo_hard_flags.values()):
        readiness = "TOE_STAGE_B_SOLO_CLOSED"
        verdict = "STRICT_CLOSURE_GATE_SOLO_PASS"
        required_next = "FREEZE_PREDICTIONS_AND_RUN_SM_GR_HEAD_TO_HEAD_PROTOCOL"
    else:
        readiness = "TOE_STAGE_B_SOLO_OPEN"
        verdict = "STRICT_CLOSURE_GATE_SOLO_PARTIAL_OR_FAIL"
        required_next = "FIX_FAILED_FLAGS_BEFORE_SM_GR_COMPARISON"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1927_verdict": d1927.get("verdict"),
            "q1922_verdict": d1922.get("verdict"),
            "q1918_verdict": d1918.get("verdict"),
            "q1931_verdict": d1931.get("verdict"),
            "q1932_verdict": d1932.get("verdict"),
            "q1933_verdict": d1933.get("verdict"),
        },
        "diagnostics": {
            "q1932_selected": {
                "eta": sel1932.get("eta"),
                "omega": sel1932.get("fit", {}).get("omega"),
                "phi": sel1932.get("fit", {}).get("phi"),
                "beta": sel1932.get("fit", {}).get("beta"),
                "rel_loss_vs_eta1": comp1932.get("rel_loss_vs_eta1"),
            },
            "q1933_primary_summary": d1933["datasets"]["primary"]["summary"],
            "q1933_stress_summary": d1933["datasets"]["stress"]["summary"],
        },
        "flags": solo_hard_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1934: STRICT CLOSURE GATE SOLO",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Key Inputs",
        f"- QW-1927: {d1927.get('verdict')}",
        f"- QW-1922: {d1922.get('verdict')}",
        f"- QW-1918: {d1918.get('verdict')}",
        f"- QW-1931: {d1931.get('verdict')}",
        f"- QW-1932: {d1932.get('verdict')}",
        f"- QW-1933: {d1933.get('verdict')}",
        "",
        "## Selected Reparam Candidate (QW-1932)",
        (
            "- eta/omega/phi/beta: "
            f"{sel1932.get('eta')}/{sel1932.get('fit', {}).get('omega')}/"
            f"{sel1932.get('fit', {}).get('phi')}/{sel1932.get('fit', {}).get('beta')}"
        ),
        f"- rel_loss_vs_eta1: {comp1932.get('rel_loss_vs_eta1')}",
        "",
        "## QW-1933 Median Signal",
        f"- primary corr/gain: {primary_corr_median:.4f}/{primary_gain_median:.4f}",
        f"- stress corr/gain: {stress_corr_median:.4f}/{stress_gain_median:.4f}",
        "",
        "## Strict Flags",
    ]
    for k, v in solo_hard_flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1934] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1934] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1934] readiness={readiness} verdict={verdict}")


if __name__ == "__main__":
    main()

