#!/usr/bin/env python3
"""
QW-2082: Missing-14 strict closure roadmap.

Purpose:
- translate QW-2081 frontier diagnosis into an actionable, parameter-by-parameter
  closure roadmap under strict scientific rigor (no retune, no scan),
- define required new observables and hard closure conditions for each unresolved item.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
IN_2081 = ROOT / "report_qw2081_missing14_strict_rigor_frontier.json"
OUT_JSON = ROOT / "report_qw2082_missing14_strict_closure_roadmap.json"
OUT_MD = ROOT / "RAPORT_QW2082_MISSING14_STRICT_CLOSURE_ROADMAP.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def row(
    pid: str,
    tier: str,
    required_new_observables: List[str],
    hard_closure_condition: str,
    next_gate: str,
    notes: str,
) -> Dict:
    return {
        "id": pid,
        "priority_tier": tier,
        "required_new_observables": required_new_observables,
        "hard_closure_condition": hard_closure_condition,
        "next_gate": next_gate,
        "notes": notes,
    }


def main() -> None:
    r2081 = load_json(IN_2081)
    unresolved = sorted(set(r2081["strict_unresolved_ids"]))

    roadmap_catalog = {
        "delta_cp_ckm": row(
            "delta_cp_ckm",
            "T1",
            ["independent_ckm_cp_sensitive_observables", "deterministic_complex_phase_refinement"],
            "CKM CP phase branch selected deterministically and matched within registry tolerance without scan/retune.",
            "QW_2097_CKM_CP_TARGET_REFINEMENT_GATE",
            "Now a strict target-miss item; requires precision refinement, not a loose fit.",
        ),
        "alpha_s_mz": row(
            "alpha_s_mz",
            "T1",
            ["independent_high_precision_qcd_event_shapes", "locked_nonanchor_rge_boundary"],
            "alpha_s(MZ) derived from non-anchored boundary within registry tolerance in blind holdout.",
            "QW_2087_ALPHA_S_NONANCHOR_BOUNDARY_GATE",
            "Fastest high-impact closure candidate in strong-interaction branch.",
        ),
        "g_f": row(
            "g_f",
            "T1",
            ["independent_muon_lifetime_chain", "v_higgs_nonanchor_map"],
            "G_F derived from independently closed v-map (not identity over anchored v).",
            "QW_2085_GF_NONANCHOR_LIFETIME_GATE",
            "Must break anchor-loop with independent lifetime-based closure map.",
        ),
        "m_z": row(
            "m_z",
            "T1",
            ["independent_electroweak_pole_inputs", "delta_r_full_nonanchor_chain"],
            "M_Z obtained from non-anchored electroweak closure and matches registry tolerance.",
            "QW_2086_MZ_NONANCHOR_EW_POLE_GATE",
            "Requires full non-anchored electroweak pole reconstruction.",
        ),
        "m_up": row(
            "m_up",
            "T2",
            ["nonperturbative_qcd_sumrule_or_lattice_bridge", "yukawa_running_nonanchor_inputs"],
            "m_u closed without direct anchor injection and validated on independent hadronic observables.",
            "QW_2088_LIGHT_QUARK_MASS_NONANCHOR_GATE",
            "Light-quark masses are baseline-derived today; need external non-anchor bridge.",
        ),
        "m_down": row(
            "m_down",
            "T2",
            ["nonperturbative_qcd_sumrule_or_lattice_bridge", "yukawa_running_nonanchor_inputs"],
            "m_d closed without direct anchor injection and validated on independent hadronic observables.",
            "QW_2088_LIGHT_QUARK_MASS_NONANCHOR_GATE",
            "Shared with m_u/m_s closure chain.",
        ),
        "m_strange": row(
            "m_strange",
            "T2",
            ["nonperturbative_qcd_sumrule_or_lattice_bridge", "yukawa_running_nonanchor_inputs"],
            "m_s closed without direct anchor injection and validated on independent hadronic observables.",
            "QW_2088_LIGHT_QUARK_MASS_NONANCHOR_GATE",
            "Shared with m_u/m_d closure chain.",
        ),
        "m_h": row(
            "m_h",
            "T2",
            ["higgs_self_coupling_independent_input", "nonanchor_radiative_higgs_chain"],
            "Higgs pole mass derived from closed lambda-running chain without model-assumption plug-ins.",
            "QW_2089_HIGGS_SELFCOUPLING_STRICT_GATE",
            "Current proxy depends on fixed loop assumptions; needs strict closure inputs.",
        ),
        "h0": row(
            "h0",
            "T3",
            ["late_time_distance_ladder_independent", "early_time_cmb_independent"],
            "H0 derived from internally fixed cosmological sector without external anchor in lambda-H0 pair.",
            "QW_2090_H0_LAMBDA_DECOUPLING_GATE",
            "Currently coupled-underdetermined with Lambda; must be decoupled.",
        ),
        "lambda_cosmological": row(
            "lambda_cosmological",
            "T3",
            ["vacuum_response_observable", "independent_background_expansion_chain"],
            "Lambda derived from internal vacuum response observable and cross-validated against expansion history.",
            "QW_2090_H0_LAMBDA_DECOUPLING_GATE",
            "Cannot remain fixed from external anchor in strict final claim.",
        ),
        "m_nu1": row(
            "m_nu1",
            "T3",
            ["absolute_neutrino_mass_scale_observable", "beta_decay_or_cosmo_bound_integration"],
            "Absolute neutrino mass branch closed without hierarchy-assumption fixes.",
            "QW_2091_NEUTRINO_ABSOLUTE_SCALE_GATE",
            "Current status depends on hierarchy assumptions.",
        ),
        "m_nu2": row(
            "m_nu2",
            "T3",
            ["absolute_neutrino_mass_scale_observable", "beta_decay_or_cosmo_bound_integration"],
            "Absolute neutrino mass branch closed without hierarchy-assumption fixes.",
            "QW_2091_NEUTRINO_ABSOLUTE_SCALE_GATE",
            "Current status depends on anchored oscillation deltas.",
        ),
        "m_nu3": row(
            "m_nu3",
            "T3",
            ["absolute_neutrino_mass_scale_observable", "beta_decay_or_cosmo_bound_integration"],
            "Absolute neutrino mass branch closed without hierarchy-assumption fixes.",
            "QW_2091_NEUTRINO_ABSOLUTE_SCALE_GATE",
            "Current status depends on anchored oscillation deltas.",
        ),
        "G_newton": row(
            "G_newton",
            "T4",
            ["dimensionless_to_si_gravity_bridge_observable", "independent_planck_to_ir_scale_map"],
            "G derived from internal dimensionless chain with explicit SI bridge and blind consistency checks.",
            "QW_2092_GNEWTON_SI_BRIDGE_GATE",
            "Hardest unresolved: currently underdetermined without a new bridge observable.",
        ),
    }

    unknown = sorted([pid for pid in unresolved if pid not in roadmap_catalog])
    if unknown:
        raise RuntimeError(f"QW-2082 missing roadmap definitions for IDs: {unknown}")

    roadmap = [roadmap_catalog[pid] for pid in unresolved]

    by_tier: Dict[str, int] = {}
    for r in roadmap:
        by_tier[r["priority_tier"]] = by_tier.get(r["priority_tier"], 0) + 1

    if unresolved:
        verdict = "MISSING14_STRICT_CLOSURE_ROADMAP_READY"
        required_next_step = "EXECUTE_TOP_PRIORITY_OPEN_TIERS_FIRST_THEN_PROGRESS_TO_T3_T4"
    else:
        verdict = "MISSING14_STRICT_CLOSURE_ROADMAP_ALL_CURRENTLY_CLOSED"
        required_next_step = "EXECUTE_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_frontier": IN_2081.name,
        "input_verdict": r2081["verdict"],
        "strict_unresolved_ids": unresolved,
        "roadmap": roadmap,
        "tier_counts": by_tier,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2082: MISSING-14 STRICT CLOSURE ROADMAP",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Input frontier: `{IN_2081.name}`",
        f"- Input verdict: `{out['input_verdict']}`",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Tier Counts",
    ]
    for k in sorted(by_tier.keys()):
        lines.append(f"- {k}: {by_tier[k]}")

    lines.extend(["", "## Parameter Plan"])
    for r in roadmap:
        lines.append(
            f"- {r['id']}: tier={r['priority_tier']} | next_gate={r['next_gate']}"
        )

    lines.extend(["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2082] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2082] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2082] verdict={out['verdict']}")


if __name__ == "__main__":
    main()
