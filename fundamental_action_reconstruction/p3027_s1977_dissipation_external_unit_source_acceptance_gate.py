#!/usr/bin/env python3
"""P3027/S1977: external unit/source theorem acceptance gate for P3026.

P3026 says the dissipation time-order lane can only reopen if a genuinely new
strict typed object or external physical unit/source theorem is supplied.  This
module builds that missing acceptance object and runs it on the currently
available candidate family, without promoting placeholders or internal ratios.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3026_s1976_dissipation_time_order_reconciliation_no_new_live_frontier import OUT as P3026

OUT = GEN / "p3027_s1977_dissipation_external_unit_source_acceptance_gate.json"
MD = GEN / "p3027_s1977_dissipation_external_unit_source_acceptance_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

OBLIGATIONS = [
    "new_non_replay_typed_object",
    "strict_nadsoliton_provenance",
    "nonpremise_chart_or_orientation_source",
    "absolute_physical_tick_or_action_unit",
    "scale_orbit_breaking",
    "explicit_hamiltonian_coupling",
]


def candidate_rows() -> list[dict[str, Any]]:
    raw = [
        {
            "candidate": "P3023_monotone_dissipation_chain",
            "new_non_replay_typed_object": False,
            "strict_nadsoliton_provenance": True,
            "nonpremise_chart_or_orientation_source": False,
            "absolute_physical_tick_or_action_unit": False,
            "scale_orbit_breaking": False,
            "explicit_hamiltonian_coupling": False,
            "blocker": "already reconciled; directed scaffold is chart-dependent and unitless",
        },
        {
            "candidate": "P3024_chart_orbit_anchors",
            "new_non_replay_typed_object": False,
            "strict_nadsoliton_provenance": True,
            "nonpremise_chart_or_orientation_source": False,
            "absolute_physical_tick_or_action_unit": False,
            "scale_orbit_breaking": False,
            "explicit_hamiltonian_coupling": False,
            "blocker": "four chart representatives and trivial stabilizer; anchors inherit labels",
        },
        {
            "candidate": "P3025_internal_tick_H_equals_S_over_tau_ratios",
            "new_non_replay_typed_object": False,
            "strict_nadsoliton_provenance": True,
            "nonpremise_chart_or_orientation_source": False,
            "absolute_physical_tick_or_action_unit": False,
            "scale_orbit_breaking": False,
            "explicit_hamiltonian_coupling": False,
            "blocker": "positive internal calibrations but no external physical unit or energy theorem",
        },
        {
            "candidate": "formal_imported_unit_symbol_U_ext",
            "new_non_replay_typed_object": True,
            "strict_nadsoliton_provenance": False,
            "nonpremise_chart_or_orientation_source": False,
            "absolute_physical_tick_or_action_unit": True,
            "scale_orbit_breaking": True,
            "explicit_hamiltonian_coupling": False,
            "blocker": "an imported symbol can set units by fiat but lacks strict provenance and coupling theorem",
        },
    ]
    rows = []
    for row in raw:
        passed = [name for name in OBLIGATIONS if row[name]]
        failed = [name for name in OBLIGATIONS if not row[name]]
        rows.append({**row, "passed_obligations": passed, "failed_obligations": failed, "accepted": not failed})
    return rows


def build_payload(p3026_path: Any) -> dict[str, Any]:
    read_json(p3026_path)
    rows = candidate_rows()
    return {
        "status": "P3027_DISSIPATION_EXTERNAL_UNIT_SOURCE_ACCEPTANCE_GATE_NO_ACCEPTED_SOURCE",
        "input_hashes": {"P3026": hashlib.sha256(p3026_path.read_bytes()).hexdigest() if p3026_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "DissipationExternalUnitSource_AcceptanceGate",
            "obligations": OBLIGATIONS,
            "candidate_rows": rows,
        },
        "finite_certificate": {
            "candidate_count": len(rows),
            "obligation_count": len(OBLIGATIONS),
            "accepted_candidate_count": sum(1 for row in rows if row["accepted"]),
            "best_pass_count": max(len(row["passed_obligations"]) for row in rows),
            "imported_symbol_accepted": next(row["accepted"] for row in rows if row["candidate"] == "formal_imported_unit_symbol_U_ext"),
        },
        "decision": {
            "breakthrough": "A concrete acceptance gate for reopening the P3026 time-order lane was constructed.  Current internal candidates fail as replays, and the only formal external unit symbol breaks scale by fiat but fails strict nadsoliton provenance, orientation/chart source, and Hamiltonian coupling.  Therefore no accepted external physical unit/source theorem is supplied.",
            "negative_export_flags": {k: False for k in ["external_unit_source_exported", "strict_time_order_object_exported", "physical_tick_theorem_exported", "hamiltonian_exported", "unit_bearing_action_eom_source_exported", "time_arrow_exported", "observed_physics_exported", "ltotal_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Use the P3027 acceptance gate only when a concrete new formula/source is supplied.  Without that, preserve the P3026/P3027 no-accepted-source boundary and pivot to a different typed object outside the dissipation time-order and internal-unit lanes.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3027/S1977 dissipation external unit/source acceptance gate", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- candidates: `{c['candidate_count']}`",
        f"- obligations: `{c['obligation_count']}`",
        f"- accepted candidates: `{c['accepted_candidate_count']}`",
        f"- best pass count: `{c['best_pass_count']}`",
        f"- imported symbol accepted: `{c['imported_symbol_accepted']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(P3026)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3027/S1977 dissipation external unit/source acceptance gate", "## P3027/S1977 dissipation external unit/source acceptance gate\n\n`P3027/S1977` constructs the acceptance gate required to reopen the P3026 dissipation time-order lane: a candidate must be a new non-replay typed object with strict nadsoliton provenance, nonpremise chart/orientation source, absolute physical tick or action unit, scale-orbit breaking, and explicit Hamiltonian coupling.  The four tested current candidates have zero accepted rows: P3023/P3024/P3025 candidates are replays/internal, while the formal imported unit symbol breaks scale only by fiat and lacks strict provenance and coupling.  No external unit source, time arrow, Hamiltonian, unit-bearing EOM, observed-physics export, `L_total`, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3027/S1977 external unit/source acceptance `L_total` guard", "## P3027/S1977 external unit/source acceptance `L_total` guard\n\n`P3027/S1977` adds no physical `L_total` term.  It exports an acceptance gate for future external unit/source candidates, but the current candidate family has zero accepted rows; imported unit symbols are not strict nadsoliton-sourced Hamiltonian couplings.\n")
    append_once(AGENTS, "Current dissipation external unit/source acceptance guardrail (P3027/S1977, 2026-06-22)", "## Current dissipation external unit/source acceptance guardrail (P3027/S1977, 2026-06-22)\n\n- P3027 constructs the acceptance gate required to reopen the P3026 dissipation time-order lane with an external physical unit/source theorem.\n- A candidate must pass six obligations: new non-replay typed object, strict nadsoliton provenance, nonpremise chart/orientation source, absolute physical tick or action unit, scale-orbit breaking, and explicit Hamiltonian coupling.\n- Current candidates have zero accepted rows: P3023/P3024/P3025 are replays/internal, and an imported unit symbol is not strict provenance or Hamiltonian coupling.\n- Do not promote external unit placeholders or internal candidates to physical time arrow, unit-bearing EOM/Hamiltonian, observed-physics, `L_total`, selector, bridge/role-transfer, or ToE closure; use the gate only when a concrete new formula/source is supplied.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
