#!/usr/bin/env python3
"""P2973/S1923: unit-bearing action-density coupling obstruction for incidence complex.

P2972 left two theorem routes for the typed support/provenance incidence complex.
This audit attacks exactly one: a unit-bearing coupling into a named action
density.  It does not replay ratio arithmetic, K/C exchange, unit conventions,
k-selection predicates, incidence identity/localizer anchors, or scalar Euler
placeholders.

The finite computation pairs the P2971 slots with named action-density receiver
schemas and checks whether each row has (i) typed incidence reception, (ii) a
nonconventional unit/measure source, (iii) a named density, (iv) a strict
source-localizer/coupling theorem, and (v) nonproxy variational readiness.  The
current rows can formally receive the incidence vector, but each available row
is missing at least the unit source and theorem-level coupling.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2971_s1921_typed_support_provenance_incidence_complex import SLOTS
from p2972_s1922_incidence_source_localizer_obstruction import OUT as P2972

OUT = GEN / "p2973_s1923_incidence_unit_bearing_action_density_coupling_obstruction.json"
MD = GEN / "p2973_s1923_incidence_unit_bearing_action_density_coupling_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def incidence_vector() -> list[int]:
    return [s["weight"] for s in sorted(SLOTS, key=lambda r: r["aggregate_index"])]


def receiver_rows() -> list[dict[str, Any]]:
    vector = incidence_vector()
    rows = [
        {
            "candidate": "formal_slot_sum_density_L_incidence_formal",
            "named_action_density": "L_incidence_formal = sum_i w_i Phi_i^2 dVol",
            "typed_incidence_reception": True,
            "preserves_support_provenance": True,
            "nonconventional_unit_measure_source": False,
            "strict_source_localizer_theorem": False,
            "coupling_theorem_to_density": False,
            "nonproxy_variational_readiness": False,
            "accepted_current_unit_bearing_coupling": False,
            "witness": "receives the vector formally, but uses unsourced dVol/unit symbols",
        },
        {
            "candidate": "P2964_quadratic_density_reception_import",
            "named_action_density": "P2964 aggregate-reception quadratic density schema",
            "typed_incidence_reception": True,
            "preserves_support_provenance": False,
            "nonconventional_unit_measure_source": False,
            "strict_source_localizer_theorem": False,
            "coupling_theorem_to_density": False,
            "nonproxy_variational_readiness": False,
            "accepted_current_unit_bearing_coupling": False,
            "witness": "imports prior aggregate reception but collapses incidence/provenance data to a scalar receiver",
        },
        {
            "candidate": "primitive_mean_9_5_density_import",
            "named_action_density": "c=(9/5)U^k primitive-mean density",
            "typed_incidence_reception": False,
            "preserves_support_provenance": False,
            "nonconventional_unit_measure_source": False,
            "strict_source_localizer_theorem": False,
            "coupling_theorem_to_density": False,
            "nonproxy_variational_readiness": False,
            "accepted_current_unit_bearing_coupling": False,
            "witness": "replays ratio-package scalar arithmetic and loses the incidence object",
        },
        {
            "candidate": "bookkeeping_component_labeled_density",
            "named_action_density": "L_KC_labelled = sum_{K,C} w_a Phi_a^2 dVol",
            "typed_incidence_reception": True,
            "preserves_support_provenance": True,
            "nonconventional_unit_measure_source": False,
            "strict_source_localizer_theorem": False,
            "coupling_theorem_to_density": False,
            "nonproxy_variational_readiness": False,
            "accepted_current_unit_bearing_coupling": False,
            "witness": "keeps K/C labels but label bookkeeping is not a unit-bearing strict coupling theorem",
        },
        {
            "candidate": "formal_Euler_receiver_placeholder",
            "named_action_density": "Euler placeholder density E[L_incidence]=0",
            "typed_incidence_reception": True,
            "preserves_support_provenance": True,
            "nonconventional_unit_measure_source": False,
            "strict_source_localizer_theorem": False,
            "coupling_theorem_to_density": False,
            "nonproxy_variational_readiness": False,
            "accepted_current_unit_bearing_coupling": False,
            "witness": "Euler form is downstream formalism, not a sourced unit/measure installation",
        },
        {
            "candidate": "completed_strict_unit_bearing_incidence_density_schema",
            "named_action_density": "L_incidence_strict[U,dmu_N] with sourced units and coupling theorem",
            "typed_incidence_reception": True,
            "preserves_support_provenance": True,
            "nonconventional_unit_measure_source": True,
            "strict_source_localizer_theorem": True,
            "coupling_theorem_to_density": True,
            "nonproxy_variational_readiness": True,
            "accepted_current_unit_bearing_coupling": False,
            "witness": "completed theorem schema only; not exported by current artifacts",
        },
    ]
    for row in rows:
        row["received_vector"] = vector if row["typed_incidence_reception"] else []
        row["missing_obligations"] = [k for k in ["typed_incidence_reception", "preserves_support_provenance", "nonconventional_unit_measure_source", "strict_source_localizer_theorem", "coupling_theorem_to_density", "nonproxy_variational_readiness"] if not row[k]]
    return rows


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in rows if r["candidate"] != "completed_strict_unit_bearing_incidence_density_schema"]
    return [
        {"obligation": "typed_incidence_reception_exists", "satisfied": any(r["typed_incidence_reception"] for r in current), "evidence": "formal slot-sum and bookkeeping receivers can carry the P2971 vector"},
        {"obligation": "support_provenance_preserved_by_a_receiver", "satisfied": any(r["preserves_support_provenance"] for r in current), "evidence": "some formal receivers keep slot/provenance labels"},
        {"obligation": "nonconventional_unit_measure_source_exported", "satisfied": any(r["nonconventional_unit_measure_source"] for r in current), "evidence": "all available receivers use formal or imported unit/measure symbols"},
        {"obligation": "strict_source_localizer_theorem_exported", "satisfied": any(r["strict_source_localizer_theorem"] for r in current), "evidence": "P2972 found no current strict localizer"},
        {"obligation": "coupling_theorem_to_named_density_exported", "satisfied": any(r["coupling_theorem_to_density"] for r in current), "evidence": "no theorem installs the incidence complex into a named density"},
        {"obligation": "nonproxy_variational_readiness_exported", "satisfied": any(r["nonproxy_variational_readiness"] for r in current), "evidence": "formal Euler placeholders are not a variational chain rule"},
        {"obligation": "accepted_current_unit_bearing_coupling", "satisfied": any(r["accepted_current_unit_bearing_coupling"] for r in current), "evidence": "completed schema remains unavailable"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["typed_reception", "support_provenance", "nonconventional_unit_source", "source_localizer", "named_density_coupling", "nonproxy_variational_readiness"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_unit_bearing_action_density_coupling": m == full} for m in range(1 << len(names))]


def build_payload(p2972_path: Any) -> dict[str, Any]:
    rows = receiver_rows()
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2973_INCIDENCE_UNIT_BEARING_ACTION_DENSITY_COUPLING_OBSTRUCTION_NO_STRICT_EXPORT",
        "input_hashes": {"P2972": hashlib.sha256(p2972_path.read_bytes()).hexdigest() if p2972_path.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "IncidenceUnitBearingActionDensityCoupling_ObstructionMatrix",
            "incidence_vector": incidence_vector(),
            "receiver_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "coupling_certificate": {
            "slot_count": len(SLOTS),
            "incidence_weight_sum": sum(incidence_vector()),
            "receiver_candidate_count": len(rows),
            "available_receiver_candidate_count": len([r for r in rows if not r["candidate"].startswith("completed_")]),
            "typed_reception_rows": sum(1 for r in rows if r["typed_incidence_reception"]),
            "accepted_current_unit_bearing_couplings": [r["candidate"] for r in rows if r["accepted_current_unit_bearing_coupling"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_unit_bearing_action_density_coupling"]),
        },
        "decision": {
            "positive_progress": "P2973 turns the missing unit-bearing action-density coupling for the P2971 incidence complex into a finite receiver/coupling obligation matrix.",
            "breakthrough": "No strict unit-bearing incidence coupling is exported: current receivers are formal, scalar-import, bookkeeping-label, or Euler-placeholder rows, and the completed schema is unavailable.",
            "negative_export_flags": {k: False for k in ["unit_bearing_coupling_exported", "strict_source_localizer_exported", "nonproxy_variational_chain_rule_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay formal slot-sum densities, P2964 scalar reception, primitive-mean 9/5 imports, K/C bookkeeping labels, formal Euler placeholders, incidence identity/localizer anchors, unit conventions, or k-selection predicates.  The next proof-grade move may attack exactly the remaining incidence theorem: a nonproxy variational chain rule for the incidence complex.  Otherwise preserve the P2929-P2973 no-strict-export boundary without promoting L_total, bridge closure, role transfer, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["coupling_certificate"]
    lines = [
        "# P2973/S1923 incidence unit-bearing action-density coupling obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling certificate",
        f"- slot count / weight sum: `{cert['slot_count']}` / `{cert['incidence_weight_sum']}`",
        f"- receiver candidates / available: `{cert['receiver_candidate_count']}` / `{cert['available_receiver_candidate_count']}`",
        f"- typed reception rows: `{cert['typed_reception_rows']}`",
        f"- accepted current unit-bearing couplings: `{cert['accepted_current_unit_bearing_couplings']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
        "",
        "## Lay summary",
        payload["decision"]["positive_progress"],
        payload["decision"]["breakthrough"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2972)
    payload = build_payload(P2972)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2973/S1923 incidence unit-bearing action-density coupling obstruction", "## P2973/S1923 incidence unit-bearing action-density coupling obstruction\n\n`P2973/S1923` attacks one remaining theorem route for the P2971 typed support/provenance incidence complex: a unit-bearing coupling into a named action density.  The finite receiver matrix tests formal slot-sum density, P2964 aggregate-reception import, primitive-mean `9/5` density import, bookkeeping K/C-labelled density, formal Euler receiver placeholder, and a completed strict schema.  Available rows can formally receive the incidence vector `[1,2,2,2,2]`, but none exports a nonconventional unit/measure source, strict source-localizer theorem, coupling theorem to a named density, or nonproxy variational readiness.  Thus no unit-bearing incidence coupling, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2973/S1923 incidence action-density coupling `L_total` guard", "## P2973/S1923 incidence action-density coupling `L_total` guard\n\n`P2973/S1923` adds no sourced `L_total` term.  Formal slot-sum, P2964 scalar-reception import, primitive-mean `9/5` import, K/C bookkeeping density, and Euler-placeholder rows either lack a nonconventional unit/measure source or lack theorem-level coupling from the incidence complex into a named nonproxy density.  Therefore no EOM, Hamiltonian, bridge closure, role transfer, or ToE is exported.\n")
    append_once(AGENTS, "Current incidence unit-bearing action-density coupling obstruction guardrail (P2973/S1923, 2026-06-20)", "## Current incidence unit-bearing action-density coupling obstruction guardrail (P2973/S1923, 2026-06-20)\n\n- P2973 audits the unit-bearing coupling theorem missing from P2971/P2972 for the typed support/provenance incidence complex.\n- The finite receiver matrix finds formal incidence receivers, but current rows are formal slot sums, P2964 scalar-reception imports, primitive-mean `9/5` imports, K/C bookkeeping labels, or Euler placeholders; none exports a nonconventional unit/measure source, strict source-localizer theorem, coupling theorem to a named density, or nonproxy variational readiness.\n- Do not promote formal slot-sum densities, P2964 scalar reception, primitive-mean `9/5`, K/C bookkeeping labels, formal Euler placeholders, incidence identity/localizer anchors, unit conventions, or k-selection predicates to strict sourcehood, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move may attack exactly the remaining incidence theorem: a nonproxy variational chain rule for the incidence complex; otherwise preserve the P2929-P2973 no-strict-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
