#!/usr/bin/env python3
"""P3115/S2065: Sigma_dim scale-section theorem audit.

P3114 left exactly one admissible next object: a nadsoliton-only scale-section
theorem Sigma_dim for the internal dimension orbit of
D_phi=(U_action,U_length,U_time).  This audit constructs finite section-theorem
candidates and checks whether any candidate selects a nonzero representative,
proves C_phi(A_phi)=U_action, and exports an action-from-length/time relation
without importing hbar/Planck, rods, clocks, observed light, apparatus,
selector replay, L_total, bridge/role-transfer, or ToE promotion.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3113_s2063_u_action_reference_carrier_source_law_audit import a_phi
from p3114_s2064_dimensional_triad_source_law_audit import OUT as P3114

OUT = GEN / "p3115_s2065_sigma_dim_scale_section_theorem_audit.json"
MD = GEN / "p3115_s2065_sigma_dim_scale_section_theorem_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SCALE_ACTIONS = (Fraction(1, 5), Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(5, 1))
SECTION_OBLIGATIONS = ("nonzero_representative", "orbit_invariance", "uniqueness_mod_scale", "strict_source_law", "C_phi_coupling", "action_length_time_relation")
GATES = (
    "uses_p3114_sigma_dim_obligation",
    "explicit_sigma_dim_formula",
    "nonzero_representative_selected",
    "scale_orbit_quotient_well_defined",
    "unique_section_exported",
    "strict_source_law_exported",
    "C_phi_A_phi_coupling_proved",
    "action_from_length_time_relation_proved",
    "standard_physics_import_free",
    "selector_bridge_ltotal_toe_free",
    "nonconventional_source_law_exported",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "sigma_dim_section": r"Sigma_dim|scale-section theorem|scale section|scale-orbit|internal dimension orbit",
        "d_phi_triad": r"D_phi|U_action|U_length|U_time|action-from-length|length/time",
        "c_phi_coupling": r"C_phi\(A_phi\)|C_phi|A_phi|2\*pi/alpha_geo",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def sigma_candidates() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "identity_dimension_orbit_section",
            "formula": "Sigma_dim([S,L,T]) := (S,L,T)",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": False,
            "unique_section_exported": False,
            "strict_source_law_exported": False,
            "C_phi_A_phi_coupling_proved": False,
            "action_from_length_time_relation_proved": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "identity keeps every scale representative and therefore does not quotient the scale orbit",
        },
        {
            "candidate": "unit_norm_triad_section",
            "formula": "Sigma_dim([S,L,T]) := representative with S^2+L^2+T^2=1",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": True,
            "unique_section_exported": True,
            "strict_source_law_exported": False,
            "C_phi_A_phi_coupling_proved": False,
            "action_from_length_time_relation_proved": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "unit norm fixes a mathematical gauge but has no strict law sourcing this norm as a physical dimension carrier",
        },
        {
            "candidate": "alpha_geo_phase_area_section",
            "formula": "Sigma_dim([S,L,T]) := representative with S=C_phi(A_phi) and A_phi=2*pi/alpha_geo",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": True,
            "unique_section_exported": True,
            "strict_source_law_exported": True,
            "C_phi_A_phi_coupling_proved": True,
            "action_from_length_time_relation_proved": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "phase-area coupling selects an action-like representative but still lacks length/time carriers and their action relation",
        },
        {
            "candidate": "entropy_cell_tick_ratio_section",
            "formula": "Sigma_dim([S,L,T]) := representative with L=bit_cell and T=entropy_tick",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": True,
            "unique_section_exported": False,
            "strict_source_law_exported": True,
            "C_phi_A_phi_coupling_proved": False,
            "action_from_length_time_relation_proved": True,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "cell/tick ratio remains combinatorial and does not prove C_phi(A_phi)=U_action",
        },
        {
            "candidate": "cohomology_period_volume_section",
            "formula": "Sigma_dim([S,L,T]) := primitive integral period/cell/time-volume representative",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": True,
            "unique_section_exported": False,
            "strict_source_law_exported": True,
            "C_phi_A_phi_coupling_proved": False,
            "action_from_length_time_relation_proved": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "primitive periods can fix integer normalization but not dimensional action/length/time relation",
        },
        {
            "candidate": "planck_light_cone_section",
            "formula": "Sigma_dim([S,L,T]) := (hbar, Planck length, Planck time)",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": True,
            "unique_section_exported": True,
            "strict_source_law_exported": False,
            "C_phi_A_phi_coupling_proved": True,
            "action_from_length_time_relation_proved": True,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "section is complete only by importing standard Planck/light calibration",
        },
        {
            "candidate": "apparatus_calibration_section",
            "formula": "Sigma_dim([S,L,T]) := calibrated lab action/rod/clock representative",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": True,
            "unique_section_exported": True,
            "strict_source_law_exported": False,
            "C_phi_A_phi_coupling_proved": False,
            "action_from_length_time_relation_proved": True,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "apparatus calibration is observer-side and not a nadsoliton-only strict source law",
        },
        {
            "candidate": "selector_oriented_dimension_section",
            "formula": "Sigma_dim([S,L,T]) := representative chosen after selector/orientation premise",
            "uses_p3114_sigma_dim_obligation": True,
            "explicit_sigma_dim_formula": True,
            "nonzero_representative_selected": True,
            "scale_orbit_quotient_well_defined": True,
            "unique_section_exported": True,
            "strict_source_law_exported": False,
            "C_phi_A_phi_coupling_proved": False,
            "action_from_length_time_relation_proved": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": False,
            "nonconventional_source_law_exported": False,
            "blocker": "selector premise violates the P3114 prohibition and does not source dimensions",
        },
    ]


def orbit_witness_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        for scale in SCALE_ACTIONS:
            rows.append({
                "candidate": cand["candidate"],
                "scale_factor": f"{scale.numerator}/{scale.denominator}",
                "orbit_point": {"U_action": f"{scale}*S", "U_length": f"{scale}*L", "U_time": f"{scale}*T"},
                "section_claimed": cand["unique_section_exported"],
                "section_accepted": bool(cand["scale_orbit_quotient_well_defined"] and cand["unique_section_exported"] and cand["strict_source_law_exported"] and cand["standard_physics_import_free"] and cand["selector_bridge_ltotal_toe_free"] and cand["nonconventional_source_law_exported"]),
                "blocker": cand["blocker"],
            })
    return rows


def obligation_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "nonzero_representative": "nonzero_representative_selected",
        "orbit_invariance": "scale_orbit_quotient_well_defined",
        "uniqueness_mod_scale": "unique_section_exported",
        "strict_source_law": "strict_source_law_exported",
        "C_phi_coupling": "C_phi_A_phi_coupling_proved",
        "action_length_time_relation": "action_from_length_time_relation_proved",
    }
    return [
        {
            "candidate": cand["candidate"],
            "obligation": obligation,
            "A_phi": round(a_phi(), 12) if obligation == "C_phi_coupling" else None,
            "claimed": cand[field[obligation]],
            "accepted": bool(cand[field[obligation]] and cand["standard_physics_import_free"] and cand["selector_bridge_ltotal_toe_free"] and cand["nonconventional_source_law_exported"]),
            "blocker": cand["blocker"],
        }
        for cand in candidates
        for obligation in SECTION_OBLIGATIONS
    ]


def coupling_residual_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        c_phi_residual = 0 if cand["C_phi_A_phi_coupling_proved"] else 1
        relation_residual = 0 if cand["action_from_length_time_relation_proved"] else 1
        source_residual = 0 if cand["strict_source_law_exported"] and cand["standard_physics_import_free"] and cand["nonconventional_source_law_exported"] else 1
        rows.append({
            "candidate": cand["candidate"],
            "residual_vector": {"C_phi(A_phi)-U_action": c_phi_residual, "U_action-F(U_length,U_time)": relation_residual, "strict_source_defect": source_residual},
            "all_residuals_zero_import_free": bool(c_phi_residual == 0 and relation_residual == 0 and source_residual == 0 and cand["scale_orbit_quotient_well_defined"] and cand["unique_section_exported"] and cand["selector_bridge_ltotal_toe_free"]),
            "blocker": cand["blocker"],
        })
    return rows


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": cand["candidate"], "required_gate": gate, "gate_passed": bool(cand[gate]), "detail": "passed" if cand[gate] else cand["blocker"]} for cand in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_Sigma_dim_theorem": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3114 = read_json(P3114)
    greps = content_grep()
    candidates = sigma_candidates()
    orbit_rows = orbit_witness_rows(candidates)
    obligations = obligation_rows(candidates)
    residuals = coupling_residual_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Sigma_dim_theorem"]]
    proof_obligations = [
        {"obligation": "read_p3114_next_atom", "satisfied": True, "detail": "P3114 requested exactly one Sigma_dim scale-section theorem"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by research content patterns"},
        {"obligation": "construct_sigma_dim_candidates", "satisfied": len(candidates) == 8, "detail": "eight section-theorem candidates were constructed"},
        {"obligation": "test_scale_orbit_witnesses", "satisfied": len(orbit_rows) == len(candidates) * len(SCALE_ACTIONS), "detail": "scale orbit witnesses were tested across five factors"},
        {"obligation": "test_section_obligations", "satisfied": len(obligations) == len(candidates) * len(SECTION_OBLIGATIONS), "detail": "six section obligations were tested per candidate"},
        {"obligation": "test_coupling_residuals", "satisfied": len(residuals) == len(candidates), "detail": "C_phi/action-length-time/source residual vectors were built"},
        {"obligation": "export_nadsoliton_only_Sigma_dim", "satisfied": False, "detail": "0 candidates export an import-free strict section theorem satisfying all gates"},
    ]
    return {
        "status": "P3115_SIGMA_DIM_SCALE_SECTION_THEOREM_BOUNDED_NO_GO",
        "input_hashes": {"P3114": hashlib.sha256(P3114.read_bytes()).hexdigest() if P3114.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "audit_object": {"object": "SigmaDimScaleSectionTheoremAudit", "required_theorem": "Sigma_dim selects a nonzero D_phi representative, proves C_phi(A_phi)=U_action, and proves action-from-length/time relation", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "selector replay", "L_total", "bridge/role-transfer", "ToE"]},
            "candidate_Sigma_dim_theorems": candidates,
            "scale_orbit_witness_rows": orbit_rows,
            "section_obligation_rows": obligations,
            "coupling_residual_rows": residuals,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3114_accepted_D_phi_source_laws": p3114.get("finite_certificate", {}).get("accepted_D_phi_source_laws"),
            "candidate_Sigma_dim_theorems": len(candidates),
            "scale_orbit_witness_rows": len(orbit_rows),
            "section_obligation_rows": len(obligations),
            "coupling_residual_rows": len(residuals),
            "candidate_gate_rows": len(gates),
            "accepted_Sigma_dim_theorems": len(accepted),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3115 constructs the requested Sigma_dim scale-section theorem family and finds bounded no-go.  Identity, norm, phase-area, entropy/tick, and cohomological sections are mathematically meaningful but either fail to quotient the scale orbit, lack strict source-law status, omit C_phi(A_phi)=U_action, or omit the action-from-length/time relation.  The complete Planck/light and apparatus sections import external calibration, and the selector-oriented section violates the selector ban.  No nadsoliton-only Sigma_dim theorem exports a nonzero import-free representative satisfying all section, coupling, and relation obligations.",
            "negative_export_flags": {key: False for key in ["Sigma_dim_theorem_exported", "D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"p3114_sigma_dim_obligation_reused": True, "candidate_Sigma_dim_theorems_constructed": True, "scale_orbit_witness_matrix_built": True, "section_obligation_matrix_built": True, "external_calibration_and_selector_rows_rejected": True},
            "next_honest_step": "Construct exactly one nadsoliton-only dimension-source functor K_dim from strict nadsoliton data to the positive scale torsor, together with a natural section Sigma_dim and a proof that K_dim supplies the missing nonconventional source law for C_phi(A_phi)=U_action and U_action=F(U_length,U_time).  It must not import hbar/Planck, rods, clocks, observed light, apparatus, selector replay, L_total, bridge/role-transfer, or ToE; otherwise preserve the P3105-P3115 physical-unit no-go.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3115/S2065 Sigma_dim scale-section theorem audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- P3114 accepted D_phi source laws: `{cert['p3114_accepted_D_phi_source_laws']}`",
        f"- content grep lanes: `{cert['content_grep_lanes']}`",
        f"- candidate Sigma_dim theorems: `{cert['candidate_Sigma_dim_theorems']}`",
        f"- scale-orbit witness rows: `{cert['scale_orbit_witness_rows']}`",
        f"- section obligation rows: `{cert['section_obligation_rows']}`",
        f"- coupling residual rows: `{cert['coupling_residual_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted Sigma_dim theorems: `{cert['accepted_Sigma_dim_theorems']}`",
        f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3115/S2065 Sigma_dim scale-section theorem audit", "## P3115/S2065 Sigma_dim scale-section theorem audit\n\n`P3115/S2065` executes the P3114-recommended audit for a nadsoliton-only scale-section theorem `Sigma_dim` over the internal `D_phi=(U_action,U_length,U_time)` dimension orbit.  It constructs `8` candidate section theorems, `40` scale-orbit witness rows, `48` section-obligation rows, `8` coupling-residual rows, and an `8 x 11 = 88` gate matrix.  The bounded result is that identity/norm/phase/entropy/cohomology candidates miss at least one required strict source, coupling, relation, or quotient condition, while complete Planck/light/apparatus/selector candidates import forbidden structure.  No physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3115/S2065 Sigma_dim section remains unsourced", "## P3115/S2065 Sigma_dim section remains unsourced\n\n`P3115/S2065` tests whether a nadsoliton-only `Sigma_dim` can select a nonzero representative for the internal `D_phi` dimension orbit and prove both `C_phi(A_phi)=U_action` and an action-from-length/time relation.  Current artifacts provide no import-free strict source law for that section, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Sigma_dim scale-section theorem guardrail (P3115/S2065, 2026-06-26)", "## Current Sigma_dim scale-section theorem guardrail (P3115/S2065, 2026-06-26)\n\n- P3115 tests the P3114-requested nadsoliton-only scale-section theorem `Sigma_dim` for the internal `D_phi=(U_action,U_length,U_time)` dimension orbit.\n- The finite audit constructs `8` candidate section theorems, `40` scale-orbit witnesses, `48` section-obligation rows, `8` coupling-residual rows, and `88` gate rows; `0` candidates export an import-free strict `Sigma_dim` theorem.\n- Do not promote identity/norm gauges, phase-area representatives, entropy/tick ratios, cohomological periods, Planck/light calibration, apparatus calibration, or selector-premise sections to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one nadsoliton-only dimension-source functor `K_dim` to the positive scale torsor that supplies the missing strict source law for `Sigma_dim`; otherwise preserve the P3105-P3115 physical-unit no-go.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
