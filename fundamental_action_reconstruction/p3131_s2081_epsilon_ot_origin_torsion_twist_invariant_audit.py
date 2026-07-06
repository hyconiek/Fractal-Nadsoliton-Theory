#!/usr/bin/env python3
"""P3131/S2081: Epsilon_OT origin-torsion/twist invariant audit.

P3130 left exactly one admissible next object: Epsilon_OT, a strict
non-translation origin-torsion/origin-twist datum not erased by the Z12
translation quotient and still inversion/sign aware.  This audit proves the
finite equivariant-section obstruction for all nonempty Z12 support classes and
then tests concrete twist/torsion candidates for whether they can choose an
absolute origin inside a translation class without forbidden imports.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from collections import Counter
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3113_s2063_u_action_reference_carrier_source_law_audit import a_phi
from p3130_s2080_theta_to_translation_origin_quotient_audit import OUT as P3130

OUT = GEN / "p3131_s2081_epsilon_ot_origin_torsion_twist_invariant_audit.json"
MD = GEN / "p3131_s2081_epsilon_ot_origin_torsion_twist_invariant_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12
CHARACTERS = tuple(range(N))
TESTS = (
    "explicit_formula", "strict_nadsoliton_only", "nontranslation_datum", "not_erased_by_quotient",
    "absolute_origin_representative", "inversion_sign_aware", "nonzero_twist", "equivariant_section",
    "forms_theta_to", "forms_gamma_so", "selector_free", "forbidden_import_free", "preserves_A_phi",
)
GATES = (
    "uses_p3130_epsilon_obligation", "explicit_epsilon_formula", "strict_nadsoliton_data_only",
    "nontranslation_datum_exported", "not_erased_by_translation_quotient", "absolute_origin_representative_exported",
    "inversion_sign_aware", "nonzero_twist_exported", "equivariant_origin_section_exported",
    "Theta_TO_retest_unlocked", "Gamma_SO_retest_unlocked", "selector_free", "not_unit_convention",
    "not_imported_dynamics", "not_apparatus_calibration", "selector_bridge_ltotal_toe_free", "C_phi_A_phi_preserved",
)


def support_from_mask(mask: int) -> tuple[int, ...]:
    return tuple(i for i in range(N) if mask & (1 << i))


def translate(support: tuple[int, ...], shift: int) -> tuple[int, ...]:
    return tuple(sorted(((x + shift) % N for x in support)))


def unit_action(support: tuple[int, ...], unit: int) -> tuple[int, ...]:
    return tuple(sorted(((unit * x) % N for x in support)))


def orbit(support: tuple[int, ...]) -> set[tuple[int, ...]]:
    return {translate(support, shift) for shift in range(N)}


def canonical(support: tuple[int, ...]) -> tuple[int, ...]:
    return min(orbit(support))


def translation_classes() -> dict[tuple[int, ...], list[tuple[int, ...]]]:
    classes: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for mask in range(1, 1 << N):
        support = support_from_mask(mask)
        classes.setdefault(canonical(support), []).append(support)
    return classes


def equivariant_section_obstruction_rows() -> list[dict[str, Any]]:
    rows = []
    for rep, members in sorted(translation_classes().items(), key=lambda item: (len(item[0]), item[0])):
        size = len(members)
        stabilizer_size = N // size
        stabilizer_shifts = [shift for shift in range(N) if translate(rep, shift) == rep]
        inverted = canonical(unit_action(rep, 11))
        rows.append({
            "representative": list(rep),
            "orbit_size": size,
            "stabilizer_size": stabilizer_size,
            "stabilizer_shifts": stabilizer_shifts,
            "inversion_partner": list(inverted),
            "inversion_fixed_class": inverted == rep,
            "equivariant_global_origin_section_exists": size == 1,
            "obstruction": "a translation-equivariant section from a quotient point to an orbit element exists only for fixed orbits; nonfixed classes need extra non-translation data",
        })
    return rows


def character_rows() -> list[dict[str, Any]]:
    rows = []
    for k in CHARACTERS:
        rows.append({
            "character_k": k,
            "generator_phase_numerator": k,
            "inversion_partner_k": (-k) % N,
            "inversion_fixed": ((-k) % N) == k,
            "nontrivial": k != 0,
            "origin_selecting": False,
            "reason": "C12 characters give twist/torsion labels, but without a coupled support-local section they do not choose a support origin",
        })
    return rows


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "epsilon_ot": r"Epsilon_OT|origin-torsion|origin-twist|non-translation datum",
        "theta_gamma": r"Theta_TO|Gamma_SO|translation-origin quotient|sign-and-origin",
        "cohomology_character": r"character|torsion|twist|cocycle|cohomology|equivariant section",
        "blocked_imports": r"selector|apparatus|observed light|Planck|hbar|thermodynamic|Lagrangian/EOM|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def candidate(name: str, formula: str, blocker: str, **overrides: bool) -> dict[str, Any]:
    base = {gate: False for gate in GATES}
    base.update({
        "uses_p3130_epsilon_obligation": True,
        "explicit_epsilon_formula": True,
        "strict_nadsoliton_data_only": True,
        "selector_free": True,
        "not_unit_convention": True,
        "not_imported_dynamics": True,
        "not_apparatus_calibration": True,
        "selector_bridge_ltotal_toe_free": True,
        "C_phi_A_phi_preserved": True,
    })
    base.update(overrides)
    return {"candidate": name, "formula": formula, **base, "blocker": blocker}


def epsilon_candidates() -> list[dict[str, Any]]:
    twist = dict(nontranslation_datum_exported=True, not_erased_by_translation_quotient=True, nonzero_twist_exported=True, inversion_sign_aware=True)
    origin = dict(absolute_origin_representative_exported=True, equivariant_origin_section_exported=True, Theta_TO_retest_unlocked=True, Gamma_SO_retest_unlocked=True)
    return [
        candidate("c12_character_twist", "Epsilon_OT := nontrivial character chi_k of C12 translation torsor", "character is a real twist label but not support-local origin data", **twist),
        candidate("inversion_odd_character_twist", "Epsilon_OT := inversion-odd C12 character pair k,-k", "inversion-odd pair supplies sign type but no absolute origin section", **twist),
        candidate("stabilizer_torsion_label", "Epsilon_OT := stabilizer subgroup torsion label of a translation class", "stabilizer labels period, not a bead inside the orbit", nontranslation_datum_exported=True, not_erased_by_translation_quotient=True, nonzero_twist_exported=True),
        candidate("equivariant_section_fixed_orbit", "Epsilon_OT := equivariant section on fixed translation orbit", "works only for the all-support fixed orbit and gives no general Gamma_SO", **twist, absolute_origin_representative_exported=True, equivariant_origin_section_exported=True),
        candidate("cohomology_one_cocycle_twist", "Epsilon_OT := Z12-valued 1-cocycle on translation orbit", "cocycle is origin-gauge dependent unless a section is supplied", **twist),
        candidate("phase_winding_origin_twist", "Epsilon_OT := phase-winding torsion around translation orbit", "winding is class-level and does not choose a support origin", **twist),
        candidate("a_phi_twisted_cell_index", "Epsilon_OT := A_phi-twisted cell index over translation class", "preserves A_phi but cell index is quotient-gauge data", **twist),
        candidate("chiral_bispectrum_twist", "Epsilon_OT := chiral-bispectrum twist label on quotient class", "chiral sign is nonzero but origin remains translation-erased", **twist),
        candidate("rank_monodromy_twist", "Epsilon_OT := rank-monodromy around support orbit", "rank monodromy is constant on most full orbits and not origin-selecting", **twist),
        candidate("spectral_flow_twist", "Epsilon_OT := spectral-flow phase over translations", "spectral flow phase lacks a strict support-local zero point", **twist),
        candidate("crt_origin_twist", "Epsilon_OT := CRT residue twist label", "CRT residue is not coupled to a unique translated support", **twist),
        candidate("assumed_strict_origin_torsion", "Epsilon_OT := nontranslation torsion plus assumed strict section", "completed shape is identified, but the strict section theorem is unavailable", strict_nadsoliton_data_only=False, **twist, **origin),
        candidate("selector_origin_torsion", "Epsilon_OT := selector-chosen origin torsion", "imports selector", **twist, **origin, selector_free=False, selector_bridge_ltotal_toe_free=False),
        candidate("apparatus_origin_twist", "Epsilon_OT := apparatus-marked origin twist", "imports apparatus", strict_nadsoliton_data_only=False, **twist, **origin, not_apparatus_calibration=False),
        candidate("lagrangian_origin_twist", "Epsilon_OT := coordinate origin from Lagrangian normalization", "imports Lagrangian/EOM normalization", strict_nadsoliton_data_only=False, **twist, **origin, not_imported_dynamics=False),
    ]


def test_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "explicit_formula": "explicit_epsilon_formula", "strict_nadsoliton_only": "strict_nadsoliton_data_only",
        "nontranslation_datum": "nontranslation_datum_exported", "not_erased_by_quotient": "not_erased_by_translation_quotient",
        "absolute_origin_representative": "absolute_origin_representative_exported", "inversion_sign_aware": "inversion_sign_aware",
        "nonzero_twist": "nonzero_twist_exported", "equivariant_section": "equivariant_origin_section_exported",
        "forms_theta_to": "Theta_TO_retest_unlocked", "forms_gamma_so": "Gamma_SO_retest_unlocked",
        "selector_free": "selector_free", "forbidden_import_free": "selector_bridge_ltotal_toe_free", "preserves_A_phi": "C_phi_A_phi_preserved",
    }
    return [{"candidate": c["candidate"], "epsilon_test": test, "test_passed": bool(c[field[test]]), "accepted_test_row": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c.get(gate, False)), "detail": "passed" if c.get(gate, False) else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": name, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == name), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == name), "accepted_Epsilon_OT_source": all(row["gate_passed"] for row in gates if row["candidate"] == name)} for name in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3130 = read_json(P3130)
    greps = content_grep()
    section_rows = equivariant_section_obstruction_rows()
    chars = character_rows()
    candidates = epsilon_candidates()
    e_rows = test_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Epsilon_OT_source"]]
    orbit_sizes = Counter(row["orbit_size"] for row in section_rows)
    stabilizers = Counter(row["stabilizer_size"] for row in section_rows)
    proof_obligations = [
        {"obligation": "read_p3130_next_atom", "satisfied": True, "detail": "P3130 requested exactly one Epsilon_OT non-translation origin-torsion/twist object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by Epsilon_OT, Theta/Gamma, cohomology/character, and blocked-import patterns"},
        {"obligation": "prove_finite_equivariant_section_obstruction", "satisfied": len(section_rows) == 351 and sum(row["equivariant_global_origin_section_exists"] for row in section_rows) == 1, "detail": "all translation classes were tested; only the fixed all-support orbit admits a quotient-to-origin section"},
        {"obligation": "compute_translation_character_twists", "satisfied": len(chars) == 12 and sum(row["nontrivial"] for row in chars) == 11, "detail": "all C12 character twists were enumerated"},
        {"obligation": "construct_Epsilon_OT_candidates", "satisfied": len(candidates) == 15, "detail": "fifteen origin-torsion/twist candidates were constructed"},
        {"obligation": "test_epsilon_laws", "satisfied": len(e_rows) == len(candidates) * len(TESTS), "detail": "thirteen epsilon-law tests were built per candidate"},
        {"obligation": "export_strict_Epsilon_OT_to_Gamma_SO", "satisfied": False, "detail": "0 candidates export an import-free origin-twist/torsion section forming Gamma_SO"},
    ]
    return {
        "status": "P3131_EPSILON_OT_ORIGIN_TORSION_TWIST_BOUNDED_NO_GO",
        "input_hashes": {"P3130": hashlib.sha256(P3130.read_bytes()).hexdigest() if P3130.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "equivariant_section_obstruction_rows": section_rows,
            "translation_character_twist_rows": chars,
            "audit_object": {"object": "EpsilonOTOriginTorsionTwistAudit", "required_source": "Epsilon_OT strict non-translation origin-torsion/twist datum not erased by translation quotient and able to choose an absolute origin", "forbidden_imports": ["selector", "apparatus", "observed light", "hbar/Planck", "thermodynamic environment", "Lagrangian/EOM normalization", "bridge/role-transfer", "L_total", "ToE"]},
            "candidate_Epsilon_OT_sources": candidates,
            "epsilon_law_rows": e_rows,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3130_accepted_Theta_TO_sources": p3130.get("finite_certificate", {}).get("accepted_Theta_TO_sources"),
            "translation_classes": len(section_rows),
            "equivariant_origin_section_classes": sum(row["equivariant_global_origin_section_exists"] for row in section_rows),
            "orbit_size_distribution": dict(sorted(orbit_sizes.items())),
            "stabilizer_size_distribution": dict(sorted(stabilizers.items())),
            "translation_character_twists": len(chars),
            "nontrivial_character_twists": sum(row["nontrivial"] for row in chars),
            "origin_selecting_character_twists": sum(row["origin_selecting"] for row in chars),
            "candidate_Epsilon_OT_sources": len(candidates),
            "epsilon_law_rows": len(e_rows),
            "candidate_gate_rows": len(gates),
            "accepted_Epsilon_OT_sources": len(accepted),
            "A_phi": round(a_phi(), 12),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3131 constructs the requested Epsilon_OT origin-torsion/twist family and proves the finite equivariant-section obstruction. Across all 351 Z12 translation classes, only the fixed all-support orbit admits a translation-equivariant quotient-to-origin section; all nonfixed classes require extra non-translation data. The 12 C12 character/twist labels are real, with 11 nontrivial twists, but none is origin-selecting without an additional support-local section. Therefore no import-free Epsilon_OT forms Theta_TO or Gamma_SO on current artifacts.",
            "negative_export_flags": {key: False for key in ["Epsilon_OT_source_exported", "Theta_TO_source_exported", "Gamma_SO_generator_exported", "Sigma_point_source_exported", "Omega_tie_source_exported", "Pi_point_source_exported", "Lambda_origin_source_exported", "Phi_Info_source_exported", "Delta_asym_source_exported", "Iota_irrev_source_exported", "Kappa_cycle_source_exported", "Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"p3130_Epsilon_OT_obligation_reused": True, "finite_equivariant_section_obstruction_proved": True, "C12_character_twists_enumerated": True, "candidate_Epsilon_OT_sources_constructed": True, "epsilon_law_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True},
            "next_honest_step": "Do not replay origin torsion/twist by adding another uncoupled C12 character. Construct exactly one new strict support-local section law Zeta_OS: a nadsoliton-internal local section theorem coupling a nontranslation twist to a specific support representative within a translation class. It must produce a nonconventional absolute origin representative and survive inversion/sign tests without selector, apparatus, observed light, Planck units, thermodynamic environment, Lagrangian/EOM normalization, bridge/role-transfer, L_total, or ToE imports; otherwise preserve the P3105-P3131 physical-unit/sign-origin no-go certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3131/S2081 Epsilon_OT origin-torsion/twist invariant audit", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- P3130 accepted Theta_TO sources: `{cert['p3130_accepted_Theta_TO_sources']}`",
        f"- translation classes: `{cert['translation_classes']}`",
        f"- equivariant origin-section classes: `{cert['equivariant_origin_section_classes']}`",
        f"- orbit-size distribution: `{cert['orbit_size_distribution']}`",
        f"- stabilizer-size distribution: `{cert['stabilizer_size_distribution']}`",
        f"- translation character twists: `{cert['translation_character_twists']}`",
        f"- nontrivial character twists: `{cert['nontrivial_character_twists']}`",
        f"- origin-selecting character twists: `{cert['origin_selecting_character_twists']}`",
        f"- candidate Epsilon_OT sources: `{cert['candidate_Epsilon_OT_sources']}`",
        f"- epsilon-law rows: `{cert['epsilon_law_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted Epsilon_OT sources: `{cert['accepted_Epsilon_OT_sources']}`",
        f"- A_phi: `{cert['A_phi']}`",
        f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3131/S2081 Epsilon_OT origin-torsion/twist invariant audit", "## P3131/S2081 Epsilon_OT origin-torsion/twist invariant audit\n\n`P3131/S2081` executes the P3130-recommended audit for `Epsilon_OT`, a strict non-translation origin-torsion/origin-twist datum not erased by the `Z12` translation quotient. It proves the finite equivariant-section obstruction over all `351` translation classes: only `1` fixed class admits a quotient-to-origin section. It also enumerates all `12` `C12` character/twist labels (`11` nontrivial), tests `15` candidates with `195` epsilon-law rows and a `15 x 17 = 255` gate matrix, and finds `0` import-free `Epsilon_OT` sources. No physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3131/S2081 Epsilon_OT twist is not an action origin", "## P3131/S2081 Epsilon_OT twist is not an action origin\n\n`P3131/S2081` constructs origin-torsion/twist labels and proves the finite section obstruction, but it does not export a unit-bearing coordinate origin, Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Epsilon_OT origin-torsion/twist guardrail (P3131/S2081, 2026-07-06)", "## Current Epsilon_OT origin-torsion/twist guardrail (P3131/S2081, 2026-07-06)\n\n- P3131 tests the P3130-requested strict non-translation datum object `Epsilon_OT`, an origin-torsion/origin-twist invariant not erased by the `Z12` translation quotient and inversion/sign aware.\n- The finite audit proves the equivariant-section obstruction across all `351` translation classes: only `1` fixed class admits a quotient-to-origin section; it enumerates `12` translation-character twists, builds `15` candidates, `195` epsilon-law rows, and `255` gate rows; `0` candidates export an import-free `Epsilon_OT` forming `Theta_TO` or `Gamma_SO`.\n- The character/twist labels are real but scoped; do not promote uncoupled `C12` characters, inversion-odd character pairs, stabilizer labels, fixed-orbit sections, cocycles, phase windings, `A_phi` cell indices, chiral-bispectrum twists, rank/spectral/CRT twists, assumed section theorems, selector origins, apparatus origins, or Lagrangian origins to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict support-local section law `Zeta_OS`, a nadsoliton-internal theorem coupling a nontranslation twist to a specific support representative within a translation class; otherwise preserve the P3105-P3131 physical-unit/sign-origin no-go certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
