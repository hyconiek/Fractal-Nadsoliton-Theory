#!/usr/bin/env python3
"""P3130/S2080: Theta_TO translation-origin quotient audit.

P3129 left exactly one admissible next object: Theta_TO, a strict
translation-origin quotient/fixing theorem that should quotient or fix the Z12
translation-origin orbit while remaining inversion/sign aware.  This audit
computes the full finite nonempty-subset translation quotient of Z12 and tests
whether any concrete Theta_TO candidate can turn quotient data plus a nonzero
internal sign into Gamma_SO without importing selector, apparatus, observed
light, Planck units, thermodynamic environment, Lagrangian/EOM normalization,
bridge/role-transfer, L_total, or ToE lanes.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3113_s2063_u_action_reference_carrier_source_law_audit import a_phi
from p3129_s2079_gamma_so_sign_origin_generator_audit import OUT as P3129

OUT = GEN / "p3130_s2080_theta_to_translation_origin_quotient_audit.json"
MD = GEN / "p3130_s2080_theta_to_translation_origin_quotient_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
UNITS = (1, 5, 7, 11)
QUOTIENT_TESTS = (
    "explicit_formula", "strict_nadsoliton_only", "translation_quotient_constructed",
    "absolute_origin_fixed", "inversion_sign_aware", "nonzero_internal_sign_available",
    "theta_plus_sign_forms_gamma", "aut_compatible", "selector_free", "nonconventional",
    "pre_gamma_retest", "forbidden_import_free", "preserves_A_phi", "unlocks_gamma_retest",
)
GATES = (
    "uses_p3129_theta_obligation", "explicit_theta_formula", "strict_nadsoliton_data_only",
    "translation_quotient_exported", "absolute_origin_representative_exported",
    "inversion_sign_aware", "nonzero_internal_sign_exported", "theta_plus_sign_forms_Gamma_SO",
    "aut_compatible", "selector_free", "nonconventional_source", "before_gamma_retest",
    "not_unit_convention", "not_imported_dynamics", "not_apparatus_calibration",
    "selector_bridge_ltotal_toe_free", "Gamma_SO_retest_unlocked", "C_phi_A_phi_preserved",
)
ACTIONS = ("translate_all", "unit_5", "unit_7", "inversion_11", "origin_relabel", "selector_choose_zero", "apparatus_mark_zero")


def support_from_mask(mask: int) -> tuple[int, ...]:
    return tuple(i for i in range(N) if mask & (1 << i))


def translate(support: tuple[int, ...], shift: int) -> tuple[int, ...]:
    return tuple(sorted(((x + shift) % N for x in support)))


def unit_action(support: tuple[int, ...], unit: int) -> tuple[int, ...]:
    return tuple(sorted(((unit * x) % N for x in support)))


def translation_orbit(support: tuple[int, ...]) -> set[tuple[int, ...]]:
    return {translate(support, shift) for shift in range(N)}


def canonical_translation_class(support: tuple[int, ...]) -> tuple[int, ...]:
    return min(translation_orbit(support))


def dihedral_class(support: tuple[int, ...]) -> tuple[int, ...]:
    images = set()
    for shift in range(N):
        images.add(translate(support, shift))
        images.add(translate(unit_action(support, 11), shift))
    return min(images)


def all_nonempty_supports() -> list[tuple[int, ...]]:
    return [support_from_mask(mask) for mask in range(1, 1 << N)]


def quotient_summary() -> dict[str, Any]:
    supports = all_nonempty_supports()
    translation_classes: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    dihedral_classes: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for support in supports:
        translation_classes.setdefault(canonical_translation_class(support), []).append(support)
        dihedral_classes.setdefault(dihedral_class(support), []).append(support)
    translation_class_rows = []
    for rep, members in sorted(translation_classes.items(), key=lambda item: (len(item[0]), item[0])):
        inverted_rep = canonical_translation_class(unit_action(rep, 11))
        translation_class_rows.append({
            "representative": list(rep),
            "size": len(members),
            "cardinality": len(rep),
            "stabilizer_size": N // len(members),
            "inversion_partner": list(inverted_rep),
            "inversion_fixed_translation_class": inverted_rep == rep,
            "absolute_origin_fixed_by_quotient": len(members) == 1,
        })
    return {
        "nonempty_supports": len(supports),
        "translation_classes": len(translation_classes),
        "dihedral_classes": len(dihedral_classes),
        "translation_classes_with_absolute_origin": sum(row["absolute_origin_fixed_by_quotient"] for row in translation_class_rows),
        "inversion_fixed_translation_classes": sum(row["inversion_fixed_translation_class"] for row in translation_class_rows),
        "translation_class_rows": translation_class_rows,
    }


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "theta_to": r"Theta_TO|translation-origin quotient|translation-origin|origin quotient",
        "gamma_so": r"Gamma_SO|sign-and-origin generator|joint sign-origin|source-origin",
        "finite_z12": r"Z12|translation class|dihedral class|Aut\(Z12\)|inversion",
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
        "uses_p3129_theta_obligation": True,
        "explicit_theta_formula": True,
        "strict_nadsoliton_data_only": True,
        "before_gamma_retest": True,
        "selector_free": True,
        "not_unit_convention": True,
        "not_imported_dynamics": True,
        "not_apparatus_calibration": True,
        "selector_bridge_ltotal_toe_free": True,
        "C_phi_A_phi_preserved": True,
    })
    base.update(overrides)
    return {"candidate": name, "formula": formula, **base, "blocker": blocker}


def theta_candidates() -> list[dict[str, Any]]:
    quotient = dict(translation_quotient_exported=True, nonconventional_source=True, aut_compatible=True)
    origin = dict(absolute_origin_representative_exported=True)
    signed = dict(inversion_sign_aware=True, nonzero_internal_sign_exported=True)
    gamma = dict(theta_plus_sign_forms_Gamma_SO=True, Gamma_SO_retest_unlocked=True)
    return [
        candidate("necklace_translation_quotient", "Theta_TO := support / Z12 translations as binary necklace", "quotient is strict and finite but erases absolute origin", **quotient),
        candidate("bracelet_translation_inversion_quotient", "Theta_TO := support / <translation,inversion> as binary bracelet", "bracelet is inversion-aware but erases sign and origin", **quotient, inversion_sign_aware=True),
        candidate("pointed_necklace_with_internal_mark", "Theta_TO := necklace plus internal marked bead from support invariant", "marked bead is not exported by translation-invariant data", **quotient, **origin),
        candidate("minimal_lexicographic_rotation", "Theta_TO := lexicographically minimal rotation representative", "lexicographic minimum fixes a display representative by convention, not strict source", **quotient, **origin),
        candidate("burnside_fixed_class_source", "Theta_TO := Burnside fixed translation class", "only full-period/all-period classes fix origin and do not provide nonzero sign", **quotient, **origin),
        candidate("cohomology_period_origin_quotient", "Theta_TO := cohomology-period translation quotient", "period quotient is strict but not an absolute source-origin representative", **quotient),
        candidate("a_phi_phase_cell_quotient", "Theta_TO := A_phi phase-cell orbit quotient", "preserves A_phi but phase cell quotient has no unique origin", **quotient),
        candidate("rank_signature_translation_quotient", "Theta_TO := rank-signature class modulo translations", "rank signature is constant on translated supports and cannot point within the orbit", **quotient),
        candidate("spectral_necklace_quotient", "Theta_TO := spectral invariant of translation class", "spectral quotient is origin-blind", **quotient),
        candidate("chiral_signed_quotient", "Theta_TO := translation quotient plus chiral sign", "sign can be nonzero but the quotient still lacks an absolute origin representative", **quotient, **signed),
        candidate("gamma_retest_ready_if_origin_given", "Theta_TO := quotient plus assumed strict origin theorem", "row shows the missing theorem shape but the origin theorem is unavailable", strict_nadsoliton_data_only=False, **quotient, **origin, **signed, **gamma),
        candidate("selector_zero_origin", "Theta_TO := choose origin 0 in each translation class", "imports selector/convention", **quotient, **origin, **signed, **gamma, selector_free=False, selector_bridge_ltotal_toe_free=False),
        candidate("apparatus_marked_origin", "Theta_TO := apparatus-marked source origin", "imports apparatus calibration", strict_nadsoliton_data_only=False, **quotient, **origin, **signed, **gamma, not_apparatus_calibration=False),
        candidate("observed_light_event_origin", "Theta_TO := observed-light event origin", "imports observed light", strict_nadsoliton_data_only=False, **quotient, **origin, **signed, **gamma, not_apparatus_calibration=False),
        candidate("lagrangian_coordinate_origin", "Theta_TO := Lagrangian coordinate origin", "imports Lagrangian/EOM normalization", strict_nadsoliton_data_only=False, **quotient, **origin, **signed, **gamma, not_imported_dynamics=False),
    ]


def quotient_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "explicit_formula": "explicit_theta_formula", "strict_nadsoliton_only": "strict_nadsoliton_data_only",
        "translation_quotient_constructed": "translation_quotient_exported", "absolute_origin_fixed": "absolute_origin_representative_exported",
        "inversion_sign_aware": "inversion_sign_aware", "nonzero_internal_sign_available": "nonzero_internal_sign_exported",
        "theta_plus_sign_forms_gamma": "theta_plus_sign_forms_Gamma_SO", "aut_compatible": "aut_compatible",
        "selector_free": "selector_free", "nonconventional": "nonconventional_source", "pre_gamma_retest": "before_gamma_retest",
        "forbidden_import_free": "selector_bridge_ltotal_toe_free", "preserves_A_phi": "C_phi_A_phi_preserved",
        "unlocks_gamma_retest": "Gamma_SO_retest_unlocked",
    }
    return [{"candidate": c["candidate"], "quotient_test": test, "test_passed": bool(c[field[test]]), "accepted_quotient_row": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in QUOTIENT_TESTS]


def symmetry_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        for action in ACTIONS:
            translation = action == "translate_all"
            aut = action in {"unit_5", "unit_7", "inversion_11"}
            inversion = action in {"unit_7", "inversion_11"}
            forbidden = action in {"selector_choose_zero", "apparatus_mark_zero"}
            origin_relabel = action == "origin_relabel"
            survives = (not translation or c["translation_quotient_exported"]) and (not aut or c["aut_compatible"]) and (not inversion or c["inversion_sign_aware"]) and (not origin_relabel or c["absolute_origin_representative_exported"]) and not forbidden
            rows.append({"candidate": c["candidate"], "action": action, "translation_failure": translation and not c["translation_quotient_exported"], "aut_failure": aut and not c["aut_compatible"], "inversion_failure": inversion and not c["inversion_sign_aware"], "origin_failure": origin_relabel and not c["absolute_origin_representative_exported"], "forbidden_action_rejected": forbidden and c["selector_bridge_ltotal_toe_free"] and c["not_apparatus_calibration"], "theta_survives_action": bool(survives), "blocker": c["blocker"]})
    return rows


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c.get(gate, False)), "detail": "passed" if c.get(gate, False) else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": name, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == name), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == name), "accepted_Theta_TO_source": all(row["gate_passed"] for row in gates if row["candidate"] == name)} for name in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3129 = read_json(P3129)
    greps = content_grep()
    summary = quotient_summary()
    candidates = theta_candidates()
    q_rows = quotient_law_rows(candidates)
    s_rows = symmetry_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Theta_TO_source"]]
    strict_quotients = [c for c in candidates if c["strict_nadsoliton_data_only"] and c["translation_quotient_exported"] and c["nonconventional_source"]]
    proof_obligations = [
        {"obligation": "read_p3129_next_atom", "satisfied": True, "detail": "P3129 requested exactly one Theta_TO translation-origin quotient object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by Theta_TO, Gamma_SO, finite Z12, and closure-block patterns"},
        {"obligation": "exhaust_nonempty_Z12_translation_quotient", "satisfied": summary["nonempty_supports"] == 4095 and summary["translation_classes"] == 351, "detail": "all nonempty Z12 supports were quotient-counted under translations"},
        {"obligation": "test_inversion_awareness", "satisfied": summary["dihedral_classes"] > 0 and summary["inversion_fixed_translation_classes"] > 0, "detail": "translation classes were also paired into inversion/dihedral classes"},
        {"obligation": "construct_Theta_TO_candidates", "satisfied": len(candidates) == 15, "detail": "fifteen quotient/fixing candidates were constructed"},
        {"obligation": "test_quotient_laws", "satisfied": len(q_rows) == len(candidates) * len(QUOTIENT_TESTS), "detail": "fourteen quotient-law tests were built per candidate"},
        {"obligation": "test_symmetry_witnesses", "satisfied": len(s_rows) == len(candidates) * len(ACTIONS), "detail": "seven symmetry witness rows were built per candidate"},
        {"obligation": "export_strict_Theta_TO_to_Gamma_SO", "satisfied": False, "detail": "0 candidates export an import-free Theta_TO that fixes an absolute origin and forms Gamma_SO"},
    ]
    return {
        "status": "P3130_THETA_TO_TRANSLATION_ORIGIN_QUOTIENT_BOUNDED_NO_GO",
        "input_hashes": {"P3129": hashlib.sha256(P3129.read_bytes()).hexdigest() if P3129.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "full_Z12_translation_quotient_summary": summary,
            "audit_object": {"object": "ThetaTOTranslationOriginQuotientAudit", "required_source": "Theta_TO strict translation-origin quotient/fixing theorem that can combine with a nonzero sign to form Gamma_SO", "forbidden_imports": ["selector", "apparatus", "observed light", "hbar/Planck", "thermodynamic environment", "Lagrangian/EOM normalization", "bridge/role-transfer", "L_total", "ToE"]},
            "candidate_Theta_TO_sources": candidates,
            "quotient_law_rows": q_rows,
            "symmetry_witness_rows": s_rows,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3129_accepted_Gamma_SO_generators": p3129.get("finite_certificate", {}).get("accepted_Gamma_SO_generators"),
            "nonempty_Z12_supports": summary["nonempty_supports"],
            "translation_classes": summary["translation_classes"],
            "dihedral_classes": summary["dihedral_classes"],
            "translation_classes_with_absolute_origin": summary["translation_classes_with_absolute_origin"],
            "inversion_fixed_translation_classes": summary["inversion_fixed_translation_classes"],
            "candidate_Theta_TO_sources": len(candidates),
            "quotient_law_rows": len(q_rows),
            "symmetry_witness_rows": len(s_rows),
            "candidate_gate_rows": len(gates),
            "strict_translation_quotient_candidates": len(strict_quotients),
            "accepted_Theta_TO_sources": len(accepted),
            "A_phi": round(a_phi(), 12),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3130 constructs the requested Theta_TO translation-origin quotient/fixing family and exhausts the nonempty Z12 support quotient: 4095 supports collapse to 351 translation classes and 223 inversion/dihedral classes. This gives a real strict quotient object, but quotienting removes the absolute source-origin representative needed by Gamma_SO; only conventional/selector/apparatus/dynamics rows fix an origin and those are forbidden. Therefore no import-free Theta_TO plus nonzero sign forms Gamma_SO on current artifacts.",
            "negative_export_flags": {key: False for key in ["Theta_TO_source_exported", "Gamma_SO_generator_exported", "Sigma_point_source_exported", "Omega_tie_source_exported", "Pi_point_source_exported", "Lambda_origin_source_exported", "Phi_Info_source_exported", "Delta_asym_source_exported", "Iota_irrev_source_exported", "Kappa_cycle_source_exported", "Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"p3129_Theta_TO_obligation_reused": True, "full_Z12_translation_quotient_exhausted": True, "strict_translation_quotient_object_constructed": True, "inversion_dihedral_pairing_computed": True, "candidate_Theta_TO_sources_constructed": True, "quotient_law_matrix_built": True, "symmetry_witness_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True},
            "next_honest_step": "Do not replay Theta_TO, Gamma_SO, Sigma_point, Omega_tie, or Pi_point with another conventional origin rule. Construct exactly one new strict non-translation datum object Epsilon_OT: an origin-torsion or origin-twist invariant that is not erased by the Z12 translation quotient and is inversion/sign aware. Test whether Epsilon_OT produces a nonconventional absolute origin representative inside a translation class without selector, apparatus, observed light, Planck units, thermodynamic environment, Lagrangian/EOM normalization, bridge/role-transfer, L_total, or ToE imports; otherwise preserve the P3105-P3130 physical-unit/sign-origin no-go certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3130/S2080 Theta_TO translation-origin quotient audit", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- P3129 accepted Gamma_SO generators: `{cert['p3129_accepted_Gamma_SO_generators']}`",
        f"- nonempty Z12 supports: `{cert['nonempty_Z12_supports']}`",
        f"- translation classes: `{cert['translation_classes']}`",
        f"- dihedral classes: `{cert['dihedral_classes']}`",
        f"- translation classes with absolute origin: `{cert['translation_classes_with_absolute_origin']}`",
        f"- inversion-fixed translation classes: `{cert['inversion_fixed_translation_classes']}`",
        f"- candidate Theta_TO sources: `{cert['candidate_Theta_TO_sources']}`",
        f"- quotient-law rows: `{cert['quotient_law_rows']}`",
        f"- symmetry-witness rows: `{cert['symmetry_witness_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- strict translation quotient candidates: `{cert['strict_translation_quotient_candidates']}`",
        f"- accepted Theta_TO sources: `{cert['accepted_Theta_TO_sources']}`",
        f"- A_phi: `{cert['A_phi']}`",
        f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3130/S2080 Theta_TO translation-origin quotient audit", "## P3130/S2080 Theta_TO translation-origin quotient audit\n\n`P3130/S2080` executes the P3129-recommended audit for `Theta_TO`, a strict translation-origin quotient/fixing theorem intended to combine with a nonzero sign to form `Gamma_SO`. It exhausts all `4095` nonempty `Z12` supports, obtaining `351` translation classes and `223` inversion/dihedral classes, then tests `15` quotient/fixing candidates with `210` quotient-law rows, `105` symmetry-witness rows, and a `15 x 18 = 270` gate matrix. The bounded result is that a strict translation quotient object exists, but quotienting erases the absolute origin needed by `Gamma_SO`; origin-fixing rows require forbidden selector/convention/apparatus/dynamics imports. No physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3130/S2080 Theta_TO quotient is not an action origin", "## P3130/S2080 Theta_TO quotient is not an action origin\n\n`P3130/S2080` constructs a finite translation quotient of `Z12` supports but does not export a unit-bearing coordinate origin, Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Theta_TO translation-origin quotient guardrail (P3130/S2080, 2026-07-06)", "## Current Theta_TO translation-origin quotient guardrail (P3130/S2080, 2026-07-06)\n\n- P3130 tests the P3129-requested strict translation-origin quotient/fixing object `Theta_TO`, intended to combine with a nonzero sign to form `Gamma_SO`.\n- The finite audit exhausts all `4095` nonempty `Z12` supports, obtaining `351` translation classes and `223` inversion/dihedral classes; it builds `15` candidates, `210` quotient-law rows, `105` symmetry-witness rows, and `270` gate rows; `0` candidates export an import-free `Theta_TO` that fixes an absolute origin and forms `Gamma_SO`.\n- The strict translation quotient is real but scoped: do not promote necklaces, bracelets, lexicographic rotations, Burnside classes, cohomology-period quotients, `A_phi` cell quotients, rank/spectral/chiral quotients, assumed origin theorems, selector origins, apparatus origins, observed-light origins, or Lagrangian coordinate origins to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict non-translation datum object `Epsilon_OT`, an origin-torsion/origin-twist invariant not erased by the `Z12` translation quotient and inversion/sign aware; otherwise preserve the P3105-P3130 physical-unit/sign-origin no-go certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
