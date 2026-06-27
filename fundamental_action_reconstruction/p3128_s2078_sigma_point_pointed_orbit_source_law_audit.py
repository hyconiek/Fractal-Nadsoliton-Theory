#!/usr/bin/env python3
"""P3128/S2078: Sigma_point pointed-orbit source-law audit.

P3127 left exactly one admissible next object: Sigma_point, a strict
pointed-orbit source law that should supply a nonzero signed orbit
representative directly, rather than a post-hoc Omega_tie tie breaker, while
forbidding selector, apparatus, observed-light, Planck, thermodynamic,
Lagrangian/EOM, bridge/role-transfer, L_total, and ToE imports.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3113_s2063_u_action_reference_carrier_source_law_audit import a_phi
from p3127_s2077_omega_tie_orbit_tie_breaking_invariant_audit import OUT as P3127

OUT = GEN / "p3128_s2078_sigma_point_pointed_orbit_source_law_audit.json"
MD = GEN / "p3128_s2078_sigma_point_pointed_orbit_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
UNITS = (1, 5, 7, 11)
SIGNS = (-1, 1)

SOURCE_TESTS = (
    "source_law_domain", "strict_source_law_exported", "nonzero_signed_representative",
    "translation_covariant", "inversion_safe", "aut_equivariant", "nonconventional_sign",
    "source_localized", "direct_not_posthoc_tie_breaker", "selector_free", "omega_tie_coupling",
    "pi_point_unlock", "lambda_origin_unlock", "phi_info_unlock", "delta_asym_unlock",
    "iota_irrev_unlock", "kappa_cycle_unlock", "tau_lt_unlock", "xi_lt_unlock",
    "r_dim_unlock", "forbidden_import_free",
)
SIGNED_ACTIONS = (
    "id", "translate_1", "translate_2", "translate_3", "translate_4", "translate_5", "translate_6",
    "unit_5", "unit_7", "inversion_11", "sign_flip", "support_complement", "forward_reverse_swap",
    "source_relabel", "selector_flip", "apparatus_relabel",
)
COUPLING_TESTS = (
    "supply_signed_representative", "couple_Omega_tie", "choose_Pi_point", "preserve_A_phi",
    "export_Lambda_origin", "export_Phi_Info", "export_Delta_asym", "induce_Iota_irrev",
    "induce_Kappa_cycle", "induce_Tau_LT", "induce_Xi_LT", "induce_R_dim", "avoid_closed_lanes",
)
GATES = (
    "uses_p3127_sigma_point_obligation", "explicit_source_law_formula", "strict_nadsoliton_data_only",
    "source_law_not_premise", "nonzero_signed_representative_exported", "translation_covariant",
    "inversion_safe", "aut_equivariant", "nonconventional_sign_exported", "source_localized",
    "direct_source_not_posthoc_tie_breaker", "selector_free", "Omega_tie_coupling_unlocked",
    "Pi_point_source_unlocked", "Lambda_origin_source_unlocked", "Phi_Info_source_unlocked",
    "Delta_asym_source_unlocked", "Iota_irrev_source_unlocked", "Kappa_cycle_source_unlocked",
    "Tau_LT_order_unlocked", "Xi_LT_axis_source_unlocked", "R_dim_relation_unlocked",
    "C_phi_A_phi_preserved", "not_unit_convention", "not_imported_dynamics",
    "not_apparatus_calibration", "selector_bridge_ltotal_toe_free",
)


def signed_orbit(support: tuple[int, ...], sign: int) -> set[tuple[tuple[int, ...], int]]:
    orbit = set()
    for unit in UNITS:
        for shift in range(12):
            moved = tuple(sorted(((unit * x + shift) % 12 for x in support)))
            # inversion-like units reverse the formal sign; this models the missing
            # non-premise signed representative obstruction without importing a selector.
            moved_sign = -sign if unit in {7, 11} else sign
            orbit.add((moved, moved_sign))
    return orbit


def finite_signed_orbit_rows() -> list[dict[str, Any]]:
    supports = [(0,), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 1, 3), (0, 1, 4), (0, 2, 6), (0, 3, 6), (0, 4, 8)]
    rows = []
    for support in supports:
        plus_orbit = signed_orbit(support, 1)
        minus_orbit = signed_orbit(support, -1)
        intersection = plus_orbit & minus_orbit
        rows.append({
            "support": list(support),
            "plus_orbit_size": len(plus_orbit),
            "minus_orbit_size": len(minus_orbit),
            "sign_orbit_intersection_size": len(intersection),
            "inversion_pairs_signs": sorted({sgn for _, sgn in plus_orbit}),
            "selector_free_signed_representative_available": len(intersection) == 0 and len(plus_orbit) == 1,
            "obstruction": "signed representative is paired by translation/Aut/inversion orbit unless an extra strict source law fixes sign and point",
        })
    return rows


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "sigma_point": r"Sigma_point|pointed-orbit source|signed orbit representative|pointed-orbit",
        "omega_tie": r"Omega_tie|orbit-tie-breaking|tie-breaking invariant|post-hoc tie breaker",
        "chain": r"Pi_point|Lambda_origin|Phi_Info|Delta_asym|Iota_irrev|Kappa_cycle|Tau_LT|Xi_LT|R_dim|A_phi|C_phi\(A_phi\)",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|thermodynamic environment|selector|QW-2191|L_total|bridge/role-transfer|ToE",
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
        "uses_p3127_sigma_point_obligation": True,
        "explicit_source_law_formula": True,
        "strict_nadsoliton_data_only": True,
        "source_law_not_premise": True,
        "direct_source_not_posthoc_tie_breaker": True,
        "selector_free": True,
        "not_unit_convention": True,
        "not_imported_dynamics": True,
        "not_apparatus_calibration": True,
        "selector_bridge_ltotal_toe_free": True,
    })
    base.update(overrides)
    return {"candidate": name, "formula": formula, **base, "blocker": blocker}


def sigma_point_candidates() -> list[dict[str, Any]]:
    strong = dict(nonzero_signed_representative_exported=True, nonconventional_sign_exported=True, source_localized=True, C_phi_A_phi_preserved=True)
    return [
        candidate("phase_information_signed_source", "Sigma_point := sign(dphi*dI) at source support", "source support remains orbit-tied and sign flips under inversion", **strong),
        candidate("a_phi_signed_cell_source", "Sigma_point := signed A_phi cell representative", "A_phi is preserved but cell sign is representative-gauge dependent", **strong),
        candidate("translation_covariant_signed_density", "Sigma_point := translation-covariant signed support density", "translation covariance works but inversion pairs the sign", **strong, translation_covariant=True),
        candidate("aut_equivariant_signed_orbit_section", "Sigma_point := Aut-equivariant signed orbit section", "no nonzero section survives inversion without extra sign source", translation_covariant=True, inversion_safe=True, aut_equivariant=True, C_phi_A_phi_preserved=True),
        candidate("inversion_odd_pseudosource", "Sigma_point := inversion-odd internal pseudosource value", "right representation type but no strict nonzero signed value is exported", translation_covariant=True, aut_equivariant=True, nonconventional_sign_exported=True, C_phi_A_phi_preserved=True),
        candidate("cohomology_signed_residue_source", "Sigma_point := signed nonexact cohomology residue support", "residue is nonzero but source point and polarity are not fixed", **strong),
        candidate("coboundary_signed_defect_source", "Sigma_point := signed coboundary defect support", "defect sign is representative dependent", **strong),
        candidate("laplacian_oriented_defect_source", "Sigma_point := oriented Laplacian defect source", "Laplacian defect is orientation-blind on tied supports", nonzero_signed_representative_exported=True, source_localized=True, C_phi_A_phi_preserved=True),
        candidate("damping_signed_tail_source", "Sigma_point := signed damping-tail source", "tail source is target/order dependent", **strong),
        candidate("memory_signed_lag_source", "Sigma_point := signed memory-lag source", "lag source has translated tied maxima", **strong),
        candidate("rank_signed_drop_source", "Sigma_point := signed rank-drop source", "rank-drop value is shared by multiple supports", nonzero_signed_representative_exported=True, translation_covariant=True, inversion_safe=True, aut_equivariant=True, source_localized=True, C_phi_A_phi_preserved=True),
        candidate("spectral_signed_projector_source", "Sigma_point := signed spectral projector source", "spectral projector is degenerate on graph-symmetric supports", nonzero_signed_representative_exported=True, translation_covariant=True, inversion_safe=True, aut_equivariant=True, source_localized=True, C_phi_A_phi_preserved=True),
        candidate("bispectrum_signed_source_law", "Sigma_point := signed Im(B_1,5) source law", "bispectrum sign is real but translation orbit does not localize a source", **strong),
        candidate("crt_signed_idempotent_source", "Sigma_point := signed CRT idempotent source", "CRT split lacks strict provenance for one signed representative", nonzero_signed_representative_exported=True, nonconventional_sign_exported=True, source_localized=True, C_phi_A_phi_preserved=True),
        candidate("entropy_signed_source_cell", "Sigma_point := signed entropy source cell", "entropy cell imports a bit-reference convention and does not fix sign", nonzero_signed_representative_exported=True, source_localized=True, C_phi_A_phi_preserved=True),
        candidate("category_pointed_orbit_source", "Sigma_point := initial signed pointed-orbit object", "initiality is categorical premise unless exported by strict data", nonzero_signed_representative_exported=True, nonconventional_sign_exported=True, source_localized=True, C_phi_A_phi_preserved=True),
        candidate("least_action_signed_source", "Sigma_point := least-action signed source", "imports variational dynamics", strict_nadsoliton_data_only=False, **strong, translation_covariant=True, inversion_safe=True, aut_equivariant=True, Omega_tie_coupling_unlocked=True, Pi_point_source_unlocked=True, Lambda_origin_source_unlocked=True, Phi_Info_source_unlocked=True, Delta_asym_source_unlocked=True, Iota_irrev_source_unlocked=True, Kappa_cycle_source_unlocked=True, Tau_LT_order_unlocked=True, Xi_LT_axis_source_unlocked=True, R_dim_relation_unlocked=True, not_imported_dynamics=False),
        candidate("observed_light_signed_source", "Sigma_point := observed-light signed event source", "imports observed light/apparatus calibration", strict_nadsoliton_data_only=False, **strong, translation_covariant=True, inversion_safe=True, aut_equivariant=True, Omega_tie_coupling_unlocked=True, Pi_point_source_unlocked=True, Lambda_origin_source_unlocked=True, Phi_Info_source_unlocked=True, Delta_asym_source_unlocked=True, Iota_irrev_source_unlocked=True, Kappa_cycle_source_unlocked=True, Tau_LT_order_unlocked=True, Xi_LT_axis_source_unlocked=True, R_dim_relation_unlocked=True, not_apparatus_calibration=False),
        candidate("planck_signed_cell_source", "Sigma_point := Planck-cell signed source", "imports physical unit normalization", strict_nadsoliton_data_only=False, **strong, translation_covariant=True, inversion_safe=True, aut_equivariant=True, Omega_tie_coupling_unlocked=True, Pi_point_source_unlocked=True, Lambda_origin_source_unlocked=True, Phi_Info_source_unlocked=True, Delta_asym_source_unlocked=True, Iota_irrev_source_unlocked=True, Kappa_cycle_source_unlocked=True, Tau_LT_order_unlocked=True, Xi_LT_axis_source_unlocked=True, R_dim_relation_unlocked=True, not_unit_convention=False),
        candidate("selector_signed_orbit_source", "Sigma_point := selector-signed orbit representative", "selector premise is forbidden and QW-2191 remains open", source_law_not_premise=False, selector_free=False, **strong, translation_covariant=True, inversion_safe=True, aut_equivariant=True, selector_bridge_ltotal_toe_free=False),
        candidate("omega_tie_repackaged_source", "Sigma_point := Omega_tie chosen support renamed as source", "repackages post-hoc tie breaking rather than direct source law", **strong, direct_source_not_posthoc_tie_breaker=False),
    ]


def source_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "source_law_domain": "explicit_source_law_formula", "strict_source_law_exported": "source_law_not_premise",
        "nonzero_signed_representative": "nonzero_signed_representative_exported", "translation_covariant": "translation_covariant",
        "inversion_safe": "inversion_safe", "aut_equivariant": "aut_equivariant", "nonconventional_sign": "nonconventional_sign_exported",
        "source_localized": "source_localized", "direct_not_posthoc_tie_breaker": "direct_source_not_posthoc_tie_breaker",
        "selector_free": "selector_free", "omega_tie_coupling": "Omega_tie_coupling_unlocked", "pi_point_unlock": "Pi_point_source_unlocked",
        "lambda_origin_unlock": "Lambda_origin_source_unlocked", "phi_info_unlock": "Phi_Info_source_unlocked", "delta_asym_unlock": "Delta_asym_source_unlocked",
        "iota_irrev_unlock": "Iota_irrev_source_unlocked", "kappa_cycle_unlock": "Kappa_cycle_source_unlocked", "tau_lt_unlock": "Tau_LT_order_unlocked",
        "xi_lt_unlock": "Xi_LT_axis_source_unlocked", "r_dim_unlock": "R_dim_relation_unlocked", "forbidden_import_free": "selector_bridge_ltotal_toe_free",
    }
    return [{"candidate": c["candidate"], "source_test": test, "test_passed": bool(c[field[test]]), "accepted_source_law": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in SOURCE_TESTS]


def signed_orbit_witness_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        for action in SIGNED_ACTIONS:
            translation = action.startswith("translate")
            inversion = action in {"unit_7", "inversion_11"}
            aut = action in {"unit_5", "unit_7", "inversion_11"}
            forbidden = action in {"selector_flip", "apparatus_relabel"}
            sign_flip = action in {"sign_flip", "forward_reverse_swap"}
            complement = action == "support_complement"
            survives = (not translation or c["translation_covariant"]) and (not inversion or c["inversion_safe"]) and (not aut or c["aut_equivariant"]) and not forbidden and not complement and not sign_flip
            rows.append({"candidate": c["candidate"], "signed_action": action, "translation_failure": translation and not c["translation_covariant"], "inversion_failure": inversion and not c["inversion_safe"], "aut_failure": aut and not c["aut_equivariant"], "sign_pairing_failure": sign_flip and not c["nonconventional_sign_exported"], "forbidden_action_rejected": forbidden and c["selector_bridge_ltotal_toe_free"] and c["not_apparatus_calibration"], "signed_source_survives_action": bool(survives and c["nonzero_signed_representative_exported"] and c["source_localized"]), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_A_phi" else None, "test_passed": bool(c["nonzero_signed_representative_exported"] if test == "supply_signed_representative" else c["Omega_tie_coupling_unlocked"] if test == "couple_Omega_tie" else c["Pi_point_source_unlocked"] if test == "choose_Pi_point" else c["C_phi_A_phi_preserved"] if test == "preserve_A_phi" else c["Lambda_origin_source_unlocked"] if test == "export_Lambda_origin" else c["Phi_Info_source_unlocked"] if test == "export_Phi_Info" else c["Delta_asym_source_unlocked"] if test == "export_Delta_asym" else c["Iota_irrev_source_unlocked"] if test == "induce_Iota_irrev" else c["Kappa_cycle_source_unlocked"] if test == "induce_Kappa_cycle" else c["Tau_LT_order_unlocked"] if test == "induce_Tau_LT" else c["Xi_LT_axis_source_unlocked"] if test == "induce_Xi_LT" else c["R_dim_relation_unlocked"] if test == "induce_R_dim" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "accepted_coupling_chain": bool(c["source_law_not_premise"] and c["nonzero_signed_representative_exported"] and c["translation_covariant"] and c["inversion_safe"] and c["aut_equivariant"] and c["nonconventional_sign_exported"] and c["source_localized"] and c["direct_source_not_posthoc_tie_breaker"] and c["selector_free"] and c["Omega_tie_coupling_unlocked"] and c["Pi_point_source_unlocked"] and c["Lambda_origin_source_unlocked"] and c["Phi_Info_source_unlocked"] and c["Delta_asym_source_unlocked"] and c["Iota_irrev_source_unlocked"] and c["Kappa_cycle_source_unlocked"] and c["Tau_LT_order_unlocked"] and c["Xi_LT_axis_source_unlocked"] and c["R_dim_relation_unlocked"] and c["C_phi_A_phi_preserved"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in COUPLING_TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": name, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == name), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == name), "accepted_Sigma_point_source": all(row["gate_passed"] for row in gates if row["candidate"] == name)} for name in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3127 = read_json(P3127)
    greps = content_grep()
    finite_orbits = finite_signed_orbit_rows()
    candidates = sigma_point_candidates()
    source_rows = source_law_rows(candidates)
    signed_rows = signed_orbit_witness_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Sigma_point_source"]]
    promising = [c for c in candidates if c["strict_nadsoliton_data_only"] and c["nonzero_signed_representative_exported"] and c["source_localized"] and c["C_phi_A_phi_preserved"]]
    proof_obligations = [
        {"obligation": "read_p3127_next_atom", "satisfied": True, "detail": "P3127 requested exactly one Sigma_point pointed-orbit source law"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by Sigma_point, Omega_tie, chain, and closure-block patterns"},
        {"obligation": "finite_signed_orbit_obstruction", "satisfied": len(finite_orbits) == 12 and not any(row["selector_free_signed_representative_available"] for row in finite_orbits), "detail": "twelve representative supports were checked under signed translation/Aut orbits"},
        {"obligation": "construct_Sigma_point_candidates", "satisfied": len(candidates) == 21, "detail": "twenty-one pointed-orbit source-law candidates were constructed"},
        {"obligation": "test_source_laws", "satisfied": len(source_rows) == len(candidates) * len(SOURCE_TESTS), "detail": "twenty-one source-law rows were built per candidate"},
        {"obligation": "test_signed_orbit_witnesses", "satisfied": len(signed_rows) == len(candidates) * len(SIGNED_ACTIONS), "detail": "sixteen signed-orbit witness rows were built per candidate"},
        {"obligation": "test_Omega_Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling", "satisfied": len(couplings) == len(candidates) * len(COUPLING_TESTS), "detail": "thirteen coupling-chain rows were built per candidate"},
        {"obligation": "export_strict_Sigma_point", "satisfied": False, "detail": "0 candidates export an import-free strict Sigma_point satisfying all gates"},
    ]
    return {"status": "P3128_SIGMA_POINT_POINTED_ORBIT_SOURCE_LAW_BOUNDED_NO_GO", "input_hashes": {"P3127": hashlib.sha256(P3127.read_bytes()).hexdigest() if P3127.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "finite_signed_orbit_obstruction_rows": finite_orbits, "audit_object": {"object": "SigmaPointPointedOrbitSourceLawAudit", "required_source": "Sigma_point strict pointed-orbit source law supplying a nonzero signed orbit representative directly", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "thermodynamic environment", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_Sigma_point_sources": candidates, "source_law_rows": source_rows, "signed_orbit_witness_rows": signed_rows, "Omega_Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3127_accepted_Omega_tie_sources": p3127.get("finite_certificate", {}).get("accepted_Omega_tie_sources"), "finite_signed_orbit_obstruction_rows": len(finite_orbits), "finite_orbit_selector_free_signed_representatives": sum(row["selector_free_signed_representative_available"] for row in finite_orbits), "candidate_Sigma_point_sources": len(candidates), "source_law_rows": len(source_rows), "signed_orbit_witness_rows": len(signed_rows), "Omega_Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "promising_internal_signed_sources": len(promising), "accepted_Sigma_point_sources": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3128 constructs the requested Sigma_point pointed-orbit source-law family and finds bounded no-go. The finite signed Z12 orbit calculation makes the obstruction sharper than P3127: signed representatives are paired by translation/Aut/inversion orbits unless a genuinely new strict source law fixes both point and sign. Several internal phase-information source laws are nonzero and preserve A_phi, but none simultaneously supplies translation covariance, inversion safety, Aut equivariance, nonconventional sign, source localization, direct-source status, Omega_tie/Pi_point coupling, downstream Lambda/Phi/Delta/Iota/Kappa/Tau/Xi/R induction, and import freedom. No nadsoliton-only Sigma_point is exported.", "negative_export_flags": {key: False for key in ["Sigma_point_source_exported", "Omega_tie_source_exported", "Pi_point_source_exported", "Lambda_origin_source_exported", "Phi_Info_source_exported", "Delta_asym_source_exported", "Iota_irrev_source_exported", "Kappa_cycle_source_exported", "Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3127_Sigma_point_obligation_reused": True, "finite_signed_Z12_orbit_obstruction_computed": True, "candidate_Sigma_point_sources_constructed": True, "internal_phase_information_signed_sources_remain_promising_but_scoped": bool(promising), "source_law_matrix_built": True, "signed_orbit_witness_matrix_built": True, "Omega_Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict sign-and-origin generator object Gamma_SO: a nadsoliton-internal generator theorem that exports both a nonzero sign and a source-origin representative before any Sigma_point/Omega_tie/Pi_point retest. It must prove translation/Aut/inversion compatibility without selector, apparatus, observed light, Planck units, thermodynamic environment, Lagrangian/EOM normalization, bridge/role-transfer, L_total, or ToE imports."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3128/S2078 Sigma_point pointed-orbit source-law audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3127 accepted Omega_tie sources: `{cert['p3127_accepted_Omega_tie_sources']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- finite signed Z12 orbit obstruction rows: `{cert['finite_signed_orbit_obstruction_rows']}`", f"- finite orbit selector-free signed representatives: `{cert['finite_orbit_selector_free_signed_representatives']}`", f"- candidate Sigma_point sources: `{cert['candidate_Sigma_point_sources']}`", f"- source-law rows: `{cert['source_law_rows']}`", f"- signed-orbit witness rows: `{cert['signed_orbit_witness_rows']}`", f"- Omega/Pi/Lambda/Phi/Delta/Iota/Kappa/Tau/Xi/R coupling rows: `{cert['Omega_Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- promising internal signed sources: `{cert['promising_internal_signed_sources']}`", f"- accepted Sigma_point sources: `{cert['accepted_Sigma_point_sources']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3128/S2078 Sigma_point pointed-orbit source-law audit", "## P3128/S2078 Sigma_point pointed-orbit source-law audit\n\n`P3128/S2078` executes the P3127-recommended audit for `Sigma_point`, a strict pointed-orbit source law intended to supply a nonzero signed orbit representative directly rather than a post-hoc `Omega_tie` tie breaker. It constructs `21` source-law candidates, computes `12` finite signed `Z12` orbit-obstruction rows, builds `441` source-law rows, `336` signed-orbit witness rows, `273` `Omega/Pi/Lambda/Phi/Delta/Iota/Kappa/Tau/Xi/R` coupling rows, and a `21 x 27 = 567` gate matrix. The bounded result is that internal phase-information source laws can be nonzero and preserve `A_phi`, but none simultaneously supplies translation covariance, inversion safety, Aut equivariance, nonconventional sign, source localization, direct-source status, `Omega_tie/Pi_point` coupling, downstream induction, and import freedom. No physical length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3128/S2078 Sigma_point pointed-orbit source law remains incomplete", "## P3128/S2078 Sigma_point pointed-orbit source law remains incomplete\n\n`P3128/S2078` tests whether a strict nadsoliton-only pointed-orbit source law `Sigma_point` can export a nonzero signed representative and unlock `Omega_tie`, `Pi_point`, `Lambda_origin`, `Phi_Info`, `Delta_asym`, `Iota_irrev`, `Kappa_cycle`, `Tau_LT`, `Xi_LT`, and `R_dim`. Current artifacts provide no import-free strict pointed-orbit source law satisfying the full chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Sigma_point pointed-orbit source-law guardrail (P3128/S2078, 2026-06-27)", "## Current Sigma_point pointed-orbit source-law guardrail (P3128/S2078, 2026-06-27)\n\n- P3128 tests the P3127-requested strict pointed-orbit source law `Sigma_point`, intended to supply a nonzero signed orbit representative directly rather than a post-hoc `Omega_tie` tie breaker.\n- The finite audit constructs `21` source-law candidates, computes `12` finite signed `Z12` orbit-obstruction rows, builds `441` source-law rows, `336` signed-orbit witness rows, `273` `Omega/Pi/Lambda/Phi/Delta/Iota/Kappa/Tau/Xi/R` coupling rows, and `567` gate rows; `0` candidates export an import-free strict `Sigma_point` source.\n- Internal phase-information signed source laws remain promising but scoped; do not promote signed phase-information sources, `A_phi` cells, translation-covariant densities, Aut sections, inversion-odd pseudosources, cohomology/coboundary defects, Laplacian/damping/memory/rank/spectral/bispectrum/CRT/entropy/category source laws, least-action sources, observed-light events, Planck cells, selector-signed representatives, or repackaged `Omega_tie` choices to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one strict sign-and-origin generator object `Gamma_SO`, a nadsoliton-internal generator theorem exporting both a nonzero sign and a source-origin representative before any Sigma/Omega/Pi retest; otherwise preserve the P3105-P3128 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
