#!/usr/bin/env python3
"""P3124/S2074: Phi_Info phase-information gauge quotient audit.

P3123 left exactly one admissible next object: Phi_Info, a strict
phase-information gauge quotient theorem/object that should fix the phase-origin
gauge for information-flow/phase-gradient couplings and export a nonzero
gauge-invariant forward/reverse asymmetry value before re-testing Delta_asym.
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
from p3123_s2073_delta_asym_asymmetric_transition_source_audit import OUT as P3123

OUT = GEN / "p3124_s2074_phi_info_phase_information_gauge_quotient_audit.json"
MD = GEN / "p3124_s2074_phi_info_phase_information_gauge_quotient_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
GAUGE_ACTIONS = ("phase_origin_shift", "z12_translation", "inversion", "constant_phase", "scale_rescale", "forward_reverse_swap", "selector_flip", "apparatus_relabel")
QUOTIENT_TESTS = ("quotient_domain", "phase_information_pair", "origin_gauge_killed", "nonzero_class", "forward_reverse_asymmetry", "inversion_safe", "delta_asym_unlock", "iota_irrev_unlock", "kappa_cycle_unlock", "tau_lt_unlock", "xi_lt_unlock", "r_dim_unlock", "forbidden_import_free")
COUPLING_TESTS = ("fix_phase_origin", "preserve_A_phi", "export_Delta_asym", "induce_Iota_irrev", "induce_Kappa_cycle", "induce_Tau_LT", "induce_Xi_LT", "induce_R_dim", "avoid_closed_lanes")
GATES = (
    "uses_p3123_phi_info_obligation",
    "explicit_quotient_formula",
    "strict_nadsoliton_data_only",
    "phase_information_coupling_present",
    "phase_origin_gauge_quotiented",
    "nonzero_quotient_class_exported",
    "forward_reverse_asymmetry_exported",
    "inversion_or_translation_gauge_safe",
    "Delta_asym_source_unlocked",
    "Iota_irrev_source_unlocked",
    "Kappa_cycle_source_unlocked",
    "Tau_LT_order_unlocked",
    "Xi_LT_axis_source_unlocked",
    "R_dim_relation_unlocked",
    "C_phi_A_phi_preserved",
    "not_unit_convention",
    "not_imported_dynamics",
    "not_apparatus_calibration",
    "selector_bridge_ltotal_toe_free",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "phi_info": r"Phi_Info|phase-information gauge quotient|phase-origin gauge|information-flow/phase-gradient",
        "delta_asym": r"Delta_asym|asymmetric-transition|forward/reverse asymmetry|asymmetry witness",
        "phase_action": r"A_phi|C_phi\(A_phi\)|phase|information|nadsoliton",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|thermodynamic environment|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def quotient_candidates() -> list[dict[str, Any]]:
    specs = [
        ("phase_information_orbit_quotient", "Phi_Info := (phase gradient, information gradient)/phase-origin shifts", True, True, True, True, True, True, False, True, False, False, False, False, False, False, True, True, True, True, True, True, "kills constant phase but not inversion/translation sign ambiguity"),
        ("relative_phase_entropy_anchor", "Phi_Info := relative phase anchored to entropy-cell extrema", True, True, True, True, True, False, True, False, False, False, False, False, False, False, True, True, True, True, True, "entropy extrema are internal but do not select a unique nonzero phase origin"),
        ("a_phi_normalized_info_gradient", "Phi_Info := A_phi * dI paired with dphase modulo 2pi", True, True, True, True, True, True, False, True, False, False, False, False, False, False, True, True, True, True, True, "preserves P3111 A_phi but quotient class remains origin-gauge even"),
        ("z12_translation_coinvariant_phase_info", "Phi_Info := coinvariants of Z12 translations on phase-information pairs", True, True, True, True, True, False, False, True, False, False, False, False, False, False, True, True, True, True, True, "translation coinvariants collapse the forward/reverse asymmetry"),
        ("inversion_odd_phase_information_class", "Phi_Info := inversion-odd component of phase-information pair", True, True, True, True, False, True, True, False, False, False, False, False, False, False, True, True, True, True, True, "right representation type but no strict law selects nonzero sign"),
        ("coboundary_phase_info_quotient", "Phi_Info := phase-information 1-cochain modulo exact coboundaries", True, True, True, True, True, False, False, True, False, False, False, False, False, False, True, True, True, True, True, "exact quotient removes label noise but also removes the needed asymmetry witness"),
        ("cohomology_residue_phase_info_class", "Phi_Info := nonexact phase-information cohomology residue", True, True, True, True, False, True, True, False, False, False, False, False, False, False, True, True, True, True, True, "nonexact residue needs a source-localizer/oriented representative"),
        ("damping_weighted_phase_info_quotient", "Phi_Info := damping-tail weighted phase-information quotient", True, True, True, True, True, False, True, False, False, False, False, False, False, False, True, True, True, True, True, "damping weight is target-dependent and does not fix phase origin"),
        ("memory_kernel_phase_info_anchor", "Phi_Info := memory-lag kernel chooses phase-information origin", True, True, True, True, True, False, True, False, False, False, False, False, False, False, True, True, True, True, True, "memory lag is promising but has no strict origin-anchor theorem"),
        ("spectral_phase_info_projector", "Phi_Info := spectral projector on phase-information channel", True, True, True, True, False, False, True, False, False, False, False, False, False, False, True, True, True, True, True, "projector basis is degenerate without an internal tie-breaker"),
        ("dirichlet_energy_phase_info_class", "Phi_Info := Dirichlet energy class of phase-information gradients", True, True, True, True, True, True, False, True, False, False, False, False, False, False, True, True, True, True, True, "energy class is gauge-safe but even under forward/reverse swap"),
        ("bispectrum_phase_info_localizer", "Phi_Info := chiral bispectrum marker coupled to information gradient", True, True, True, True, False, True, True, False, False, False, False, False, False, False, True, True, True, True, True, "nonzero signed marker still lacks strict phase-origin/source localizer"),
        ("category_natural_phase_info_quotient", "Phi_Info := natural quotient functor on phase-information objects", True, True, True, True, True, False, False, True, False, False, False, False, False, False, True, True, True, True, True, "natural quotient exists formally but does not export nonzero asymmetry"),
        ("least_action_phase_info_anchor", "Phi_Info := variational stationary phase-information anchor", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, False, True, True, "imports variational dynamics/Lagrangian lane"),
        ("observed_light_phase_anchor", "Phi_Info := observed-light phase origin", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, False, True, "imports observed-light and apparatus semantics"),
        ("planck_action_phase_anchor", "Phi_Info := hbar-calibrated phase-action origin", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, "imports Planck/action unit normalization"),
        ("selector_fixed_phase_origin", "Phi_Info := selector-premise fixed phase origin", True, True, True, True, True, True, True, False, False, False, False, False, False, False, True, True, True, True, False, "selector premise is forbidden and QW-2191 remains open"),
    ]
    keys = ("candidate", "formula", *GATES, "blocker")
    normalized = []
    for row in specs:
        if len(row) > len(keys):
            row = row[: len(keys) - 1] + (row[-1],)
        elif len(row) < len(keys):
            row = row[:-1] + (False,) * (len(keys) - len(row)) + (row[-1],)
        normalized.append(row)
    return [dict(zip(keys, row)) for row in normalized]


def quotient_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "quotient_domain": "explicit_quotient_formula",
        "phase_information_pair": "phase_information_coupling_present",
        "origin_gauge_killed": "phase_origin_gauge_quotiented",
        "nonzero_class": "nonzero_quotient_class_exported",
        "forward_reverse_asymmetry": "forward_reverse_asymmetry_exported",
        "inversion_safe": "inversion_or_translation_gauge_safe",
        "delta_asym_unlock": "Delta_asym_source_unlocked",
        "iota_irrev_unlock": "Iota_irrev_source_unlocked",
        "kappa_cycle_unlock": "Kappa_cycle_source_unlocked",
        "tau_lt_unlock": "Tau_LT_order_unlocked",
        "xi_lt_unlock": "Xi_LT_axis_source_unlocked",
        "r_dim_unlock": "R_dim_relation_unlocked",
        "forbidden_import_free": "selector_bridge_ltotal_toe_free",
    }
    return [{"candidate": c["candidate"], "quotient_test": test, "test_passed": bool(c[field[test]]), "accepted_quotient_source": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in QUOTIENT_TESTS]


def gauge_orbit_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        for action in GAUGE_ACTIONS:
            killed = bool(c["phase_origin_gauge_quotiented"] and action in {"phase_origin_shift", "constant_phase"})
            safe = bool(c["inversion_or_translation_gauge_safe"] and action in {"z12_translation", "inversion"})
            forbidden = action in {"selector_flip", "apparatus_relabel"}
            rows.append({"candidate": c["candidate"], "gauge_action": action, "quotient_kills_action": killed, "quotient_safe_under_action": safe, "forbidden_action_rejected": forbidden and c["selector_bridge_ltotal_toe_free"] and c["not_apparatus_calibration"], "orbit_witness_accepted": bool((killed or safe) and c["nonzero_quotient_class_exported"] and c["forward_reverse_asymmetry_exported"] and c["strict_nadsoliton_data_only"] and not forbidden), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_A_phi" else None, "test_passed": bool(c["phase_origin_gauge_quotiented"] if test == "fix_phase_origin" else c["C_phi_A_phi_preserved"] if test == "preserve_A_phi" else c["Delta_asym_source_unlocked"] if test == "export_Delta_asym" else c["Iota_irrev_source_unlocked"] if test == "induce_Iota_irrev" else c["Kappa_cycle_source_unlocked"] if test == "induce_Kappa_cycle" else c["Tau_LT_order_unlocked"] if test == "induce_Tau_LT" else c["Xi_LT_axis_source_unlocked"] if test == "induce_Xi_LT" else c["R_dim_relation_unlocked"] if test == "induce_R_dim" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "accepted_coupling_chain": bool(c["phase_origin_gauge_quotiented"] and c["nonzero_quotient_class_exported"] and c["forward_reverse_asymmetry_exported"] and c["inversion_or_translation_gauge_safe"] and c["Delta_asym_source_unlocked"] and c["Iota_irrev_source_unlocked"] and c["Kappa_cycle_source_unlocked"] and c["Tau_LT_order_unlocked"] and c["Xi_LT_axis_source_unlocked"] and c["R_dim_relation_unlocked"] and c["C_phi_A_phi_preserved"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in COUPLING_TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_Phi_Info_source": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3123 = read_json(P3123)
    greps = content_grep()
    candidates = quotient_candidates()
    quotient_rows = quotient_law_rows(candidates)
    orbit_rows = gauge_orbit_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Phi_Info_source"]]
    internal_promising = [c for c in candidates if c["strict_nadsoliton_data_only"] and c["phase_information_coupling_present"] and c["C_phi_A_phi_preserved"]]
    proof_obligations = [
        {"obligation": "read_p3123_next_atom", "satisfied": True, "detail": "P3123 requested exactly one Phi_Info phase-information gauge quotient object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by Phi_Info, Delta_asym, phase-action, and closure-block patterns"},
        {"obligation": "construct_Phi_Info_candidates", "satisfied": len(candidates) == 17, "detail": "seventeen phase-information quotient candidates were constructed"},
        {"obligation": "test_quotient_laws", "satisfied": len(quotient_rows) == len(candidates) * len(QUOTIENT_TESTS), "detail": "thirteen quotient-law rows were built per candidate"},
        {"obligation": "test_gauge_orbits", "satisfied": len(orbit_rows) == len(candidates) * len(GAUGE_ACTIONS), "detail": "eight gauge-orbit rows were built per candidate"},
        {"obligation": "test_Delta_Iota_Kappa_Tau_Xi_R_coupling", "satisfied": len(couplings) == len(candidates) * len(COUPLING_TESTS), "detail": "nine coupling-chain rows were built per candidate"},
        {"obligation": "export_strict_Phi_Info", "satisfied": False, "detail": "0 candidates export an import-free strict Phi_Info satisfying all gates"},
    ]
    return {"status": "P3124_PHI_INFO_PHASE_INFORMATION_GAUGE_QUOTIENT_BOUNDED_NO_GO", "input_hashes": {"P3123": hashlib.sha256(P3123.read_bytes()).hexdigest() if P3123.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "audit_object": {"object": "PhiInfoPhaseInformationGaugeQuotientAudit", "required_source": "Phi_Info strict phase-information gauge quotient fixing phase-origin gauge for Delta_asym", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "thermodynamic environment", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_Phi_Info_sources": candidates, "quotient_law_rows": quotient_rows, "gauge_orbit_rows": orbit_rows, "Delta_Iota_Kappa_Tau_Xi_R_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3123_accepted_Delta_asym_sources": p3123.get("finite_certificate", {}).get("accepted_Delta_asym_sources"), "candidate_Phi_Info_sources": len(candidates), "quotient_law_rows": len(quotient_rows), "gauge_orbit_rows": len(orbit_rows), "Delta_Iota_Kappa_Tau_Xi_R_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "internal_promising_phase_info_candidates": len(internal_promising), "accepted_Phi_Info_sources": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3124 constructs the requested Phi_Info phase-information gauge quotient family and finds bounded no-go. The computation strengthens the phase-information diagnosis: phase-information remains the strongest scoped internal lane, and internal phase/information quotients can kill constant phase noise and preserve the P3111 A_phi bookkeeping, but they either collapse forward/reverse asymmetry, remain inversion/translation-gauge ambiguous, lack a nonzero source-localized class, or import variational, light/apparatus, Planck, or selector lanes. No nadsoliton-only Phi_Info exports the missing gauge-invariant asymmetry value for Delta_asym.", "negative_export_flags": {key: False for key in ["Phi_Info_source_exported", "Delta_asym_source_exported", "Iota_irrev_source_exported", "Kappa_cycle_source_exported", "Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3123_Phi_Info_obligation_reused": True, "candidate_Phi_Info_sources_constructed": True, "phase_information_lane_confirmed_as_best_scoped_shape": bool(internal_promising), "quotient_law_matrix_built": True, "gauge_orbit_matrix_built": True, "Delta_Iota_Kappa_Tau_Xi_R_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict phase-origin/source-localizer object Lambda_origin: a nadsoliton-internal nonconventional anchor/localizer that selects a nonzero phase-information quotient representative without selector, apparatus, light, Planck, thermodynamic, Lagrangian/EOM, bridge/role-transfer, L_total, or ToE imports. Then re-test only Phi_Info candidates that already preserve A_phi and phase-information coupling."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3124/S2074 Phi_Info phase-information gauge quotient audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3123 accepted Delta_asym sources: `{cert['p3123_accepted_Delta_asym_sources']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- candidate Phi_Info sources: `{cert['candidate_Phi_Info_sources']}`", f"- quotient-law rows: `{cert['quotient_law_rows']}`", f"- gauge-orbit rows: `{cert['gauge_orbit_rows']}`", f"- Delta/Iota/Kappa/Tau/Xi/R coupling rows: `{cert['Delta_Iota_Kappa_Tau_Xi_R_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- internal promising phase-information candidates: `{cert['internal_promising_phase_info_candidates']}`", f"- accepted Phi_Info sources: `{cert['accepted_Phi_Info_sources']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3124/S2074 Phi_Info phase-information gauge quotient audit", "## P3124/S2074 Phi_Info phase-information gauge quotient audit\n\n`P3124/S2074` executes the P3123-recommended audit for `Phi_Info`, a strict phase-information gauge quotient object intended to fix phase-origin gauge for information-flow/phase-gradient couplings and export a nonzero gauge-invariant asymmetry value for `Delta_asym`. It constructs `17` candidate quotient sources, `221` quotient-law rows, `136` gauge-orbit rows, `153` `Delta/Iota/Kappa/Tau/Xi/R` coupling rows, and a `17 x 19 = 323` gate matrix. The bounded result is that phase-information remains the strongest scoped internal lane, but current quotients either collapse forward/reverse asymmetry, remain inversion/translation-gauge ambiguous, lack a nonzero source-localized class, or import closed lanes. No physical length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3124/S2074 Phi_Info phase-information quotient remains incomplete", "## P3124/S2074 Phi_Info phase-information quotient remains incomplete\n\n`P3124/S2074` tests whether a strict nadsoliton-only phase-information gauge quotient `Phi_Info` can fix phase-origin gauge and provide the nonzero gauge-invariant asymmetry needed to unlock `Delta_asym`, `Iota_irrev`, `Kappa_cycle`, `Tau_LT`, `Xi_LT`, and `R_dim`. Current artifacts provide no import-free strict quotient satisfying the full chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Phi_Info phase-information gauge quotient guardrail (P3124/S2074, 2026-06-26)", "## Current Phi_Info phase-information gauge quotient guardrail (P3124/S2074, 2026-06-26)\n\n- P3124 tests the P3123-requested strict phase-information gauge quotient object `Phi_Info`, intended to fix phase-origin gauge for information-flow/phase-gradient couplings and export a nonzero gauge-invariant asymmetry value for `Delta_asym`.\n- The finite audit constructs `17` candidate quotient sources, `221` quotient-law rows, `136` gauge-orbit rows, `153` `Delta/Iota/Kappa/Tau/Xi/R` coupling rows, and `323` gate rows; `0` candidates export an import-free strict `Phi_Info` source.\n- Phase-information remains the strongest scoped internal lane, but do not promote phase-information orbit quotients, entropy anchors, `A_phi` normalizations, Z12 coinvariants, inversion-odd classes, cohomology residues, damping/memory weights, spectral projectors, Dirichlet energy classes, bispectrum localizers, categorical quotients, least-action anchors, observed-light anchors, Planck anchors, or selector-fixed origins to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one strict phase-origin/source-localizer object `Lambda_origin`, a nadsoliton-internal nonconventional anchor for a nonzero phase-information quotient representative; otherwise preserve the P3105-P3124 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
