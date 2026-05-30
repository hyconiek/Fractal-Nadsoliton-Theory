#!/usr/bin/env python3
"""Scratch probe: legacy-kernel reactivation candidate from DIAGRAMS evidence.

This probe addresses the author's current instruction that the legacy kernel
should not be treated as dead if DIAGRAMS_KERNEL_TRANSFORMATION.md already
records the nadsoliton-characteristic content of K_total -> K(d), and if the
remaining bridge gap is now concentrated around the one-bit chi_11/unit-axis
source.

The output deliberately does *not* assert a completed theorem.  It upgrades the
previous negative-only audit into a guarded positive theorem-target:

    K_strict_gate as candidate completed/renormalized K_legacy_ont
    after proving the missing one-bit torsion/orientation source and the
    remaining parameter-renormalization maps.

No false pass: no silent kernel identity, no unproved beta_tors -> chi_11
collapse, no legacy physical role transfer, no QW-2191 discharge, and no ToE
closure are claimed.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_legacy_kernel_reactivation_from_diagrams_candidate_report.json"
OUT_MD = HERE / "bridge_legacy_kernel_reactivation_from_diagrams_candidate_report.md"

DIAGRAMS = ROOT / "DIAGRAMS_KERNEL_TRANSFORMATION.md"
K1 = FAR / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
S2 = FAR / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
F253 = FAR / "F253_FIRST_CURRENT_LEGACY_TO_STRICT_BRIDGE_CLOSURE_WITNESS_TARGET_PACKET.md"
TORSION_OPINION = HERE / "bridge_legacy_torsion_chi11_opinion_audit_report.json"
REYNOLDS = HERE / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.json"
META_CHARACTER = HERE / "bridge_strict_alpha_unit_character_meta_symmetry_uniqueness_obstruction_report.json"
COHOMOLOGY = HERE / "bridge_strict_alpha_unit_character_cohomology_obstruction_audit_report.json"
CONVEX_SELECTOR = HERE / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_report.json"
STRICT_FRACTAL = HERE / "bridge_strict_alpha_fractal_dimension_report.json"

TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
BALANCED_LEDGER = [2, 2, 2, 1, 1]
AUTHOR_SESSION_UPDATE = "2026-05-30 user reactivation instruction: legacy kernel may be restored as a bridge theorem-target, not silently identified."


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing text file: {path}")
    return path.read_text(encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def has_all(text: str, terms: list[str]) -> bool:
    return all(term in text for term in terms)


def diagrams_evidence(diagrams_text: str) -> list[dict[str, Any]]:
    return [
        {
            "evidence_id": "k_total_four_mechanisms",
            "claim": "DIAGRAMS records a universal K_total with geometric, resonance, torsion, and topological components.",
            "found": has_all(diagrams_text, ["K_total = K_geo × K_res × (1+0.2K_tors) × K_topo", "FOUR MECHANISMS OF K_total"]),
            "bridge_relevance": "Supports using the legacy lane as a compact nadsoliton-characteristic carrier, not merely a random fit.",
        },
        {
            "evidence_id": "effective_kernel_compression",
            "claim": "DIAGRAMS records K_total ≈ K(d) and an effective K(d)=alpha_geo*cos/(1+beta_tors*d) compression.",
            "found": has_all(diagrams_text, ["K_total ≈ K(d) with 95% accuracy", "K(d) = α_geo · cos(ω·d + φ) / (1 + β_tors·d)"]),
            "bridge_relevance": "Supports the user's statement that legacy K(d) encodes many nadsoliton characteristics as a nadsoliton-characteristic carrier before the strict completion layer.",
        },
        {
            "evidence_id": "torsion_cosine_fingerprint",
            "claim": "DIAGRAMS ties the cosine term to K_tors oscillation plus resonance.",
            "found": has_all(diagrams_text, ["cos(ωd+φ)", "From K_tors oscillation + K_res resonance"]),
            "bridge_relevance": "Makes torsion/orientation a plausible place to search for the missing chi_11 source, but not yet a proof of it.",
        },
        {
            "evidence_id": "hyperbolic_damping_fingerprint",
            "claim": "DIAGRAMS ties 1/(1+beta*d) to K_geo transformation through fractal topology and topological path summation.",
            "found": has_all(diagrams_text, ["1/(1+βd)", "fractal", "topological path summation"]),
            "bridge_relevance": "Makes the beta_tors -> strict beta*d^eta renormalization a concrete bridge slot.",
        },
        {
            "evidence_id": "beta_tors_sensitivity",
            "claim": "DIAGRAMS records beta_tors as physically sensitive in the legacy inverse hierarchy window.",
            "found": has_all(diagrams_text, ["β_tors variation", "Optimal range: β_tors ∈ [0.01, 0.08] for physics"]),
            "bridge_relevance": "Supports restoring beta_tors as an archival characteristic to be mapped, not discarded, if a bridge theorem is pursued.",
        },
    ]


def bridge_slot_table(reports: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "slot": "legacy_characteristic_carrier",
            "legacy_side": "K_total -> K_legacy_ont carries geometry/resonance/torsion/topology fingerprints.",
            "strict_side": "K_strict_gate can be treated as a candidate completed/renormalized continuation only after explicit bridge maps.",
            "current_status": "REACTIVATED_AS_THEOREM_TARGET_BY_AUTHOR_INSTRUCTION",
            "missing_object": "formal completion map B_completed_legacy_to_strict_v2",
        },
        {
            "slot": "amplitude_alpha_geo",
            "legacy_side": "explicit alpha_geo amplitude in K_legacy_ont",
            "strict_side": "strict-side alpha_geo_strict_derived_v1 := 4 ln 2 exists as strategic/source-upgrade value in S2 wording",
            "current_status": "PARTIAL_STRICT_SIDE_SOURCE_AVAILABLE_NOT_ROLE_TRANSFER",
            "missing_object": "role-safe amplitude normalization theorem",
        },
        {
            "slot": "torsion_orientation_bit",
            "legacy_side": "K_tors / beta_tors torsion fingerprint may contain orientation/chirality information",
            "strict_side": "missing chi_11/unit-axis bit has reduced algebraic support",
            "current_status": "MOST_IMPORTANT_ONE_BIT_BRIDGE_TARGET",
            "missing_object": "T_beta_tors_orientation_exports_chi11_or_nonbridge",
        },
        {
            "slot": "character_uniqueness",
            "legacy_side": "torsion source would need to prefer chi_11 over the other nontrivial U12 characters",
            "strict_side": reports["meta_character"]["interpretation"]["negative_result"],
            "current_status": "OPEN_UNIQUENESS_BLOCKER",
            "missing_object": "unit-label theorem selecting kernel {1,11}",
        },
        {
            "slot": "reynolds_obstruction_escape",
            "legacy_side": "reactivated torsion source must be orientation-bearing, not a full-Aut invariant scalar",
            "strict_side": reports["reynolds"]["interpretation"]["honest_negative"],
            "current_status": "MUST_BE_NON_FULL_AUT_OR_EXPLICITLY_ORIENTED",
            "missing_object": "proof that legacy torsion source is not annihilated by full-Aut Reynolds averaging",
        },
        {
            "slot": "eta_compression_pipeline",
            "legacy_side": "legacy hyperbolic beta_tors*d damping lacks strict eta compression",
            "strict_side": f"conditional target replay q^5={TARGET_Q_POWER}, eta={TARGET_ETA}, ledger={BALANCED_LEDGER}",
            "current_status": "CONDITIONAL_DOWNSTREAM_CERTIFICATES_EXIST",
            "missing_object": "damping/compression renormalization theorem beta_tors*d -> beta*d^eta",
        },
    ]


def build_payload() -> dict[str, Any]:
    diagrams_text = read_text(DIAGRAMS)
    k1_text = read_text(K1)
    s2_text = read_text(S2)
    f253_text = read_text(F253)
    reports = {
        "torsion_opinion": load_json(TORSION_OPINION),
        "reynolds": load_json(REYNOLDS),
        "meta_character": load_json(META_CHARACTER),
        "cohomology": load_json(COHOMOLOGY),
        "convex_selector": load_json(CONVEX_SELECTOR),
        "strict_fractal": load_json(STRICT_FRACTAL),
    }
    evidence_rows = diagrams_evidence(diagrams_text)
    slots = bridge_slot_table(reports)

    return {
        "result_kind": "SCRATCH_LEGACY_KERNEL_REACTIVATION_FROM_DIAGRAMS_CANDIDATE__THEOREM_TARGET_NOT_CLOSURE",
        "status": "legacy-kernel-restored-as-guarded-bridge-theorem-target-not-silent-identity",
        "author_session_update": AUTHOR_SESSION_UPDATE,
        "source_cross_checks": {
            "diagrams_evidence_all_found": all(row["found"] for row in evidence_rows),
            "k1_kernel_split_note_present": has_all(k1_text, ["K_legacy_ont", "K_strict_gate", "not yet rigorously identified"]),
            "s2_retirement_recorded_but_reinterpreted_by_current_author_instruction": has_all(s2_text, ["RETIRED to archival status", "strict-only ToE closure"]),
            "f253_future_bridge_target_exists": has_all(f253_text, ["Omega_legacy_strict_bridge_closure_witness_target_v1", "future-only"]),
            "chi11_reduced_bit_supported": reports["cohomology"]["character_table_audit"]["required_character"]["name"] == "chi_11_required_d5_unit_axis",
            "chi11_meta_uniqueness_still_open": "three-character ambiguity" in reports["meta_character"]["interpretation"]["direct_result"],
            "reynolds_annihilator_still_blocks_full_aut_source": reports["reynolds"]["matrix_certificate"]["reynolds_times_chi11_numerator_is_zero_matrix"],
            "conditional_balanced_ledger": reports["convex_selector"]["target_eta_9_5_certificate"]["balanced_ledger"],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "reactivation_decision": {
            "previous_s2_reading": "legacy kernel retired to archival status for strict-only work discipline",
            "current_author_update_reading": "legacy kernel may be restored as a positive bridge theorem-target because DIAGRAMS records K_total -> K(d) as a compact nadsoliton-characteristic carrier",
            "new_working_position": "Treat K_strict_gate as candidate completed/renormalized K_legacy_ont only inside an explicit bridge-target lane, with chi_11 as the leading missing bit.",
            "not_allowed": "Do not state K_legacy_ont == K_strict_gate, beta_tors == chi_11, or eta=9/5 follows from legacy torsion until the named bridge theorems are proved.",
        },
        "diagrams_evidence_rows": evidence_rows,
        "bridge_slot_table": slots,
        "one_bit_frontier": {
            "frontier_name": "T_beta_tors_orientation_exports_chi11_or_nonbridge",
            "positive_branch": "legacy torsion/topology fingerprint exports a non-full-Aut orientation source whose U12 character is chi_11",
            "negative_branch": "legacy torsion/topology fingerprint remains only archival/provenance and cannot select chi_11 without an extra axiom; this is the nonbridge theorem branch",
            "why_this_is_now_the_main_bridge_subproblem": "All loaded chi_11 audits concentrate the missing strict selector datum into a nontrivial unit character, while DIAGRAMS supplies a legacy torsion/topology place where such an orientation could be searched for explicitly.",
        },
        "theorem_target_requirements": [
            "Define the exact legacy torsion/orientation observable extracted from K_total or K_legacy_ont.",
            "Prove this observable is not full-Aut invariant, or prove a controlled symmetry reduction before applying it.",
            "Prove the resulting U12 character is uniquely chi_11, not chi_5 or chi_7.",
            "Prove the damping/compression renormalization from beta_tors*d to beta*d^eta with eta=9/5 or record a nonbridge.",
            "Prove that any use of alpha_geo/beta_tors is a bridge theorem, not silent legacy physical-role transfer.",
            "Only after those steps may K_strict_gate be called a completed K_legacy_ont rather than a guarded candidate completion.",
        ],
        "exact_proof_certificate": {
            "positive_upgrade": "Compared with the previous opinion audit, DIAGRAMS justifies restoring legacy K(d) as an active bridge theorem-target because it records K_total -> K(d) as a four-mechanism nadsoliton compression.",
            "one_bit_focus": "The leading open bridge subproblem is now sharply formulated as beta_tors/K_tors/topology orientation -> chi_11 or a nonbridge theorem.",
            "still_open": "The current repo still lacks the orientation map, chi_11 uniqueness theorem, Reynolds-obstruction escape, and eta compression renormalization theorem.",
            "guarded_language": "Use 'candidate completed/renormalized legacy kernel' for K_strict_gate until the bridge target is proved; do not use unqualified identity language.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state; K_total/K_legacy_ont may be read as a compact internal characteristic map of that nadsoliton in the bridge-target lane.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No unqualified identity K_legacy_ont == K_strict_gate is asserted.",
            "No beta_tors -> chi_11 theorem is asserted yet; it is promoted to the named main bridge theorem-target.",
            "No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.",
            "No theorem derives chi_11 uniqueness, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.",
            "No theorem derives eta=9/5 or q^5=256/243 from legacy torsion yet.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    checks = payload["source_cross_checks"]
    decision = payload["reactivation_decision"]
    frontier = payload["one_bit_frontier"]
    lines = [
        "# Legacy kernel reactivation candidate from DIAGRAMS",
        "",
        f"Status: `{payload['status']}`",
        f"Author/session update: {payload['author_session_update']}",
        "",
        "## Cross-checks",
        "",
    ]
    for key, value in checks.items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        "## Reactivation decision",
        "",
        f"- Previous S2 reading: {decision['previous_s2_reading']}",
        f"- Current author update reading: {decision['current_author_update_reading']}",
        f"- New working position: {decision['new_working_position']}",
        f"- Not allowed: {decision['not_allowed']}",
        "",
        "## DIAGRAMS evidence",
        "",
    ])
    for row in payload["diagrams_evidence_rows"]:
        lines.extend([
            f"### {row['evidence_id']}",
            f"- Claim: {row['claim']}",
            f"- Found: `{row['found']}`",
            f"- Bridge relevance: {row['bridge_relevance']}",
            "",
        ])
    lines.extend(["## Bridge slot table", ""])
    for row in payload["bridge_slot_table"]:
        lines.extend([
            f"### {row['slot']}",
            f"- Legacy side: {row['legacy_side']}",
            f"- Strict side: {row['strict_side']}",
            f"- Current status: `{row['current_status']}`",
            f"- Missing object: `{row['missing_object']}`",
            "",
        ])
    lines.extend([
        "## One-bit frontier",
        "",
        f"- Frontier: `{frontier['frontier_name']}`",
        f"- Positive branch: {frontier['positive_branch']}",
        f"- Negative branch: {frontier['negative_branch']}",
        f"- Why main bridge subproblem: {frontier['why_this_is_now_the_main_bridge_subproblem']}",
        "",
        "## Theorem-target requirements",
        "",
    ])
    lines.extend(f"- {item}" for item in payload["theorem_target_requirements"])
    lines.extend(["", "## Proof certificate", ""])
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
