#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2612_s1562_p2607_gf2_physical_origin_obstruction.json"
MD = GEN / "p2612_s1562_p2607_gf2_physical_origin_obstruction.md"

SOURCE_FILES = {
    "P2607_SOURCE": ROOT / "p2607_s1557_strict_phase_topological_selector_bridge_completion.py",
    "P2607_PAYLOAD": GEN / "p2607_s1557_strict_phase_topological_selector_bridge_completion.json",
    "P2610_REVALIDATION": GEN / "p2610_s1560_p2601_p2608_critical_revalidation_audit.json",
    "P2611_ROLE_SEMANTICS": GEN / "p2611_s1561_ltotal_role_semantics_acceptance_predicate.json",
}

PHYSICAL_ORIGIN_REQUIREMENTS = [
    "chiral_current_conservation_equations",
    "topological_charge_quantization_map",
    "winding_number_to_phase_bit_boundary_operator",
    "node_label_invariance_or_perturbation_invariance_test",
]

PHYSICAL_KEYWORDS = [
    "chiral",
    "current",
    "topological charge",
    "winding",
    "homology",
    "cohomology",
    "boundary operator",
    "orientation class",
]

NEGATIVE_EXPORT_FLAGS = [
    "p2607_quarantine_lifted",
    "gf2_physical_origin_theorem_exported",
    "bridge_completion_revalidated",
    "p2608_role_bearing_ltotal_reenabled",
    "strict_damping_source_reexported",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
    "apd_source_exported",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2612|S1562|GF\\(2\\) physical-origin obstruction|GF2 physical-origin obstruction|P2607.*physical-origin obstruction",
        "intended_research_nonduplication": "chiral current obstruction|winding number obstruction|topological charge obstruction|P2601.*monoid action uniqueness",
        "p2607_origin_terms": "P2607|S1557|GF\\(2\\)|chiral current|topological charge|winding number|phase/topological selector",
        "guardrails": "K_legacy_ont|K_strict_gate|QW-2191|role-bearing L_total|ToE closure|APD source",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def gf2_rank(matrix: list[list[int]]) -> int:
    rows = [row[:] for row in matrix]
    if not rows:
        return 0
    rank = 0
    col = 0
    while rank < len(rows) and col < len(rows[0]):
        pivot = next((idx for idx in range(rank, len(rows)) if rows[idx][col] == 1), None)
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for idx in range(len(rows)):
            if idx != rank and rows[idx][col] == 1:
                rows[idx] = [(a ^ b) for a, b in zip(rows[idx], rows[rank])]
        rank += 1
        col += 1
    return rank


def mat_vec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(a & b for a, b in zip(row, vector)) % 2 for row in matrix]


def is_lower_triangular_unit(matrix: list[list[int]]) -> bool:
    return all(
        row[i] == 1 and all(row[j] == 0 for j in range(i + 1, len(row)))
        for i, row in enumerate(matrix)
    )


def row_supports(matrix: list[list[int]]) -> list[list[int]]:
    return [[idx for idx, value in enumerate(row) if value] for row in matrix]


def source_origin_audit(source_text: str) -> dict[str, Any]:
    lower_triangular_construction_markers = {
        "unit_diagonal_assignment": "row[i] = 1" in source_text,
        "nearest_neighbor_index_rule": "if i > 0" in source_text and "row[i - 1] = 1" in source_text,
        "mod_three_extra_coupling_rule": "i % 3 == 0" in source_text and "row[i - 2] = 1" in source_text,
    }
    keyword_hits = {keyword: (keyword.lower() in source_text.lower()) for keyword in PHYSICAL_KEYWORDS}
    return {
        "lower_triangular_construction_markers": lower_triangular_construction_markers,
        "physical_keyword_hits_in_p2607_source": keyword_hits,
        "physical_keyword_hit_count": sum(keyword_hits.values()),
        "all_algorithmic_index_markers_present": all(lower_triangular_construction_markers.values()),
    }


def physical_data_audit(payload_text: str) -> dict[str, Any]:
    requirement_hits = {}
    for req in PHYSICAL_ORIGIN_REQUIREMENTS:
        words = req.replace("_", ".*")
        requirement_hits[req] = re.search(words, payload_text, flags=re.IGNORECASE | re.DOTALL) is not None
    return {
        "requirement_hits": requirement_hits,
        "missing_requirements": [req for req, hit in requirement_hits.items() if not hit],
        "all_physical_requirements_present": all(requirement_hits.values()),
    }


def nullity_under_rhs_sweep(matrix: list[list[int]], limit: int = 16) -> dict[str, Any]:
    n = len(matrix)
    samples = []
    for bits in product([0, 1], repeat=min(n, 4)):
        vector = list(bits) + [0] * (n - len(bits))
        samples.append({
            "phase_bits_prefix": list(bits),
            "rhs_prefix": mat_vec_gf2(matrix, vector)[:4],
            "rank_over_gf2": gf2_rank(matrix),
            "nullity_over_gf2": n - gf2_rank(matrix),
        })
        if len(samples) >= limit:
            break
    return {
        "sample_count": len(samples),
        "rank_independent_of_rhs_samples": all(sample["rank_over_gf2"] == n for sample in samples),
        "samples": samples,
    }


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2607_payload = load_json(SOURCE_FILES["P2607_PAYLOAD"])
    p2610_payload = load_json(SOURCE_FILES["P2610_REVALIDATION"])
    p2611_payload = load_json(SOURCE_FILES["P2611_ROLE_SEMANTICS"])
    p2607_source = SOURCE_FILES["P2607_SOURCE"].read_text(encoding="utf-8")
    p2607_theorem = theorem(p2607_payload, "strict_phase_topological_selector_bridge_completion")
    p2610_theorem = theorem(p2610_payload, "p2601_p2608_critical_revalidation_audit")
    p2611_theorem = theorem(p2611_payload, "ltotal_role_semantics_acceptance_predicate")
    matrix = p2607_theorem["phase_topological_selector_audit"]["gf2_selector_matrix"]
    rank = gf2_rank(matrix)
    nullity = len(matrix) - rank
    origin = source_origin_audit(p2607_source)
    physical = physical_data_audit(json.dumps(p2607_theorem, sort_keys=True) + "\n" + p2607_source)
    rhs_sweep = nullity_under_rhs_sweep(matrix)
    obstruction = {
        "obstruction_statement": (
            "The P2607 GF(2) matrix cannot be lifted out of quarantine from the current repository evidence. Its full rank follows from an algorithmic lower-triangular unit construction over node indices, while no derivation from chiral-current conservation, topological-charge quantization, or nadsoliton winding-number boundary data is present. Therefore the P2607 bridge-completion export remains obstructed."
        ),
        "matrix_size": [len(matrix), len(matrix[0]) if matrix else 0],
        "rank_over_gf2_recomputed": rank,
        "nullity_over_gf2_recomputed": nullity,
        "lower_triangular_unit_matrix": is_lower_triangular_unit(matrix),
        "row_supports": row_supports(matrix),
        "rank_explained_by_unit_diagonal": is_lower_triangular_unit(matrix) and rank == len(matrix),
        "rhs_rank_independence_audit": rhs_sweep,
        "source_origin_audit": origin,
        "physical_origin_requirement_audit": physical,
        "p2610_quarantine_inherited": {
            "p2607_quarantined": "P2607" in set(p2610_theorem.get("quarantined_packet_ids_after_revalidation", [])),
            "p2608_quarantined": "P2608" in set(p2610_theorem.get("quarantined_packet_ids_after_revalidation", [])),
        },
        "p2611_role_semantics_inherited": {
            "role_semantics_defined": p2611_theorem.get("role_semantics_definition_exported") is True,
            "current_ltotal_assignment_accepts": p2611_theorem.get("acceptance_truth_table", {}).get("current_assignment_accepts") is True,
        },
        "strict_alternative_proposal": {
            "preferred_path": "P2601_monoid_action_uniqueness",
            "why_preferred_over_p2602_now": "The repository already frames T_1 as zero RG time/no-op transport; a strict monoid proof needs fewer new physical assumptions than a prime spectral-gap theorem selecting primes over all other generator families.",
            "required_steps": [
                "Define the nadsoliton RG transport action T_d on the admissible state space with composition T_de = T_d o T_e.",
                "Prove identity uniqueness: T_1 = Id is the unique neutral element of that action.",
                "Define damping character y_d as the logarithm of the multiplicative attenuation character and prove y_de = y_d + y_e.",
                "Use identity uniqueness to prove y_1 = 0 and show the result is independent of D_f.",
            ],
        },
    }
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2612_T1_p2607_gf2_physical_origin_obstruction",
        "audited_chain": ["P2607/S1557", "P2610/S1560", "P2611/S1561"],
        "gf2_physical_origin_obstruction_exported": True,
        "obstruction_certificate": obstruction,
        "recommended_next_honest_step": (
            "Close the P2607 GF(2) bridge path until a physical-origin theorem is supplied; pursue the P2601 monoid action uniqueness proof as the strict alternative."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2607_rank_recomputed": theorem_export["obstruction_certificate"]["rank_over_gf2_recomputed"] == 11,
        "unit_diagonal_explains_rank": theorem_export["obstruction_certificate"]["rank_explained_by_unit_diagonal"],
        "physical_requirements_missing": not theorem_export["obstruction_certificate"]["physical_origin_requirement_audit"]["all_physical_requirements_present"],
        "p2607_quarantine_inherited": theorem_export["obstruction_certificate"]["p2610_quarantine_inherited"]["p2607_quarantined"],
        "p2608_quarantine_inherited": theorem_export["obstruction_certificate"]["p2610_quarantine_inherited"]["p2608_quarantined"],
        "role_semantics_does_not_accept_current_ltotal": theorem_export["obstruction_certificate"]["p2611_role_semantics_inherited"]["current_ltotal_assignment_accepts"] is False,
        "obstruction_exported": theorem_export["gf2_physical_origin_obstruction_exported"],
        "no_p2607_quarantine_lift": theorem_export["p2607_quarantine_lifted"] is False,
        "no_bridge_revalidation": theorem_export["bridge_completion_revalidated"] is False,
        "no_p2608_ltotal_reenable": theorem_export["p2608_role_bearing_ltotal_reenabled"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2612",
        "stage_id": "S1562",
        "status": "P2612_P2607_GF2_PHYSICAL_ORIGIN_OBSTRUCTION_CLOSES_BRIDGE_LIFT_PATH_PROPOSES_P2601_MONOID_UNIQUENESS_ALTERNATIVE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "p2607_gf2_physical_origin_obstruction": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2607_SOURCE": sha256_text(p2607_source),
                "P2607_PAYLOAD": sha256_json(p2607_payload),
                "P2610_REVALIDATION": sha256_json(p2610_payload),
                "P2611_ROLE_SEMANTICS": sha256_json(p2611_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2607_gf2_physical_origin_obstruction"]["theorem_export"]
    obs = t["obstruction_certificate"]
    lines = [
        "# P2612/S1562 P2607 GF(2) physical-origin obstruction", "",
        f"Status: `{payload['status']}`", "", "## Obstruction", "",
        obs["obstruction_statement"], "", "## Computed audit", "",
        f"- Recomputed GF(2) rank: `{obs['rank_over_gf2_recomputed']}`.",
        f"- Recomputed nullity: `{obs['nullity_over_gf2_recomputed']}`.",
        f"- Lower-triangular unit matrix: `{obs['lower_triangular_unit_matrix']}`.",
        f"- Rank explained by unit diagonal: `{obs['rank_explained_by_unit_diagonal']}`.",
        f"- Physical requirements all present: `{obs['physical_origin_requirement_audit']['all_physical_requirements_present']}`.",
        f"- Missing physical requirements: `{obs['physical_origin_requirement_audit']['missing_requirements']}`.",
        f"- P2607 quarantined inherited: `{obs['p2610_quarantine_inherited']['p2607_quarantined']}`.",
        f"- P2608 quarantined inherited: `{obs['p2610_quarantine_inherited']['p2608_quarantined']}`.", "",
        "## Alternative", "",
        f"Preferred path: `{obs['strict_alternative_proposal']['preferred_path']}`.",
        obs["strict_alternative_proposal"]["why_preferred_over_p2602_now"], "",
        "## Scope guards", "",
        "P2612 does not lift P2607 quarantine, does not revalidate bridge completion, does not re-enable P2608 role-bearing L_total, and does not export QW-2191 discharge, APD sourcehood, legacy physical-role transfer, or ToE closure.", "", "## Fingerprint", "",
        f"`{payload['p2607_gf2_physical_origin_obstruction']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2612/S1562 P2607 GF(2) physical-origin obstruction

`P2612/S1562` attempts the requested physical-origin lift for the P2607 GF(2) phase/selector matrix and closes the path with an obstruction: the audited full rank follows from an index-built lower-triangular unit matrix, while no derivation from chiral currents, topological charge, winding-number boundary data, or node-label invariance is present.  P2607 bridge completion and P2608 role-bearing `L_total` therefore remain quarantined; the recommended alternative is a strict P2601 monoid action uniqueness proof.
""".strip()
    lag_section = """
## P2612/S1562 GF(2)-origin obstruction Ltotal guard

`P2612/S1562` keeps `L_total` role-bearing acceptance blocked: the P2607 GF(2) matrix has no audited physical origin from nadsoliton topology in the current evidence, so P2608 cannot be re-enabled.  The next admissible proof route is the P2601 monoid action uniqueness theorem, not another combinatorial bridge lift.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2612/S1562 P2607 GF(2) physical-origin obstruction", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2612/S1562 GF(2)-origin obstruction Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
