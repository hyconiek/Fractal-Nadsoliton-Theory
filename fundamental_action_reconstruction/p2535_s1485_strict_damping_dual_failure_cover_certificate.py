#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2534_s1484_strict_damping_boolean_prime_implicant_certificate import (
    ANY,
    KEYS,
    KEY_LABELS,
    VARIABLES,
    cube_degree,
    cube_label,
    cube_literals,
    cube_matches,
    valid_assignments,
)

GEN = ROOT / "generated"
OUT = GEN / "p2535_s1485_strict_damping_dual_failure_cover_certificate.json"
MD = GEN / "p2535_s1485_strict_damping_dual_failure_cover_certificate.md"

SOURCE_FILES = {
    "P2534_BOOLEAN_PRIME_IMPLICANT": GEN / "p2534_s1484_strict_damping_boolean_prime_implicant_certificate.json",
}


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2535|S1485|dual failure cover|strict damping failure cover|strict-theorem failure cover|Boolean failure implicant",
        "intended_research_nonduplication": "failure cover|failure implicant|dual cover|one-literal relaxation|missing strict theorem cover|irredundant cover",
        "precursor_packets": "P2534|S1484|Boolean prime implicant|P2533|ternary generating polynomial|P2532|strictization distance",
        "source_obligation_language": "strict theorem|axiom-augmented|absent source|strict_damping_beta_eta_source|source-boundary",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def status_tuple(assignment: dict[str, int | str]) -> dict[str, str]:
    return {KEY_LABELS[key]: str(assignment[f"status_{key}"]) for key in KEYS}


def strict_missing_cube(key: str) -> tuple[int | str, ...]:
    cube = [ANY] * len(VARIABLES)
    cube[VARIABLES.index(f"t_{key}_strict")] = 0
    return tuple(cube)


def present_absent_cube(key: str) -> tuple[int | str, ...]:
    cube = [ANY] * len(VARIABLES)
    cube[VARIABLES.index(f"p_{key}_present")] = 0
    return tuple(cube)


def covered_by_cube(cube: tuple[int | str, ...], assignments: list[dict[str, int | str]]) -> list[dict[str, int | str]]:
    return [assignment for assignment in assignments if cube_matches(cube, assignment)]


def is_failure_implicant(cube: tuple[int | str, ...], assignments: list[dict[str, int | str]]) -> bool:
    covered = covered_by_cube(cube, assignments)
    return bool(covered) and all(assignment["strict_accept"] == 0 for assignment in covered)


def remove_literal(cube: tuple[int | str, ...], index: int) -> tuple[int | str, ...]:
    relaxed = list(cube)
    relaxed[index] = ANY
    return tuple(relaxed)


def is_prime_failure_implicant(cube: tuple[int | str, ...], assignments: list[dict[str, int | str]]) -> bool:
    if not is_failure_implicant(cube, assignments):
        return False
    for index, value in enumerate(cube):
        if value != ANY and is_failure_implicant(remove_literal(cube, index), assignments):
            return False
    return True


def all_cubes() -> list[tuple[int | str, ...]]:
    from itertools import product
    return [tuple(cube) for cube in product([ANY, 0, 1], repeat=len(VARIABLES))]


def strict_theorem_cover_rows(assignments: list[dict[str, int | str]]) -> list[dict[str, Any]]:
    rows = []
    for key in KEYS:
        cube = strict_missing_cube(key)
        covered = covered_by_cube(cube, assignments)
        rows.append({
            "cover_literal": f"t_{key}_strict=0",
            "cover_key": KEY_LABELS[key],
            "covered_failure_rows": len(covered),
            "covered_axiom_rows": sum(assignment[f"status_{key}"] == "axiom" for assignment in covered),
            "covered_absent_rows": sum(assignment[f"status_{key}"] == "absent" for assignment in covered),
            "nearest_axiom_witness": status_tuple(next(assignment for assignment in covered if assignment[f"status_{key}"] == "axiom" and sum(assignment[f"t_{other}_strict"] == 0 for other in KEYS) == 1)),
            "nearest_absent_witness": status_tuple(next(assignment for assignment in covered if assignment[f"status_{key}"] == "absent" and sum(assignment[f"t_{other}_strict"] == 0 for other in KEYS) == 1)),
        })
    return rows


def cover_intersections(assignments: list[dict[str, int | str]]) -> list[dict[str, Any]]:
    rows = []
    for size in range(1, len(KEYS) + 1):
        for subset in combinations(KEYS, size):
            covered = [assignment for assignment in assignments if all(assignment[f"t_{key}_strict"] == 0 for key in subset)]
            rows.append({
                "missing_strict_subset": [KEY_LABELS[key] for key in subset],
                "subset_size": size,
                "intersection_failure_rows": len(covered),
                "closed_form_count": (2 ** size) * (3 ** (len(KEYS) - size)),
            })
    return rows


def inclusion_exclusion_union_count(intersections: list[dict[str, Any]]) -> int:
    total = 0
    for row in intersections:
        sign = 1 if row["subset_size"] % 2 == 1 else -1
        total += sign * row["intersection_failure_rows"]
    return total


def exact_prime_failure_implicants(assignments: list[dict[str, int | str]]) -> list[dict[str, Any]]:
    primes = [cube for cube in all_cubes() if is_prime_failure_implicant(cube, assignments)]
    return [{"cube": cube_label(cube), "degree": cube_degree(cube), "literals": cube_literals(cube)} for cube in primes]


def build_failure_cover_certificate() -> dict[str, Any]:
    assignments = valid_assignments()
    true_rows = [assignment for assignment in assignments if assignment["strict_accept"] == 1]
    failure_rows = [assignment for assignment in assignments if assignment["strict_accept"] == 0]
    cover_rows = strict_theorem_cover_rows(assignments)
    intersections = cover_intersections(assignments)
    primes = exact_prime_failure_implicants(assignments)
    strict_missing_prime_cubes = {cube_label(strict_missing_cube(key)) for key in KEYS}
    absent_prime_cubes = {cube_label(present_absent_cube(key)) for key in KEYS}
    return {
        "valid_ternary_assignment_count": len(assignments),
        "strict_accept_assignment_count": len(true_rows),
        "failure_assignment_count": len(failure_rows),
        "strict_theorem_failure_cover_dnf": " OR ".join(f"(t_{key}_strict=0)" for key in KEYS),
        "strict_theorem_failure_cover_rows": cover_rows,
        "strict_theorem_failure_cover_union_count_by_inclusion_exclusion": inclusion_exclusion_union_count(intersections),
        "strict_theorem_failure_cover_intersections": intersections,
        "exact_prime_failure_implicants": primes,
        "exact_prime_failure_implicant_count": len(primes),
        "strict_missing_prime_cubes_present": strict_missing_prime_cubes.issubset({row["cube"] for row in primes}),
        "absent_prime_cubes_present_due_valid_embedding": absent_prime_cubes.issubset({row["cube"] for row in primes}),
        "nearest_failure_rows": [status_tuple(assignment) for assignment in failure_rows if sum(assignment[f"t_{key}_strict"] == 0 for key in KEYS) == 1],
    }


def append_doc_sections() -> None:
    eq_section = """
`P2535/S1485` adds the dual failure-cover audit for the P2534 Boolean strict-accept cube.  On the `3^4=81` valid ternary assignments, the strict-theorem failure cover is `t_M=0 OR t_P=0 OR t_A=0 OR t_O=0`; its four cover terms each cover `54` failure rows, and inclusion-exclusion over intersections `2^k*3^(4-k)` gives exactly `80` non-accepting rows.  Exhaustive prime-failure enumeration also records the valid-embedding absent-source primes `p_i=0`, but these are already contained in the stricter `t_i=0` failure cover.  Thus every missing strict theorem blocks acceptance; this is a dual obstruction certificate, not a source theorem.
""".strip()
    lag_section = """
`P2535/S1485` records the dual failure-cover normal form for the four-key strict damping boundary.  Any one non-strict key (`t_i=0`) blocks strict source acceptance, and the four-term cover exhausts all `80` failure rows.  The cover exports no key source, beta-eta source, `m=2` operator source, nonlinear compression-flow source, role-bearing `L_total`, bridge completion, role-transfer theorem, QW-2191 discharge, or ToE closure.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2535/S1485", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2535/S1485", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2534 = theorem(sources["P2534_BOOLEAN_PRIME_IMPLICANT"], "strict_damping_boolean_prime_implicant_certificate")
    cert = build_failure_cover_certificate()
    theorem_export = {
        "frontier_atom_under_attack": "four_key_strict_accept_dual_failure_cover",
        "p2534_boolean_prime_implicant_inherited": bool(p2534.get("boolean_prime_implicant_certificate_exported", False)),
        "source_keys": [KEY_LABELS[key] for key in KEYS],
        "valid_ternary_assignment_count": cert["valid_ternary_assignment_count"],
        "strict_accept_assignment_count": cert["strict_accept_assignment_count"],
        "failure_assignment_count": cert["failure_assignment_count"],
        "strict_theorem_failure_cover_dnf": cert["strict_theorem_failure_cover_dnf"],
        "strict_theorem_failure_cover_term_count": len(cert["strict_theorem_failure_cover_rows"]),
        "each_strict_theorem_cover_term_covers_54_rows": all(row["covered_failure_rows"] == 54 for row in cert["strict_theorem_failure_cover_rows"]),
        "each_strict_theorem_cover_term_splits_27_axiom_27_absent": all(row["covered_axiom_rows"] == 27 and row["covered_absent_rows"] == 27 for row in cert["strict_theorem_failure_cover_rows"]),
        "cover_intersections_match_closed_form": all(row["intersection_failure_rows"] == row["closed_form_count"] for row in cert["strict_theorem_failure_cover_intersections"]),
        "strict_theorem_failure_cover_union_count_by_inclusion_exclusion": cert["strict_theorem_failure_cover_union_count_by_inclusion_exclusion"],
        "strict_theorem_failure_cover_exhausts_all_failures": cert["strict_theorem_failure_cover_union_count_by_inclusion_exclusion"] == cert["failure_assignment_count"] == 80,
        "exact_prime_failure_implicant_count": cert["exact_prime_failure_implicant_count"],
        "strict_missing_prime_cubes_present": cert["strict_missing_prime_cubes_present"],
        "absent_prime_cubes_present_due_valid_embedding": cert["absent_prime_cubes_present_due_valid_embedding"],
        "nearest_failure_row_count": len(cert["nearest_failure_rows"]),
        "dual_failure_cover_certificate_exported": True,
        "axiom_promotion_to_strict_exported": False,
        "multiplicative_character_law_source_exported": False,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "strict_damping_dual_failure_cover_certificate": cert,
    }
    gatekeepers = {
        "p2534_inherited": theorem_export["p2534_boolean_prime_implicant_inherited"],
        "four_term_cover_verified": theorem_export["strict_theorem_failure_cover_term_count"] == 4 and theorem_export["each_strict_theorem_cover_term_covers_54_rows"] and theorem_export["each_strict_theorem_cover_term_splits_27_axiom_27_absent"],
        "cover_exhausts_all_failures": theorem_export["cover_intersections_match_closed_form"] and theorem_export["strict_theorem_failure_cover_exhausts_all_failures"],
        "prime_failure_implicants_audited": theorem_export["exact_prime_failure_implicant_count"] == 8 and theorem_export["strict_missing_prime_cubes_present"] and theorem_export["absent_prime_cubes_present_due_valid_embedding"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "axiom_promotion_to_strict_exported",
            "multiplicative_character_law_source_exported",
            "prime_log_proportionality_source_exported",
            "slope_value_or_prime_anchor_source_exported",
            "beta_eta_numeric_source_exported",
            "m2_operator_signature_source_exported",
            "strict_damping_beta_eta_source_exported",
            "damping_compression_bridge_component_ready",
            "full_bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "role_bearing_ltotal_exported",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2535",
        "stage_id": "S1485",
        "status": "STRICT_DAMPING_DUAL_FAILURE_COVER_CERTIFICATE_ALL_NONSTRICT_KEYS_BLOCK_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_dual_failure_cover_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_dual_failure_cover_certificate"]["theorem_export"]
    lines = [
        "# P2535/S1485 strict damping dual failure cover certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2534 Boolean prime implicant inherited: `{t['p2534_boolean_prime_implicant_inherited']}`.",
        f"- Failure assignments: `{t['failure_assignment_count']}`.",
        f"- Strict-theorem failure cover: `{t['strict_theorem_failure_cover_dnf']}`.",
        f"- Cover term count: `{t['strict_theorem_failure_cover_term_count']}`.",
        f"- Cover exhausts all failures: `{t['strict_theorem_failure_cover_exhausts_all_failures']}`.",
        f"- Exact prime failure implicant count: `{t['exact_prime_failure_implicant_count']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a dual failure-cover obstruction certificate for the conditional four-key source boundary. It does not source any key, promote axioms to strict theorems, export bridge completion, export a role-transfer theorem, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_dual_failure_cover_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_dual_failure_cover_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
