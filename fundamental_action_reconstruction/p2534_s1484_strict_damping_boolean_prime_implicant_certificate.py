#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2534_s1484_strict_damping_boolean_prime_implicant_certificate.json"
MD = GEN / "p2534_s1484_strict_damping_boolean_prime_implicant_certificate.md"

SOURCE_FILES = {
    "P2533_TERNARY_POLYNOMIAL": GEN / "p2533_s1483_strict_damping_ternary_generating_polynomial_certificate.json",
}

KEYS = ["M", "P", "A", "O"]
KEY_LABELS = {
    "M": "multiplicative_character_law_source",
    "P": "prime_log_proportionality_source",
    "A": "slope_value_or_prime_anchor_source",
    "O": "m2_operator_signature_source",
}
VARIABLES = [f"p_{key}_present" for key in KEYS] + [f"t_{key}_strict" for key in KEYS]
ANY = "*"


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
        "new_packet": "P2534|S1484|Boolean prime implicant|strict accept prime implicant|four-key Boolean minimization|strict damping Boolean boundary",
        "intended_research_nonduplication": "Boolean minimization|prime implicant|minimal implicant|strict accept cube|valid ternary embedding|essential literal",
        "precursor_packets": "P2533|S1483|ternary generating polynomial|P2532|strictization distance|P2531|four-key axiom boundary",
        "source_obligation_language": "strict theorem|axiom-augmented|absent source|strict_damping_beta_eta_source|source-boundary",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def assignment_from_statuses(statuses: tuple[str, ...]) -> dict[str, int | str]:
    assignment: dict[str, int | str] = {}
    for key, status in zip(KEYS, statuses):
        present = int(status != "absent")
        strict = int(status == "strict")
        assignment[f"p_{key}_present"] = present
        assignment[f"t_{key}_strict"] = strict
        assignment[f"status_{key}"] = status
    assignment["strict_accept"] = int(all(assignment[f"t_{key}_strict"] == 1 for key in KEYS))
    assignment["present_axiom_augmented"] = int(all(assignment[f"p_{key}_present"] == 1 for key in KEYS) and not assignment["strict_accept"])
    assignment["blocked_missing_key"] = int(any(assignment[f"p_{key}_present"] == 0 for key in KEYS))
    return assignment


def valid_assignments() -> list[dict[str, int | str]]:
    return [assignment_from_statuses(statuses) for statuses in product(["absent", "axiom", "strict"], repeat=len(KEYS))]


def cube_literals(cube: tuple[int | str, ...]) -> dict[str, int]:
    return {variable: value for variable, value in zip(VARIABLES, cube) if value != ANY}


def cube_label(cube: tuple[int | str, ...]) -> str:
    literals = cube_literals(cube)
    if not literals:
        return "TRUE"
    return " AND ".join(f"{name}={value}" for name, value in literals.items())


def cube_matches(cube: tuple[int | str, ...], assignment: dict[str, int | str]) -> bool:
    return all(value == ANY or assignment[variable] == value for variable, value in zip(VARIABLES, cube))


def covered_assignments(cube: tuple[int | str, ...], assignments: list[dict[str, int | str]]) -> list[dict[str, int | str]]:
    return [assignment for assignment in assignments if cube_matches(cube, assignment)]


def cube_degree(cube: tuple[int | str, ...]) -> int:
    return sum(value != ANY for value in cube)


def is_strict_implicant(cube: tuple[int | str, ...], assignments: list[dict[str, int | str]]) -> bool:
    covered = covered_assignments(cube, assignments)
    return bool(covered) and all(assignment["strict_accept"] == 1 for assignment in covered)


def remove_literal(cube: tuple[int | str, ...], index: int) -> tuple[int | str, ...]:
    relaxed = list(cube)
    relaxed[index] = ANY
    return tuple(relaxed)


def is_prime_implicant(cube: tuple[int | str, ...], assignments: list[dict[str, int | str]]) -> bool:
    if not is_strict_implicant(cube, assignments):
        return False
    for index, value in enumerate(cube):
        if value != ANY and is_strict_implicant(remove_literal(cube, index), assignments):
            return False
    return True


def all_cubes() -> list[tuple[int | str, ...]]:
    return [tuple(cube) for cube in product([ANY, 0, 1], repeat=len(VARIABLES))]


def status_tuple(assignment: dict[str, int | str]) -> dict[str, str]:
    return {KEY_LABELS[key]: str(assignment[f"status_{key}"]) for key in KEYS}


def omitted_literal_witnesses(prime_cube: tuple[int | str, ...], assignments: list[dict[str, int | str]]) -> list[dict[str, Any]]:
    witnesses = []
    for key in KEYS:
        index = VARIABLES.index(f"t_{key}_strict")
        relaxed = remove_literal(prime_cube, index)
        covered = covered_assignments(relaxed, assignments)
        false_rows = [assignment for assignment in covered if assignment["strict_accept"] == 0]
        axiom_witness = next(assignment for assignment in false_rows if assignment[f"status_{key}"] == "axiom")
        absent_witness = next(assignment for assignment in false_rows if assignment[f"status_{key}"] == "absent")
        witnesses.append({
            "omitted_literal": VARIABLES[index],
            "relaxed_cube": cube_label(relaxed),
            "false_rows_after_relaxation_count": len(false_rows),
            "axiom_augmented_false_witness": status_tuple(axiom_witness),
            "blocked_missing_key_false_witness": status_tuple(absent_witness),
        })
    return witnesses


def build_boolean_certificate() -> dict[str, Any]:
    assignments = valid_assignments()
    cubes = all_cubes()
    implicants = [cube for cube in cubes if is_strict_implicant(cube, assignments)]
    prime_implicants = [cube for cube in implicants if is_prime_implicant(cube, assignments)]
    minimal_degree = min(cube_degree(cube) for cube in implicants)
    minimal_implicants = [cube for cube in implicants if cube_degree(cube) == minimal_degree]
    expected_prime = tuple([ANY] * len(KEYS) + [1] * len(KEYS))
    p_only_implicants = [cube for cube in implicants if all(cube[VARIABLES.index(f"t_{key}_strict")] == ANY for key in KEYS)]
    nearest_false_rows = [assignment for assignment in assignments if sum(assignment[f"t_{key}_strict"] == 0 for key in KEYS) == 1 and all(assignment[f"p_{key}_present"] == 1 for key in KEYS)]
    return {
        "variable_order": VARIABLES,
        "valid_ternary_assignment_count": len(assignments),
        "all_boolean_cube_count": len(cubes),
        "strict_true_assignment_count": sum(int(assignment["strict_accept"]) for assignment in assignments),
        "strict_implicant_count": len(implicants),
        "strict_implicant_degree_histogram": {str(degree): sum(1 for cube in implicants if cube_degree(cube) == degree) for degree in range(len(VARIABLES) + 1)},
        "minimal_strict_implicant_degree": minimal_degree,
        "minimal_strict_implicants": [{"cube": cube_label(cube), "degree": cube_degree(cube), "literals": cube_literals(cube)} for cube in minimal_implicants],
        "prime_strict_implicants": [{"cube": cube_label(cube), "degree": cube_degree(cube), "literals": cube_literals(cube)} for cube in prime_implicants],
        "expected_prime_cube": cube_label(expected_prime),
        "expected_prime_cube_literal_count": cube_degree(expected_prime),
        "p_only_strict_implicant_count": len(p_only_implicants),
        "omitted_strict_literal_witnesses": omitted_literal_witnesses(expected_prime, assignments),
        "nearest_axiom_false_rows_count": len(nearest_false_rows),
        "nearest_axiom_false_rows": [status_tuple(assignment) for assignment in nearest_false_rows],
    }


def append_doc_sections() -> None:
    eq_section = """
`P2534/S1484` converts the P2533 ternary source-boundary polynomial into a Boolean minimization certificate over present bits `p_i` and strict-theorem bits `t_i` under the validity rule `t_i => p_i`.  Exhaustive enumeration of all `3^8=6561` Boolean cubes relative to the `3^4=81` valid ternary assignments proves that strict damping acceptance has exactly one prime/minimal implicant: `t_M=t_P=t_A=t_O=1`.  No present-only cube implies strict acceptance, and relaxing any one strict literal immediately admits explicit axiom-augmented and missing-source false witnesses.  Thus all four strict theorem literals are essential; the certificate is a minimization proof only, not a source theorem.
""".strip()
    lag_section = """
`P2534/S1484` records the Boolean prime-implicant normal form for the four-key strict damping boundary: strict source acceptance is the unique four-literal cube requiring all four source keys as strict theorems.  Present/axiom bookkeeping cannot replace any strict theorem literal, and every one-literal relaxation has false witnesses.  No key source, beta-eta source, `m=2` operator source, nonlinear compression-flow source, role-bearing `L_total`, bridge completion, role-transfer theorem, QW-2191 discharge, or ToE closure is exported.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "P2534/S1484", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "P2534/S1484", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2533 = theorem(sources["P2533_TERNARY_POLYNOMIAL"], "strict_damping_ternary_generating_polynomial_certificate")
    cert = build_boolean_certificate()
    theorem_export = {
        "frontier_atom_under_attack": "four_key_strict_accept_boolean_prime_implicant_minimization",
        "p2533_polynomial_certificate_inherited": bool(p2533.get("polynomial_certificate_exported", False)),
        "source_keys": [KEY_LABELS[key] for key in KEYS],
        "variable_order": cert["variable_order"],
        "valid_ternary_assignment_count": cert["valid_ternary_assignment_count"],
        "all_boolean_cube_count": cert["all_boolean_cube_count"],
        "strict_true_assignment_count": cert["strict_true_assignment_count"],
        "strict_implicant_count": cert["strict_implicant_count"],
        "minimal_strict_implicant_degree": cert["minimal_strict_implicant_degree"],
        "minimal_strict_implicant_degree_is_four": cert["minimal_strict_implicant_degree"] == 4,
        "minimal_strict_implicants": cert["minimal_strict_implicants"],
        "prime_strict_implicants": cert["prime_strict_implicants"],
        "unique_prime_strict_implicant": len(cert["prime_strict_implicants"]) == 1,
        "prime_implicant_is_all_four_strict_theorems": cert["prime_strict_implicants"] == [{"cube": cert["expected_prime_cube"], "degree": 4, "literals": {f"t_{key}_strict": 1 for key in KEYS}}],
        "p_only_strict_implicant_count": cert["p_only_strict_implicant_count"],
        "no_present_only_cube_implies_strict_acceptance": cert["p_only_strict_implicant_count"] == 0,
        "omitted_strict_literal_witness_count": len(cert["omitted_strict_literal_witnesses"]),
        "every_strict_literal_has_axiom_and_absent_false_witnesses": len(cert["omitted_strict_literal_witnesses"]) == len(KEYS) and all(row["false_rows_after_relaxation_count"] == 2 for row in cert["omitted_strict_literal_witnesses"]),
        "nearest_axiom_false_rows_count": cert["nearest_axiom_false_rows_count"],
        "boolean_prime_implicant_certificate_exported": True,
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
        "strict_damping_boolean_prime_implicant_certificate": cert,
    }
    gatekeepers = {
        "p2533_inherited": theorem_export["p2533_polynomial_certificate_inherited"],
        "unique_all_strict_prime_implicant": theorem_export["unique_prime_strict_implicant"] and theorem_export["prime_implicant_is_all_four_strict_theorems"],
        "no_weaker_present_only_or_lower_degree_shortcut": theorem_export["minimal_strict_implicant_degree_is_four"] and theorem_export["no_present_only_cube_implies_strict_acceptance"],
        "essential_literal_witnesses_verified": theorem_export["every_strict_literal_has_axiom_and_absent_false_witnesses"],
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
        "packet_id": "P2534",
        "stage_id": "S1484",
        "status": "STRICT_DAMPING_BOOLEAN_PRIME_IMPLICANT_CERTIFICATE_UNIQUE_ALL_STRICT_CUBE_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_boolean_prime_implicant_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_boolean_prime_implicant_certificate"]["theorem_export"]
    lines = [
        "# P2534/S1484 strict damping Boolean prime implicant certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2533 polynomial certificate inherited: `{t['p2533_polynomial_certificate_inherited']}`.",
        f"- Valid ternary assignments: `{t['valid_ternary_assignment_count']}`.",
        f"- Boolean cubes enumerated: `{t['all_boolean_cube_count']}`.",
        f"- Minimal strict implicant degree: `{t['minimal_strict_implicant_degree']}`.",
        f"- Unique prime strict implicant: `{t['unique_prime_strict_implicant']}`.",
        f"- Prime implicant is all four strict theorem bits: `{t['prime_implicant_is_all_four_strict_theorems']}`.",
        f"- Present-only strict implicants: `{t['p_only_strict_implicant_count']}`.",
        f"- Omitted strict-literal witness count: `{t['omitted_strict_literal_witness_count']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a Boolean minimization/prime-implicant certificate for the conditional four-key source boundary. It does not source any key, promote axioms to strict theorems, export bridge completion, export a role-transfer theorem, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_boolean_prime_implicant_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_boolean_prime_implicant_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
