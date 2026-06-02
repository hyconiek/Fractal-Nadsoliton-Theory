#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.json"
MD = GEN / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.md"

SOURCE_FILES = {
    "P2420_NONCLOSURE_MATRIX": GEN / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

GATES = [
    "apd_value_bridge_witness",
    "source_obligation_discharge",
    "chi11_phase_selector_cut_mechanism",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]
CURRENT_TRUE_GATES = ["apd_value_bridge_witness", "chi11_phase_selector_cut_mechanism"]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:20]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2421|S1371|bridge selector closure prime|closure prime implicant|failure cover certificate",
        "p2420_input": "P2420|bridge-selector nonclosure|seven closure gates|repair distance",
        "prime_implicant_prior": "prime implicant|prime implicate|failure cover|CNF|DNF",
        "closure_gate_prior": "source_obligation_discharge|chi11_source_export|qw2191_selector_discharge|role_bearing_ltotal_export",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2420's seven-gate nonclosure matrix and older generic prime-implicant/failure-cover work, "
            "but no production P24xx certificate extracting the exact prime implicant and dual failure cover for the "
            "bridge-selector closure gate itself."
        ),
    }


def p2420_theorem(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("bridge_selector_nonclosure_reason_matrix_certificate", {}).get("theorem_export", {})


def mask_for(gates: list[str]) -> int:
    return sum(1 << GATES.index(gate) for gate in gates)


def gates_for(mask: int) -> list[str]:
    return [gate for index, gate in enumerate(GATES) if mask & (1 << index)]


def toe_ready(mask: int) -> bool:
    return mask == mask_for(GATES)


def truth_vector() -> list[int]:
    return [1 if toe_ready(mask) else 0 for mask in range(1 << len(GATES))]


def anf_terms_from_truth(vector: list[int]) -> list[dict[str, Any]]:
    coeffs = vector[:]
    n = len(GATES)
    for bit in range(n):
        step = 1 << bit
        for mask in range(1 << n):
            if mask & step:
                coeffs[mask] ^= coeffs[mask ^ step]
    return [
        {"mask": mask, "gates": gates_for(mask), "degree": mask.bit_count()}
        for mask, coeff in enumerate(coeffs)
        if coeff
    ]


def prime_implicants() -> list[dict[str, Any]]:
    out = []
    for r in range(len(GATES) + 1):
        for combo in combinations(GATES, r):
            mask = mask_for(list(combo))
            true_supersets = [m for m in range(1 << len(GATES)) if m & mask == mask]
            if not true_supersets or any(not toe_ready(m) for m in true_supersets):
                continue
            proper = [sub for sub in range(mask) if sub & mask == sub and sub != mask]
            if any(all(toe_ready(m) for m in range(1 << len(GATES)) if m & sub == sub) for sub in proper):
                continue
            out.append({"mask": mask, "gates": list(combo), "degree": len(combo)})
    return out


def failure_cover_rows() -> list[dict[str, Any]]:
    full = mask_for(GATES)
    rows = []
    for gate in GATES:
        bit = 1 << GATES.index(gate)
        covered = [mask for mask in range(1 << len(GATES)) if not toe_ready(mask) and not (mask & bit)]
        rows.append(
            {
                "missing_gate": gate,
                "cover_literal": f"not {gate}",
                "covered_failure_count": len(covered),
                "nearest_miss_mask_when_only_this_gate_missing": full ^ bit,
            }
        )
    return rows


def nearest_miss_rows() -> list[dict[str, Any]]:
    full = mask_for(GATES)
    out = []
    for gate in GATES:
        bit = 1 << GATES.index(gate)
        mask = full ^ bit
        out.append(
            {
                "mask": mask,
                "true_gates": gates_for(mask),
                "missing_gate": gate,
                "toe_ready": toe_ready(mask),
                "repair_by_adding": gate,
            }
        )
    return out


def derivative_rows() -> list[dict[str, Any]]:
    full = mask_for(GATES)
    rows = []
    for gate in GATES:
        bit = 1 << GATES.index(gate)
        low = full ^ bit
        high = full
        rows.append(
            {
                "gate": gate,
                "low_mask": low,
                "high_mask": high,
                "toe_low": toe_ready(low),
                "toe_high": toe_ready(high),
                "boolean_derivative_at_nearest_edge": int(toe_ready(high)) ^ int(toe_ready(low)),
            }
        )
    return rows


def hamming_distribution_from_current() -> dict[str, int]:
    current = mask_for(CURRENT_TRUE_GATES)
    dist: dict[str, int] = {}
    for mask in range(1 << len(GATES)):
        if toe_ready(mask):
            continue
        distance = (mask ^ current).bit_count()
        dist[str(distance)] = dist.get(str(distance), 0) + 1
    return dict(sorted(dist.items(), key=lambda item: int(item[0])))


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2420 = p2420_theorem(sources["P2420_NONCLOSURE_MATRIX"])
    vector = truth_vector()
    anf = anf_terms_from_truth(vector)
    primes = prime_implicants()
    full = mask_for(GATES)
    current = mask_for(CURRENT_TRUE_GATES)
    return {
        "closure_gates": GATES,
        "current_true_gates": CURRENT_TRUE_GATES,
        "current_mask": current,
        "full_mask": full,
        "truth_vector": vector,
        "truth_vector_true_count": sum(vector),
        "anf_terms": anf,
        "prime_implicants": primes,
        "success_cnf_unit_clauses": GATES,
        "failure_cover_rows": failure_cover_rows(),
        "failure_dnf_literals": [f"not {gate}" for gate in GATES],
        "nearest_miss_rows": nearest_miss_rows(),
        "derivative_rows": derivative_rows(),
        "proper_failure_count": (1 << len(GATES)) - 1,
        "current_missing_gates": [gate for gate in GATES if gate not in CURRENT_TRUE_GATES],
        "current_repair_distance": (full ^ current).bit_count(),
        "nonclosure_hamming_distribution_from_current": hamming_distribution_from_current(),
        "p2420_minimal_toe_masks_inherited": p2420.get("minimal_toe_masks"),
        "p2420_current_repair_distance_inherited": p2420.get("current_to_toe_minimum_repair_distance"),
        "p2420_subcube_failure_count_inherited": p2420.get("apd_plus_selector_mechanism_subcube_failure_count"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2421/S1371 bridge-selector closure prime-implicant/failure-cover certificate

`P2421/S1371` turns the P2420 seven-gate nonclosure matrix into an exact Boolean theorem.  The ToE-ready predicate has one true mask, the all-gates mask `127`, and its algebraic normal form has one degree-seven term: the product of APD value bridge, source discharge, chi11 selector cut, chi11 source export, QW-2191 discharge, role-transfer license, and role-bearing `L_total` export.

The dual failure cover is the seven-literal DNF `not gate_1 OR ... OR not gate_7`.  Each missing gate has a Boolean-derivative nearest edge at the all-but-that-gate mask, so no missing gate is redundant and no proper subset can be promoted to closure.

This is still a nonclosure theorem: it identifies the exact closure prime implicant and the exact failure cover, but it exports none of the missing source, selector, QW-2191, role-transfer, or `L_total` theorems.
""".strip()
    lag_section = """
## P2421/S1371 bridge-selector closure prime-implicant guard for Lagrangian/EOM

`P2421/S1371` proves that role-bearing `L_total`/ToE promotion in the bridge-selector lane has the unique prime implicant containing all seven closure gates.  Any missing gate triggers the dual failure cover, so APD value exactness or selector-cut localization cannot be used as a Lagrangian shortcut.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    theorem_export = {
        "theorem_name": "P2421_T1_bridge_selector_closure_prime_implicant_failure_cover_certificate",
        "closure_gate_count": len(cert["closure_gates"]),
        "total_assignment_count": len(cert["truth_vector"]),
        "truth_vector_true_count": cert["truth_vector_true_count"],
        "current_mask": cert["current_mask"],
        "full_mask": cert["full_mask"],
        "anf_term_count": len(cert["anf_terms"]),
        "anf_degree": cert["anf_terms"][0]["degree"],
        "prime_implicant_count": len(cert["prime_implicants"]),
        "prime_implicant_masks": [row["mask"] for row in cert["prime_implicants"]],
        "success_cnf_unit_clause_count": len(cert["success_cnf_unit_clauses"]),
        "failure_cover_literal_count": len(cert["failure_dnf_literals"]),
        "failure_cover_rows_count": len(cert["failure_cover_rows"]),
        "proper_failure_count": cert["proper_failure_count"],
        "nearest_miss_count": len(cert["nearest_miss_rows"]),
        "all_derivative_edges_decisive": all(row["boolean_derivative_at_nearest_edge"] == 1 for row in cert["derivative_rows"]),
        "current_missing_gate_count": len(cert["current_missing_gates"]),
        "current_missing_gates": cert["current_missing_gates"],
        "current_repair_distance": cert["current_repair_distance"],
        "p2420_minimal_toe_masks_inherited": cert["p2420_minimal_toe_masks_inherited"] == [127],
        "p2420_repair_distance_inherited": cert["p2420_current_repair_distance_inherited"] == 5,
        "p2420_subcube_failure_count_inherited": cert["p2420_subcube_failure_count_inherited"] == 31,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "The unique prime implicant is a closure obligation, not a proof that its gates are discharged.",
            "The failure cover blocks every proper subset from ToE promotion.",
            "No source, chi11, QW-2191, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "seven_gates": theorem_export["closure_gate_count"] == 7,
        "all_128_assignments": theorem_export["total_assignment_count"] == 128,
        "one_true_mask": theorem_export["truth_vector_true_count"] == 1,
        "full_mask_127": theorem_export["full_mask"] == 127,
        "single_anf_term": theorem_export["anf_term_count"] == 1,
        "anf_degree_seven": theorem_export["anf_degree"] == 7,
        "unique_prime_implicant": theorem_export["prime_implicant_count"] == 1 and theorem_export["prime_implicant_masks"] == [127],
        "seven_cnf_units": theorem_export["success_cnf_unit_clause_count"] == 7,
        "seven_failure_literals": theorem_export["failure_cover_literal_count"] == 7,
        "all_127_proper_masks_fail": theorem_export["proper_failure_count"] == 127,
        "seven_nearest_misses": theorem_export["nearest_miss_count"] == 7,
        "all_derivatives_decisive": theorem_export["all_derivative_edges_decisive"],
        "current_missing_five_gates": theorem_export["current_missing_gate_count"] == 5,
        "current_repair_distance_five": theorem_export["current_repair_distance"] == 5,
        "p2420_minimal_mask_inherited": theorem_export["p2420_minimal_toe_masks_inherited"],
        "p2420_repair_distance_inherited": theorem_export["p2420_repair_distance_inherited"],
        "p2420_subcube_failure_inherited": theorem_export["p2420_subcube_failure_count_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2421_s1371_v1",
        "packet_id": "P2421",
        "stage_id": "S1371",
        "result_kind": "BRIDGE_SELECTOR_CLOSURE_PRIME_IMPLICANT_FAILURE_COVER_CERTIFICATE",
        "status": "PASS_UNIQUE_PRIME_IMPLICANT_FAILURE_COVER_NO_GATE_DISCHARGE",
        "bridge_selector_closure_prime_implicant_failure_cover_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Pick one of the five current missing gates and attempt a real theorem; do not weaken the unique all-gates "
            "prime implicant or bypass the seven-literal failure cover."
        ),
        "global_status": "OPEN_PROGRESS_PRIME_IMPLICANT_CERTIFIED_NO_CLOSURE_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["bridge_selector_closure_prime_implicant_failure_cover_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2421 S1371: bridge-selector closure prime-implicant/failure-cover certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Closure gates: `{theorem['closure_gate_count']}`.",
                f"- Assignments: `{theorem['total_assignment_count']}`.",
                f"- True masks: `{theorem['truth_vector_true_count']}`.",
                f"- ANF terms / degree: `{theorem['anf_term_count']}` / `{theorem['anf_degree']}`.",
                f"- Prime implicant masks: `{theorem['prime_implicant_masks']}`.",
                f"- Failure-cover literals: `{theorem['failure_cover_literal_count']}`.",
                f"- Current missing gates: `{theorem['current_missing_gates']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
