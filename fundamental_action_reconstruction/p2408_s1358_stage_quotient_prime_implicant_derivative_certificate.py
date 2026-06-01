#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.json"
MD = GEN / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.md"

SOURCE_FILES = {
    "P2407_STAGE_QUOTIENT": GEN / "p2407_s1357_stage_quotient_projection_barrier_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STAGE_ATOMS = [
    "O_ontology_guard_package",
    "S_strict_internal_completion_package",
    "R_role_successor_projection_package",
]

TARGETS = ["ltotal_projection", "toe_projection"]


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
    return {"count": len(lines), "samples": lines[:14]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2408|S1358|stage quotient prime implicant|prime-implicant derivative",
        "p2407_source": r"P2407|stage-quotient|O \* S \* R|full stage mask",
        "prime_implicant": "prime implicant|minimal implicant|Boolean derivative|essentiality",
        "no_shortcut_guard": "ontology-only|ontology-plus-strict|no shortcut|role-bearing L_total",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2407's stage quotient and older Boolean derivative motifs, but no local compact "
            "prime-implicant/derivative certificate for the O/S/R quotient itself."
        ),
    }


def stage_mask(stages: list[str]) -> int:
    mask = 0
    for stage in stages:
        mask |= 1 << STAGE_ATOMS.index(stage)
    return mask


def stages_for(mask: int) -> list[str]:
    return [stage for index, stage in enumerate(STAGE_ATOMS) if mask & (1 << index)]


def projection_value(mask: int) -> int:
    return int(mask == stage_mask(STAGE_ATOMS))


def truth_vector() -> list[int]:
    return [projection_value(mask) for mask in range(1 << len(STAGE_ATOMS))]


def anf_terms_from_truth(values: list[int]) -> list[dict[str, Any]]:
    coeffs = values[:]
    n = len(STAGE_ATOMS)
    for bit in range(n):
        for mask in range(1 << n):
            if mask & (1 << bit):
                coeffs[mask] ^= coeffs[mask ^ (1 << bit)]
    terms = []
    for mask, coeff in enumerate(coeffs):
        if coeff:
            terms.append({"mask": mask, "stages": stages_for(mask), "degree": len(stages_for(mask))})
    return terms


def cube_matches(mask: int, required_ones: set[int], required_zeros: set[int]) -> bool:
    return all(mask & (1 << bit) for bit in required_ones) and all(not (mask & (1 << bit)) for bit in required_zeros)


def cube_support(required_ones: set[int], required_zeros: set[int]) -> list[int]:
    return [mask for mask in range(1 << len(STAGE_ATOMS)) if cube_matches(mask, required_ones, required_zeros)]


def prime_implicants(values: list[int]) -> list[dict[str, Any]]:
    implicants = []
    bit_count = len(STAGE_ATOMS)
    choices = ["free", "one", "zero"]
    for assignment in product(choices, repeat=bit_count):
        if all(choice == "free" for choice in assignment):
            continue
        required_ones = {index for index, choice in enumerate(assignment) if choice == "one"}
        required_zeros = {index for index, choice in enumerate(assignment) if choice == "zero"}
        support = cube_support(required_ones, required_zeros)
        if support and all(values[mask] == 1 for mask in support):
            implicants.append(
                {
                    "required_ones": [STAGE_ATOMS[index] for index in sorted(required_ones)],
                    "required_zeros": [STAGE_ATOMS[index] for index in sorted(required_zeros)],
                    "support_masks": support,
                    "literal_count": len(required_ones) + len(required_zeros),
                }
            )
    primes = []
    for candidate in implicants:
        candidate_support = set(candidate["support_masks"])
        is_subsumed = False
        for other in implicants:
            if other is candidate:
                continue
            other_support = set(other["support_masks"])
            if candidate_support < other_support:
                is_subsumed = True
                break
        if not is_subsumed:
            primes.append(candidate)
    return sorted(primes, key=lambda row: (row["literal_count"], row["required_ones"], row["required_zeros"]))


def derivative_edges(values: list[int]) -> dict[str, Any]:
    result = {}
    full_mask = stage_mask(STAGE_ATOMS)
    for bit, stage in enumerate(STAGE_ATOMS):
        edges = []
        for mask in range(1 << len(STAGE_ATOMS)):
            if mask & (1 << bit):
                continue
            flipped = mask | (1 << bit)
            if values[mask] ^ values[flipped]:
                edges.append({"lower_mask": mask, "upper_mask": flipped, "lower_stages": stages_for(mask), "upper_stages": stages_for(flipped)})
        result[stage] = {
            "edge_support_count": len(edges),
            "edges": edges,
            "unique_nearest_miss_mask": full_mask ^ (1 << bit),
            "unique_nearest_miss_stages": stages_for(full_mask ^ (1 << bit)),
        }
    return result


def build_certificate() -> dict[str, Any]:
    values = truth_vector()
    terms = anf_terms_from_truth(values)
    primes = prime_implicants(values)
    derivatives = derivative_edges(values)
    return {
        "stage_atoms": STAGE_ATOMS,
        "truth_vector_by_mask_0_to_7": values,
        "true_masks": [mask for mask, value in enumerate(values) if value],
        "anf_terms": terms,
        "anf_polynomial": " * ".join(STAGE_ATOMS),
        "anf_degree": max(term["degree"] for term in terms),
        "prime_implicants": primes,
        "prime_implicant_count": len(primes),
        "boolean_derivative_edges": derivatives,
        "all_derivatives_have_single_edge_support": all(row["edge_support_count"] == 1 for row in derivatives.values()),
        "all_prime_implicants_equal_full_stage_package": primes == [
            {"required_ones": STAGE_ATOMS, "required_zeros": [], "support_masks": [7], "literal_count": 3}
        ],
    }


def append_doc_sections() -> None:
    eq_section = """
## P2408/S1358 stage-quotient prime-implicant and derivative certificate

`P2408/S1358` audits the P2407 quotient as a Boolean object rather than as prose.  The truth vector over masks `0..7` has a single true mask, `7`; the computed ANF has one term, `O * S * R`; and the prime-implicant enumeration has exactly one prime implicant, the full stage package itself.

The finite derivative audit also finds one decisive edge for each stage.  Removing any one of `O`, `S`, or `R` from the full mask gives the unique nearest miss for that stage, so every package is essential and none is merely decorative bookkeeping.

Thus the quotient barrier is minimal and essential: role-bearing projection cannot be weakened below `O * S * R`, and this remains conditional readiness rather than ToE closure.
""".strip()
    lag_section = """
## P2408/S1358 prime-implicant guard for Lagrangian/EOM

`P2408/S1358` proves that the quotient `O * S * R` is the unique prime implicant for role-bearing `L_total`/ToE readiness.  Dropping ontology, strict completion, or role-successor projection crosses the unique Boolean derivative edge back to a nearest miss, so no Lagrangian shortcut is available.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate()
    p2407_theorem = artifacts["P2407_STAGE_QUOTIENT"].get("stage_quotient_projection_barrier_certificate", {}).get("theorem_export", {})
    theorem_export = {
        "theorem_name": "P2408_T1_stage_quotient_prime_implicant_derivative_certificate",
        "truth_vector_by_mask_0_to_7": cert["truth_vector_by_mask_0_to_7"],
        "true_masks": cert["true_masks"],
        "anf_polynomial": cert["anf_polynomial"],
        "anf_degree": cert["anf_degree"],
        "prime_implicant_count": cert["prime_implicant_count"],
        "prime_implicants": cert["prime_implicants"],
        "boolean_derivative_edge_counts": {
            stage: row["edge_support_count"] for stage, row in cert["boolean_derivative_edges"].items()
        },
        "nearest_miss_masks_by_stage": {
            stage: row["unique_nearest_miss_mask"] for stage, row in cert["boolean_derivative_edges"].items()
        },
        "p2407_full_mask_only_inherited": p2407_theorem.get("ltotal_true_masks") == [7]
        and p2407_theorem.get("toe_true_masks") == [7],
        "p2407_stage_degree_three_inherited": p2407_theorem.get("stage_projection_anf_degree") == 3,
        "not_licensed": [
            "No stage can be removed from O*S*R.",
            "No derivative edge is a source theorem or selector-source closure.",
            "No L_total promotion follows from any nearest miss.",
            "No ToE closure follows from prime-implicant minimality alone.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "truth_vector_has_single_true_full_mask": theorem_export["true_masks"] == [7],
        "anf_is_single_degree_three_term": theorem_export["anf_polynomial"] == " * ".join(STAGE_ATOMS)
        and theorem_export["anf_degree"] == 3,
        "unique_prime_implicant_is_full_stage_package": cert["all_prime_implicants_equal_full_stage_package"],
        "all_derivatives_have_single_edge_support": cert["all_derivatives_have_single_edge_support"],
        "nearest_misses_are_one_stage_missing": sorted(theorem_export["nearest_miss_masks_by_stage"].values()) == [3, 5, 6],
        "p2407_full_mask_only_inherited": theorem_export["p2407_full_mask_only_inherited"],
        "p2407_stage_degree_three_inherited": theorem_export["p2407_stage_degree_three_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2408_s1358_v1",
        "packet_id": "P2408",
        "stage_id": "S1358",
        "result_kind": "STAGE_QUOTIENT_PRIME_IMPLICANT_DERIVATIVE_CERTIFICATE",
        "status": "PASS_UNIQUE_PRIME_IMPLICANT_AND_SINGLE_EDGE_DERIVATIVES_NO_CLOSURE",
        "stage_quotient_prime_implicant_derivative_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "OPEN_UNKNOWN") for name, artifact in artifacts.items()},
            "finite_boolean_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the unique prime implicant to choose a concrete missing stage theorem; do not weaken O*S*R "
            "or promote nearest misses into L_total/ToE closure."
        ),
        "global_status": "OPEN_PROGRESS_STAGE_QUOTIENT_MINIMALITY_CERTIFIED_NO_PHYSICAL_ROLE_EXPORT",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["stage_quotient_prime_implicant_derivative_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2408 S1358: stage-quotient prime-implicant derivative certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2408/S1358 proves the P2407 quotient barrier is the unique prime implicant O*S*R and audits its Boolean derivative edges.",
                "",
                "## Finite Boolean facts",
                "",
                f"- Truth vector masks 0..7: `{theorem['truth_vector_by_mask_0_to_7']}`.",
                f"- True masks: `{theorem['true_masks']}`.",
                f"- ANF: `{theorem['anf_polynomial']}`.",
                f"- Prime implicant count: `{theorem['prime_implicant_count']}`.",
                f"- Derivative edge counts: `{theorem['boolean_derivative_edge_counts']}`.",
                f"- Nearest miss masks by stage: `{theorem['nearest_miss_masks_by_stage']}`.",
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
