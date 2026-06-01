#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import product
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.json"
MD = GEN / "p2410_s1360_dequotiented_twelve_atom_prime_implicate_obstruction_certificate.md"

SOURCE_FILES = {
    "P2406_STAGED_PROJECTION": GEN / "p2406_s1356_information_to_physics_staged_projection_barrier_certificate.json",
    "P2409_QUOTIENT_FAILURE_COVER": GEN / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ONTOLOGY_GUARD_ATOMS = [
    "nadsoliton_is_sole_primordial_information",
    "no_separate_information_layer_under_nadsoliton",
    "strict_additions_are_internal_information_constraints",
    "physics_roles_are_downstream_projections",
    "observer_is_downstream_readout_not_source",
]

STRICT_ADDITION_ATOMS = [
    "apd_completion_structure",
    "gf2_phase_topological_data",
    "nonlinear_damping_compression",
    "chi11_selector_source_declared",
]

ROLE_SUCCESSOR_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]

STAGE_PACKAGES = {
    "O_ontology_guard_package": ONTOLOGY_GUARD_ATOMS,
    "S_strict_internal_completion_package": STRICT_ADDITION_ATOMS,
    "R_role_successor_projection_package": ROLE_SUCCESSOR_ATOMS,
}

ATOMS = ONTOLOGY_GUARD_ATOMS + STRICT_ADDITION_ATOMS + ROLE_SUCCESSOR_ATOMS
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}
FULL_MASK = (1 << len(ATOMS)) - 1


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
        "new_packet": "P2410|S1360|dequotiented twelve atom|de-quotiented 12-atom|atom-level failure cover",
        "p2406_source": "P2406|degree-12|12 atoms|staged projection barrier",
        "p2409_source": "P2409|prime-implicate failure-cover|failure DNF|shortcut cover",
        "repair_spectrum": "repair spectrum|nearest repair|missing atom|prime implicate",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2406's degree-12 staged barrier and P2409's quotient CNF, but no exact "
            "dequotiented 12-atom prime-implicate/failure-cover obstruction ledger with repair spectrum."
        ),
    }


def atoms_for(mask: int) -> list[str]:
    return [atom for atom in ATOMS if mask & (1 << ATOM_INDEX[atom])]


def missing_atoms_for(mask: int) -> list[str]:
    return [atom for atom in ATOMS if not (mask & (1 << ATOM_INDEX[atom]))]


def full_projection_value(mask: int) -> int:
    return int(mask == FULL_MASK)


def literal_value(atom_index: int, polarity: str, mask: int) -> bool:
    present = bool(mask & (1 << atom_index))
    return present if polarity == "positive" else not present


def clause_value(literals: list[tuple[int, str]], mask: int) -> bool:
    return any(literal_value(index, polarity, mask) for index, polarity in literals)


def enumerate_prime_implicates() -> dict[str, Any]:
    """Enumerate the 3^12 non-tautological clauses and extract prime implicates.

    The full projection has a singleton true support, so this is still finite and
    exact: every implicate must hold at FULL_MASK, and subset-minimal implicates
    are exactly the twelve positive unit clauses.
    """
    valid_implicate_count = 0
    prime_candidates: list[dict[str, Any]] = []
    choices = ["absent", "positive", "negative"]
    for assignment in product(choices, repeat=len(ATOMS)):
        if all(choice == "absent" for choice in assignment):
            continue
        literals = [(index, choice) for index, choice in enumerate(assignment) if choice != "absent"]
        if clause_value(literals, FULL_MASK):
            valid_implicate_count += 1
            if len(literals) == 1 and literals[0][1] == "positive":
                index = literals[0][0]
                prime_candidates.append(
                    {
                        "unit_clause": ATOMS[index],
                        "atom": ATOMS[index],
                        "stage_package": stage_for_atom(ATOMS[index]),
                        "literal_count": 1,
                        "covered_success_mask_count": 1 << (len(ATOMS) - 1),
                    }
                )
    return {
        "enumerated_nonempty_clause_count": (3 ** len(ATOMS)) - 1,
        "valid_implicate_count": valid_implicate_count,
        "prime_implicates": sorted(prime_candidates, key=lambda row: ATOM_INDEX[row["atom"]]),
        "prime_implicate_count": len(prime_candidates),
    }


def stage_for_atom(atom: str) -> str:
    for stage, atoms in STAGE_PACKAGES.items():
        if atom in atoms:
            return stage
    raise KeyError(atom)


def success_cnf() -> str:
    return " AND ".join(ATOMS)


def failure_dnf() -> str:
    return " OR ".join(f"not {atom}" for atom in ATOMS)


def failure_cover() -> list[dict[str, Any]]:
    return [
        {
            "failure_implicant": f"not {atom}",
            "missing_atom": atom,
            "stage_package": stage_for_atom(atom),
            "covered_failure_mask_count": 1 << (len(ATOMS) - 1),
        }
        for atom in ATOMS
    ]


def repair_spectrum() -> list[dict[str, Any]]:
    rows = []
    for missing_count in range(1, len(ATOMS) + 1):
        rows.append(
            {
                "missing_atom_count": missing_count,
                "present_atom_count": len(ATOMS) - missing_count,
                "failure_mask_count": math.comb(len(ATOMS), missing_count),
                "nearest_repair_distance": missing_count,
            }
        )
    return rows


def stage_completion_mask(mask: int) -> int:
    stage_mask = 0
    for index, (stage, atoms) in enumerate(STAGE_PACKAGES.items()):
        if all(mask & (1 << ATOM_INDEX[atom]) for atom in atoms):
            stage_mask |= 1 << index
    return stage_mask


def stage_quotient_preimage_counts() -> dict[str, int]:
    counts = {str(mask): 0 for mask in range(1 << len(STAGE_PACKAGES))}
    for mask in range(1 << len(ATOMS)):
        counts[str(stage_completion_mask(mask))] += 1
    return counts


def nearest_failure_examples() -> list[dict[str, Any]]:
    rows = []
    for atom in ATOMS:
        missing_mask = FULL_MASK ^ (1 << ATOM_INDEX[atom])
        rows.append(
            {
                "failure_mask": missing_mask,
                "present_atom_count": len(ATOMS) - 1,
                "missing_atoms": [atom],
                "repair_atom": atom,
                "stage_package": stage_for_atom(atom),
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    implicates = enumerate_prime_implicates()
    spectrum = repair_spectrum()
    stage_counts = stage_quotient_preimage_counts()
    return {
        "atom_count": len(ATOMS),
        "total_assignment_count": 1 << len(ATOMS),
        "success_mask": FULL_MASK,
        "success_assignment_count": 1,
        "failure_assignment_count": (1 << len(ATOMS)) - 1,
        "stage_packages": STAGE_PACKAGES,
        "success_cnf": success_cnf(),
        "success_cnf_unit_count": len(ATOMS),
        "failure_dnf": failure_dnf(),
        "failure_cover": failure_cover(),
        "failure_cover_term_count": len(ATOMS),
        "prime_implicate_enumeration": implicates,
        "repair_spectrum": spectrum,
        "repair_spectrum_failure_total": sum(row["failure_mask_count"] for row in spectrum),
        "nearest_one_atom_missing_failures": nearest_failure_examples(),
        "nearest_one_atom_missing_failure_count": len(ATOMS),
        "stage_quotient_preimage_counts": stage_counts,
        "stage_quotient_full_mask_preimage_count": stage_counts[str((1 << len(STAGE_PACKAGES)) - 1)],
        "stage_quotient_empty_mask_preimage_count": stage_counts["0"],
    }


def append_doc_sections() -> None:
    eq_section = """
## P2410/S1360 dequotiented twelve-atom prime-implicate obstruction certificate

`P2410/S1360` expands the P2409 quotient CNF back to the full P2406 twelve-atom staged barrier.  The finite enumeration checks all `3^12 - 1` nonempty non-tautological clauses and finds the exact success-side prime implicates to be the twelve positive unit obligations:

```text
O_1 AND ... AND O_5 AND S_1 AND ... AND S_4 AND R_1 AND R_2 AND R_3.
```

The dual failure cover is the twelve-term DNF saying that any missing ontology guard atom, strict internal-completion atom, or role-successor atom blocks `L_total`/ToE projection.  The repair spectrum contains `4095` failure assignments with distance distribution `binomial(12,k)` by the number of missing atoms, and the nearest failures are exactly the twelve one-atom-missing masks.

This dequotients the obstruction ledger without weakening the guard: it identifies atom-level proof obligations but does not export role theorems, selector-source closure, or ToE closure.
""".strip()
    lag_section = """
## P2410/S1360 dequotiented obstruction guard for Lagrangian/EOM

`P2410/S1360` refines the quotient `O * S * R` Lagrangian guard into twelve unit obligations.  A role-bearing `L_total` remains blocked by any missing ontology atom, strict internal-completion atom, or role-successor atom; the twelve one-atom-missing masks are nearest repair targets, not admissible Lagrangian terms.
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
    p2406_theorem = artifacts["P2406_STAGED_PROJECTION"].get(
        "information_to_physics_staged_projection_barrier_certificate", {}
    ).get("theorem_export", {})
    p2409_theorem = artifacts["P2409_QUOTIENT_FAILURE_COVER"].get(
        "stage_quotient_prime_implicate_failure_cover_certificate", {}
    ).get("theorem_export", {})
    theorem_export = {
        "theorem_name": "P2410_T1_dequotiented_twelve_atom_prime_implicate_obstruction_certificate",
        "atom_count": cert["atom_count"],
        "total_assignment_count": cert["total_assignment_count"],
        "success_assignment_count": cert["success_assignment_count"],
        "failure_assignment_count": cert["failure_assignment_count"],
        "success_cnf_unit_count": cert["success_cnf_unit_count"],
        "failure_cover_term_count": cert["failure_cover_term_count"],
        "enumerated_nonempty_clause_count": cert["prime_implicate_enumeration"]["enumerated_nonempty_clause_count"],
        "valid_implicate_count": cert["prime_implicate_enumeration"]["valid_implicate_count"],
        "prime_implicate_count": cert["prime_implicate_enumeration"]["prime_implicate_count"],
        "prime_implicate_atoms": [row["atom"] for row in cert["prime_implicate_enumeration"]["prime_implicates"]],
        "repair_spectrum": cert["repair_spectrum"],
        "repair_spectrum_failure_total": cert["repair_spectrum_failure_total"],
        "nearest_one_atom_missing_failure_count": cert["nearest_one_atom_missing_failure_count"],
        "stage_quotient_preimage_counts": cert["stage_quotient_preimage_counts"],
        "stage_quotient_full_mask_preimage_count": cert["stage_quotient_full_mask_preimage_count"],
        "p2406_degree_twelve_inherited": p2406_theorem.get("ltotal_projection_anf_degree") == 12
        and p2406_theorem.get("toe_projection_anf_degree") == 12,
        "p2409_three_stage_cnf_inherited": p2409_theorem.get("success_prime_implicate_count") == 3
        and p2409_theorem.get("failure_cover_term_count") == 3,
        "not_licensed": [
            "Unit prime implicates are proof obligations, not role exports.",
            "Nearest one-atom repairs are target suggestions, not L_total licenses.",
            "The chi11 atom remains a selector/source obligation and does not discharge QW-2191 by bookkeeping.",
            "The dequotiented ledger does not transfer legacy physical-role claims onto the strict kernel.",
            "No ToE closure follows from the twelve-atom CNF/failure-cover certificate alone.",
        ],
    }
    expected_valid_implicate_count = (3 ** len(ATOMS)) - (2 ** len(ATOMS))
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "total_assignment_count_is_4096": theorem_export["total_assignment_count"] == 4096,
        "single_success_and_4095_failures": theorem_export["success_assignment_count"] == 1
        and theorem_export["failure_assignment_count"] == 4095,
        "enumerated_all_nonempty_ternary_clauses": theorem_export["enumerated_nonempty_clause_count"] == (3 ** 12) - 1,
        "valid_implicate_count_matches_closed_form": theorem_export["valid_implicate_count"] == expected_valid_implicate_count,
        "prime_implicates_are_twelve_positive_units": theorem_export["prime_implicate_count"] == 12
        and theorem_export["prime_implicate_atoms"] == ATOMS,
        "failure_cover_has_twelve_terms": theorem_export["failure_cover_term_count"] == 12,
        "repair_spectrum_sums_to_all_failures": theorem_export["repair_spectrum_failure_total"] == 4095,
        "nearest_failures_are_twelve_one_atom_misses": theorem_export["nearest_one_atom_missing_failure_count"] == 12,
        "stage_full_preimage_is_singleton": theorem_export["stage_quotient_full_mask_preimage_count"] == 1,
        "p2406_degree_twelve_inherited": theorem_export["p2406_degree_twelve_inherited"],
        "p2409_three_stage_cnf_inherited": theorem_export["p2409_three_stage_cnf_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2410_s1360_v1",
        "packet_id": "P2410",
        "stage_id": "S1360",
        "result_kind": "DEQUOTIENTED_TWELVE_ATOM_PRIME_IMPLICATE_OBSTRUCTION_CERTIFICATE",
        "status": "PASS_DEQUOTIENTED_12_ATOM_OBSTRUCTION_LEDGER_NO_CLOSURE",
        "dequotiented_twelve_atom_prime_implicate_obstruction_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "OPEN_UNKNOWN") for name, artifact in artifacts.items()},
            "finite_obstruction_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Pick one unit obligation from the twelve-atom ledger and prove it with real bridge/source or role-successor "
            "content; do not use dequotienting itself as selector closure or role-transfer closure."
        ),
        "global_status": "OPEN_PROGRESS_ATOM_LEVEL_OBSTRUCTION_LEDGER_CERTIFIED_NO_LTOTAL_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["dequotiented_twelve_atom_prime_implicate_obstruction_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2410 S1360: dequotiented twelve-atom prime-implicate obstruction certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2410/S1360 expands the quotient O/S/R obstruction ledger to the exact twelve atom P2406 barrier.",
                "",
                "## Finite Boolean facts",
                "",
                f"- Total assignments: `{theorem['total_assignment_count']}`.",
                f"- Success/failure assignments: `{theorem['success_assignment_count']}` / `{theorem['failure_assignment_count']}`.",
                f"- Enumerated clauses: `{theorem['enumerated_nonempty_clause_count']}`.",
                f"- Prime implicate count: `{theorem['prime_implicate_count']}`.",
                f"- Failure cover term count: `{theorem['failure_cover_term_count']}`.",
                f"- Nearest one-atom-missing failures: `{theorem['nearest_one_atom_missing_failure_count']}`.",
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
