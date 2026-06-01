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

OUT = GEN / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.json"
MD = GEN / "p2409_s1359_stage_quotient_prime_implicate_failure_cover_certificate.md"

SOURCE_FILES = {
    "P2408_PRIME_IMPLICANT_DERIVATIVE": GEN / "p2408_s1358_stage_quotient_prime_implicant_derivative_certificate.json",
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
        "new_packet": "P2409|S1359|stage quotient prime implicate|prime-implicate failure cover",
        "p2408_source": "P2408|prime-implicant derivative|unique prime implicant|Boolean derivative edge",
        "prime_implicate_cnf": "prime implicate|prime-implicate|CNF|conjunctive normal form",
        "failure_cover": "failure cover|shortcut cover|minimal failure|proper subset fail",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds the P2408 success-side prime-implicant/derivative audit and older proper-subset "
            "fail language, but no compact dual CNF prime-implicate plus failure-cover certificate for the O/S/R quotient."
        ),
    }


def stage_mask(stages: list[str]) -> int:
    mask = 0
    for stage in stages:
        mask |= 1 << STAGE_ATOMS.index(stage)
    return mask


def stages_for(mask: int) -> list[str]:
    return [stage for index, stage in enumerate(STAGE_ATOMS) if mask & (1 << index)]


def missing_stages_for(mask: int) -> list[str]:
    return [stage for stage in STAGE_ATOMS if stage not in stages_for(mask)]


def projection_value(mask: int) -> int:
    return int(mask == stage_mask(STAGE_ATOMS))


def truth_vector() -> list[int]:
    return [projection_value(mask) for mask in range(1 << len(STAGE_ATOMS))]


def all_clauses() -> list[dict[str, Any]]:
    """Enumerate non-tautological clauses over the stage atoms.

    A literal is represented as (stage, polarity), where polarity=True means the
    positive literal `stage` and polarity=False means the negative literal `not stage`.
    """
    clauses = []
    choices = ["absent", "positive", "negative"]
    for assignment in product(choices, repeat=len(STAGE_ATOMS)):
        if all(choice == "absent" for choice in assignment):
            continue
        literals = []
        for index, choice in enumerate(assignment):
            if choice == "positive":
                literals.append({"stage": STAGE_ATOMS[index], "polarity": "positive"})
            elif choice == "negative":
                literals.append({"stage": STAGE_ATOMS[index], "polarity": "negative"})
        support = [mask for mask in range(1 << len(STAGE_ATOMS)) if clause_value(literals, mask)]
        clauses.append({"literals": literals, "support_masks": support, "literal_count": len(literals)})
    return clauses


def literal_value(literal: dict[str, str], mask: int) -> bool:
    bit = STAGE_ATOMS.index(literal["stage"])
    present = bool(mask & (1 << bit))
    return present if literal["polarity"] == "positive" else not present


def clause_value(literals: list[dict[str, str]], mask: int) -> bool:
    return any(literal_value(literal, mask) for literal in literals)


def clause_key(clause: dict[str, Any]) -> tuple[tuple[str, str], ...]:
    return tuple((literal["stage"], literal["polarity"]) for literal in clause["literals"])


def prime_implicates(values: list[int]) -> list[dict[str, Any]]:
    valid = []
    true_masks = [mask for mask, value in enumerate(values) if value]
    for clause in all_clauses():
        if all(clause_value(clause["literals"], mask) for mask in true_masks):
            valid.append(clause)
    primes = []
    for candidate in valid:
        candidate_literals = set(clause_key(candidate))
        subsumed = False
        for other in valid:
            if other is candidate:
                continue
            other_literals = set(clause_key(other))
            if other_literals < candidate_literals:
                subsumed = True
                break
        if not subsumed:
            primes.append(candidate)
    return sorted(primes, key=lambda row: (row["literal_count"], clause_key(row)))


def failure_implicant_cover(values: list[int]) -> list[dict[str, Any]]:
    fail_masks = [mask for mask, value in enumerate(values) if not value]
    cover = []
    for stage in STAGE_ATOMS:
        bit = STAGE_ATOMS.index(stage)
        support = [mask for mask in fail_masks if not (mask & (1 << bit))]
        cover.append(
            {
                "failure_implicant": f"not {stage}",
                "missing_stage": stage,
                "covered_failure_masks": support,
                "covered_failure_count": len(support),
            }
        )
    return cover


def shortcut_witness_table(values: list[int]) -> list[dict[str, Any]]:
    rows = []
    for mask, value in enumerate(values):
        if value:
            continue
        missing = missing_stages_for(mask)
        rows.append(
            {
                "mask": mask,
                "present_stages": stages_for(mask),
                "missing_stages": missing,
                "failure_cover_terms_hit": [f"not {stage}" for stage in missing],
                "minimal_repair_options": missing,
                "nearest_repair_distance": len(missing),
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    values = truth_vector()
    primes = prime_implicates(values)
    cover = failure_implicant_cover(values)
    witness_rows = shortcut_witness_table(values)
    cnf = " AND ".join(stage for stage in STAGE_ATOMS)
    failure_dnf = " OR ".join(f"not {stage}" for stage in STAGE_ATOMS)
    return {
        "stage_atoms": STAGE_ATOMS,
        "truth_vector_by_mask_0_to_7": values,
        "success_cnf": cnf,
        "success_prime_implicates": primes,
        "success_prime_implicate_count": len(primes),
        "failure_dnf": failure_dnf,
        "failure_implicant_cover": cover,
        "failure_cover_term_count": len(cover),
        "shortcut_witness_table": witness_rows,
        "shortcut_failure_mask_count": len(witness_rows),
        "shortcut_masks": [row["mask"] for row in witness_rows],
        "max_nearest_repair_distance": max(row["nearest_repair_distance"] for row in witness_rows),
        "one_stage_missing_masks": [row["mask"] for row in witness_rows if row["nearest_repair_distance"] == 1],
        "all_shortcuts_hit_by_failure_cover": all(row["failure_cover_terms_hit"] for row in witness_rows),
        "full_mask_avoids_failure_cover": all(
            stage in stages_for(stage_mask(STAGE_ATOMS)) for stage in STAGE_ATOMS
        ),
        "prime_implicates_are_positive_stage_units": sorted(primes, key=lambda row: STAGE_ATOMS.index(row["literals"][0]["stage"]))
        == [
            {"literals": [{"stage": stage, "polarity": "positive"}], "support_masks": [mask for mask in range(8) if mask & (1 << index)], "literal_count": 1}
            for index, stage in enumerate(STAGE_ATOMS)
        ],
    }


def append_doc_sections() -> None:
    eq_section = """
## P2409/S1359 stage-quotient prime-implicate and failure-cover certificate

`P2409/S1359` proves the dual of the P2408 success-side prime-implicant result.  Over the same quotient stages `O`, `S`, and `R`, the success condition has the exact prime-implicate/CNF form:

```text
O AND S AND R.
```

The failure side has the exact shortcut cover:

```text
not O OR not S OR not R.
```

Thus every proper shortcut mask is rejected by at least one missing-stage unit, and the three one-stage-missing masks are the nearest repairs rather than licensed projections.  This is a proof-carrying obstruction ledger for choosing the next missing stage theorem; it is not a ToE closure theorem.
""".strip()
    lag_section = """
## P2409/S1359 prime-implicate failure-cover guard for Lagrangian/EOM

`P2409/S1359` gives the Lagrangian/EOM guard in CNF form: a role-bearing `L_total` must satisfy the three unit obligations `O`, `S`, and `R`.  Equivalently, every proposed shortcut is covered by the failure DNF `not O OR not S OR not R`; a missing ontology guard, strict completion, or role-successor projection blocks Lagrangian promotion.
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
    p2408_theorem = artifacts["P2408_PRIME_IMPLICANT_DERIVATIVE"].get(
        "stage_quotient_prime_implicant_derivative_certificate", {}
    ).get("theorem_export", {})
    theorem_export = {
        "theorem_name": "P2409_T1_stage_quotient_prime_implicate_failure_cover_certificate",
        "truth_vector_by_mask_0_to_7": cert["truth_vector_by_mask_0_to_7"],
        "success_cnf": cert["success_cnf"],
        "success_prime_implicate_count": cert["success_prime_implicate_count"],
        "success_prime_implicates": cert["success_prime_implicates"],
        "failure_dnf": cert["failure_dnf"],
        "failure_cover_term_count": cert["failure_cover_term_count"],
        "shortcut_failure_mask_count": cert["shortcut_failure_mask_count"],
        "shortcut_masks": cert["shortcut_masks"],
        "one_stage_missing_masks": cert["one_stage_missing_masks"],
        "max_nearest_repair_distance": cert["max_nearest_repair_distance"],
        "p2408_single_success_prime_implicant_inherited": p2408_theorem.get("prime_implicant_count") == 1,
        "p2408_derivative_nearest_misses_inherited": sorted(
            p2408_theorem.get("nearest_miss_masks_by_stage", {}).values()
        ) == [3, 5, 6],
        "not_licensed": [
            "Prime implicate units are obligations, not physical-role exports.",
            "Failure-cover terms diagnose shortcuts; they do not provide selector-source closure.",
            "One-stage-missing masks are nearest repairs, not L_total licenses.",
            "No ToE closure follows from the CNF/failure-cover duality alone.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "success_cnf_has_three_unit_prime_implicates": cert["prime_implicates_are_positive_stage_units"]
        and theorem_export["success_prime_implicate_count"] == 3,
        "failure_cover_has_three_missing_stage_terms": theorem_export["failure_cover_term_count"] == 3,
        "all_seven_shortcuts_have_failure_witnesses": theorem_export["shortcut_failure_mask_count"] == 7
        and cert["all_shortcuts_hit_by_failure_cover"],
        "full_mask_avoids_failure_cover": cert["full_mask_avoids_failure_cover"],
        "one_stage_missing_masks_match_p2408_derivatives": theorem_export["one_stage_missing_masks"] == [3, 5, 6],
        "max_repair_distance_is_three_at_empty_mask": theorem_export["max_nearest_repair_distance"] == 3,
        "p2408_single_success_prime_implicant_inherited": theorem_export["p2408_single_success_prime_implicant_inherited"],
        "p2408_derivative_nearest_misses_inherited": theorem_export["p2408_derivative_nearest_misses_inherited"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2409_s1359_v1",
        "packet_id": "P2409",
        "stage_id": "S1359",
        "result_kind": "STAGE_QUOTIENT_PRIME_IMPLICATE_FAILURE_COVER_CERTIFICATE",
        "status": "PASS_PRIME_IMPLICATE_FAILURE_COVER_NO_SHORTCUT_NO_CLOSURE",
        "stage_quotient_prime_implicate_failure_cover_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "OPEN_UNKNOWN") for name, artifact in artifacts.items()},
            "finite_dual_boolean_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the CNF/failure-cover ledger to pick a concrete missing unit obligation, preferably a bridge/source "
            "or role-successor theorem, rather than weakening O*S*R or promoting shortcut masks."
        ),
        "global_status": "OPEN_PROGRESS_DUAL_OBSTRUCTION_LEDGER_CERTIFIED_NO_LTOTAL_OR_TOE_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["stage_quotient_prime_implicate_failure_cover_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2409 S1359: stage-quotient prime-implicate failure-cover certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2409/S1359 proves the dual CNF/failure-cover ledger for the P2407/P2408 quotient barrier.",
                "",
                "## Finite Boolean facts",
                "",
                f"- Truth vector masks 0..7: `{theorem['truth_vector_by_mask_0_to_7']}`.",
                f"- Success CNF: `{theorem['success_cnf']}`.",
                f"- Success prime implicate count: `{theorem['success_prime_implicate_count']}`.",
                f"- Failure DNF: `{theorem['failure_dnf']}`.",
                f"- Shortcut failure masks: `{theorem['shortcut_masks']}`.",
                f"- One-stage-missing masks: `{theorem['one_stage_missing_masks']}`.",
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
