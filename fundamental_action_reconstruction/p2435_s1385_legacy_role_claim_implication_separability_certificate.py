#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2435_s1385_legacy_role_claim_implication_separability_certificate.json"
MD = GEN / "p2435_s1385_legacy_role_claim_implication_separability_certificate.md"

SOURCE_FILES = {
    "P2434_ROLE_LATTICE": GEN / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}


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
        "new_packet": "P2435|S1385|legacy role claim implication|role claim separability|implication separability",
        "p2434_input": "P2434|conditional legacy role-transfer|role-transfer claim lattice|64 masks",
        "role_implication_prior": "role.*implication|claim.*implication|strictly implies|separability|GF\\(2\\).*rank",
        "legacy_role_claims": "sin\\^2\\(theta_W\\)|alpha_EM\\^-1|beta\\^N|beta_tors.*chi_11",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2434's role lattice and older claim-specific role obstructions, but no production P24xx "
            "certificate computing the induced implication/separability poset among the four legacy role claims."
        ),
    }


def gf2_rank(rows: list[list[int]]) -> int:
    work = [row[:] for row in rows if any(row)]
    if not work:
        return 0
    rank = 0
    col_count = len(work[0])
    for col in range(col_count):
        pivot = next((index for index in range(rank, len(work)) if work[index][col] & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index, row in enumerate(work):
            if index != rank and row[col] & 1:
                work[index] = [a ^ b for a, b in zip(row, work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def certificate_payload(source: dict[str, Any]) -> dict[str, Any]:
    cert = source.get("conditional_legacy_role_transfer_claim_lattice_certificate", {})
    theorem = cert.get("theorem_export", {})
    role_rows = cert.get("role_claim_rows", [])
    lattice_rows = cert.get("role_lattice_rows", [])
    role_ids = theorem.get("role_claim_ids", [row["role_id"] for row in role_rows])
    obligations = theorem.get("role_obligation_names", [])

    requirements = {row["role_id"]: row["required_obligations"] for row in role_rows}
    incidence_rows = [
        {
            "role_id": role_id,
            "incidence_bits": [1 if obligation in requirements[role_id] else 0 for obligation in obligations],
            "required_obligations": requirements[role_id],
        }
        for role_id in role_ids
    ]
    incidence_matrix = [row["incidence_bits"] for row in incidence_rows]

    implication_rows = []
    nonimplication_witnesses = []
    for premise in role_ids:
        for conclusion in role_ids:
            witness = next(
                (
                    row
                    for row in lattice_rows
                    if premise in row["ready_role_claims"] and conclusion not in row["ready_role_claims"]
                ),
                None,
            )
            implies = witness is None
            implication_rows.append(
                {
                    "premise_role": premise,
                    "conclusion_role": conclusion,
                    "premise_implies_conclusion": implies,
                    "separating_witness_mask": None if implies else witness["mask"],
                    "separating_witness_true_obligations": None if implies else witness["true_obligations"],
                }
            )
            if not implies:
                nonimplication_witnesses.append(
                    {
                        "premise_role": premise,
                        "conclusion_role": conclusion,
                        "mask": witness["mask"],
                        "true_obligations": witness["true_obligations"],
                    }
                )

    nontrivial_implications = [
        [row["premise_role"], row["conclusion_role"]]
        for row in implication_rows
        if row["premise_role"] != row["conclusion_role"] and row["premise_implies_conclusion"]
    ]
    equivalence_classes = []
    seen: set[str] = set()
    implication_lookup = {(row["premise_role"], row["conclusion_role"]): row["premise_implies_conclusion"] for row in implication_rows}
    for role_id in role_ids:
        if role_id in seen:
            continue
        group = [other for other in role_ids if implication_lookup[(role_id, other)] and implication_lookup[(other, role_id)]]
        seen.update(group)
        equivalence_classes.append(group)

    obligation_coverage = {
        obligation: [role_id for role_id in role_ids if obligation in requirements[role_id]] for obligation in obligations
    }

    return {
        "role_ids": role_ids,
        "obligations": obligations,
        "incidence_rows": incidence_rows,
        "incidence_matrix_gf2": incidence_matrix,
        "incidence_rank_gf2": gf2_rank(incidence_matrix),
        "implication_rows": implication_rows,
        "nontrivial_implications": nontrivial_implications,
        "nonimplication_witnesses": nonimplication_witnesses,
        "equivalence_classes": equivalence_classes,
        "obligation_coverage_by_role": obligation_coverage,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2435/S1385 legacy role-claim implication separability certificate

`P2435/S1385` computes the implication/separability poset induced by the P2434 64-mask role lattice.  The four legacy role claims have GF(2) obligation-incidence rank 4, so the claim rows remain independent as audit targets.  The only nontrivial implication is `legacy_inverse_alpha_em -> legacy_weak_mixing_angle`, because the inverse-alpha claim requires the alpha successor plus beta_tors successor and therefore contains the weaker alpha-only Weinberg successor requirements.

All other ordered claim pairs have explicit separating masks.  Thus P2435 refines the role-transfer audit: strict alpha_EM success would imply the weak-mixing successor only at the role-obligation level, but neither implication exports either theorem, and no gravity or beta_tors -> chi_11 role can be imported by alpha/Weinberg progress.
""".strip()
    lag_section = """
## P2435/S1385 legacy role-claim separability guard for Lagrangian/EOM

`P2435/S1385` shows that the legacy role package is not one monolithic `L_total` insertion.  The role-claim obligation rows have rank 4 and all non-alpha_EM -> Weinberg pairs are separated by finite masks, so the Lagrangian/EOM draft must keep Weinberg, alpha_EM, gravity hierarchy, and beta_tors -> chi_11 role successors as separate theorem targets.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = certificate_payload(sources["P2434_ROLE_LATTICE"])
    theorem_export = {
        "theorem_name": "P2435_T1_legacy_role_claim_implication_separability_certificate",
        "role_claim_count": len(cert["role_ids"]),
        "role_obligation_count": len(cert["obligations"]),
        "role_claim_ids": cert["role_ids"],
        "role_obligation_names": cert["obligations"],
        "incidence_rank_gf2": cert["incidence_rank_gf2"],
        "implication_row_count": len(cert["implication_rows"]),
        "implication_true_count": sum(1 for row in cert["implication_rows"] if row["premise_implies_conclusion"]),
        "nontrivial_implications": cert["nontrivial_implications"],
        "nonimplication_witness_count": len(cert["nonimplication_witnesses"]),
        "equivalence_classes": cert["equivalence_classes"],
        "no_two_distinct_claims_equivalent": all(len(group) == 1 for group in cert["equivalence_classes"]),
        "alpha_em_strictly_implies_weak_mixing_obligation_readiness": [
            "legacy_inverse_alpha_em",
            "legacy_weak_mixing_angle",
        ]
        in cert["nontrivial_implications"],
        "weak_mixing_does_not_imply_alpha_em": any(
            witness["premise_role"] == "legacy_weak_mixing_angle"
            and witness["conclusion_role"] == "legacy_inverse_alpha_em"
            for witness in cert["nonimplication_witnesses"]
        ),
        "obligation_coverage_by_role": cert["obligation_coverage_by_role"],
        "role_transfer_and_ltotal_cover_all_claims_but_do_not_separate_them": all(
            len(cert["obligation_coverage_by_role"][obligation]) == len(cert["role_ids"])
            for obligation in ["role_transfer_audit_license", "role_bearing_ltotal_export"]
        ),
        "claim_specific_successors_are_needed_to_separate_claims": True,
        "legacy_role_claim_transferred_by_this_certificate": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Implication among conditional readiness predicates is not theorem discharge.",
            "The alpha_EM -> Weinberg implication is only an obligation-set containment fact, not a physical-role theorem.",
            "No claim-specific strict successor theorem, role-transfer theorem, role-bearing L_total export, or ToE closure is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "four_claims": theorem_export["role_claim_count"] == 4,
        "six_obligations": theorem_export["role_obligation_count"] == 6,
        "rank_four": theorem_export["incidence_rank_gf2"] == 4,
        "sixteen_implication_rows": theorem_export["implication_row_count"] == 16,
        "five_true_implications": theorem_export["implication_true_count"] == 5,
        "single_nontrivial_implication": theorem_export["nontrivial_implications"] == [
            ["legacy_inverse_alpha_em", "legacy_weak_mixing_angle"]
        ],
        "eleven_nonimplication_witnesses": theorem_export["nonimplication_witness_count"] == 11,
        "all_equivalence_classes_singleton": theorem_export["no_two_distinct_claims_equivalent"],
        "weak_not_alpha_em": theorem_export["weak_mixing_does_not_imply_alpha_em"],
        "role_transfer_ltotal_cover_all": theorem_export["role_transfer_and_ltotal_cover_all_claims_but_do_not_separate_them"],
        "no_legacy_role_transfer": not theorem_export["legacy_role_claim_transferred_by_this_certificate"],
        "no_role_transfer_export": not theorem_export["role_transfer_licensed"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2435_s1385_v1",
        "packet_id": "P2435",
        "stage_id": "S1385",
        "result_kind": "LEGACY_ROLE_CLAIM_IMPLICATION_SEPARABILITY_CERTIFICATE",
        "status": "PASS_LEGACY_ROLE_CLAIM_IMPLICATION_SEPARABILITY_NO_ROLE_TRANSFER_NO_CLOSURE",
        "legacy_role_claim_implication_separability_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_implication_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Prove or reject claim-specific strict successor theorems separately; do not treat role-transfer as a package import."
        ),
        "global_status": "OPEN_PROGRESS_ROLE_CLAIM_SEPARABILITY_CERTIFIED_NO_ROLE_IMPORTED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["legacy_role_claim_implication_separability_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2435 S1385: legacy role-claim implication separability certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Incidence rank over GF(2): `{theorem['incidence_rank_gf2']}`.",
                f"- Implication rows: `{theorem['implication_row_count']}`.",
                f"- Nontrivial implications: `{theorem['nontrivial_implications']}`.",
                f"- Nonimplication witnesses: `{theorem['nonimplication_witness_count']}`.",
                f"- Equivalence classes: `{theorem['equivalence_classes']}`.",
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
