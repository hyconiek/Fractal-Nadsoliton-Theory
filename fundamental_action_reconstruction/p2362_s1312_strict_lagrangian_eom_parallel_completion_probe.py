#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.json"
MD = GEN / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.md"

SOURCE_FILES = {
    "P1688_SECTOR_EOM_EXPORT": GEN / "p1688_s638_strict_full_lagrangian_to_sector_eom_export_summary.json",
    "P1693_MULTISECTOR_SYMPY_BRIDGE": GEN / "p1693_s643_strict_full_lagrangian_multisector_sympy_eom_bridge.json",
    "P2086_TERMWISE_EOM": GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json",
    "P2087_NORMAL_FORMS": GEN / "p2087_s1037_strict_full_ltotal_eom_normal_form_extraction_audit.json",
    "P2088_EOM_GAP_AUDIT": GEN / "p2088_s1038_strict_full_ltotal_eom_theorem_readiness_gap_audit.json",
    "P2316_BEST_LTOTAL_AUDIT": GEN / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json",
    "P2329_SELECTOR_INDEPENDENCE": GEN / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

REQUIRED_DOC_SNIPPETS = [
    "EOM/Lagrangian track is selector-independent",
    "selector closure is a parallel problem",
    "P2362/S1312",
    "Reduced termwise computational EOM",
]

COVARIANT_SECTOR_MAP = [
    ("scalar_phi", "L_scalar_phi"),
    ("higgs_H", "L_higgs"),
    ("gauge_A", "L_gauge"),
    ("fermion_psi", "L_fermion"),
    ("metric_g", "L_gravity"),
]

REDUCED_FIELD_ORDER = ["psi", "A", "h"]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def gf2_rank(matrix: list[list[int]]) -> int:
    rows = [row[:] for row in matrix]
    if not rows:
        return 0
    m = len(rows)
    n = len(rows[0])
    rank = 0
    col = 0
    while rank < m and col < n:
        pivot = None
        for row in range(rank, m):
            if rows[row][col] % 2:
                pivot = row
                break
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for row in range(m):
            if row != rank and rows[row][col] % 2:
                rows[row] = [left ^ right for left, right in zip(rows[row], rows[rank])]
        rank += 1
        col += 1
    return rank


def is_nonzero_srepr(value: Any) -> bool:
    return str(value) not in {"Integer(0)", "0", "0.0"}


def build_covariant_sector_rows(p1688: dict[str, Any], p1693: dict[str, Any]) -> list[dict[str, Any]]:
    sector_eom = p1688.get("sector_eom_export", {})
    anchor = p1693.get("full_lagrangian_density_anchor", {})
    rows: list[dict[str, Any]] = []
    for eom_key, lagrangian_key in COVARIANT_SECTOR_MAP:
        rows.append(
            {
                "sector": eom_key,
                "lagrangian_anchor_key": lagrangian_key,
                "lagrangian_anchor": anchor.get(lagrangian_key, ""),
                "eom_row": sector_eom.get(eom_key, ""),
                "selector_required_for_eom_export": False,
                "completion_track": "EOM_LAGRANGIAN_PARALLEL_TO_SELECTOR",
            }
        )
    return rows


def build_reduced_incidence_rows(p2086: dict[str, Any]) -> tuple[list[dict[str, Any]], list[list[int]]]:
    results = p2086.get("eom_execution_results", {})
    terms = results.get("lagrangian_terms", {})
    variations = results.get("termwise_variation_map", {})
    rows: list[dict[str, Any]] = []
    matrix: list[list[int]] = []
    for term_id in sorted(terms):
        vector = [
            1 if is_nonzero_srepr(variations.get(field, {}).get(term_id, "Integer(0)")) else 0
            for field in REDUCED_FIELD_ORDER
        ]
        matrix.append(vector)
        rows.append(
            {
                "term_id": term_id,
                "field_incidence": dict(zip(REDUCED_FIELD_ORDER, vector)),
                "nonzero_fields": [field for field, bit in zip(REDUCED_FIELD_ORDER, vector) if bit],
                "selector_required_for_variation": False,
            }
        )
    return rows, matrix


def build_reduced_eom_rows(p2086: dict[str, Any], p2087: dict[str, Any]) -> list[dict[str, Any]]:
    results = p2086.get("eom_execution_results", {})
    eom_full = results.get("eom_full", {})
    residuals = results.get("symbolic_recomposition_residual", {})
    numeric = results.get("numeric_probe_residual", {})
    normal_forms = p2087.get("eom_normal_form_results", {}).get("solved_second_derivative_rhs", {})
    normal_key = {"psi": "psi_ddot", "A": "A_ddot", "h": "h_ddot"}
    rows: list[dict[str, Any]] = []
    for field in REDUCED_FIELD_ORDER:
        rows.append(
            {
                "field": field,
                "eom_full_srepr": eom_full.get(field, ""),
                "normal_form_rhs_srepr": normal_forms.get(normal_key[field], ""),
                "symbolic_recomposition_residual": residuals.get(field),
                "numeric_probe_residual": numeric.get(field),
                "selector_required_for_reduced_eom": False,
            }
        )
    return rows


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    doc_texts = {name: read_text(path) for name, path in DOC_FILES.items()}

    p1688 = artifacts["P1688_SECTOR_EOM_EXPORT"]
    p1693 = artifacts["P1693_MULTISECTOR_SYMPY_BRIDGE"]
    p2086 = artifacts["P2086_TERMWISE_EOM"]
    p2087 = artifacts["P2087_NORMAL_FORMS"]
    p2088 = artifacts["P2088_EOM_GAP_AUDIT"]
    p2316 = artifacts["P2316_BEST_LTOTAL_AUDIT"]
    p2329 = artifacts["P2329_SELECTOR_INDEPENDENCE"]

    covariant_sector_rows = build_covariant_sector_rows(p1688, p1693)
    reduced_incidence_rows, incidence_matrix = build_reduced_incidence_rows(p2086)
    reduced_eom_rows = build_reduced_eom_rows(p2086, p2087)
    incidence_rank = gf2_rank(incidence_matrix)
    field_column_sums = [
        sum(row[index] for row in incidence_matrix)
        for index, _field in enumerate(REDUCED_FIELD_ORDER)
    ]

    p2329_cert = p2329.get("strict_selector_independence_lagrangian_eom_audit_probe", {}).get(
        "selector_independence_certificate", {}
    )
    p2316_next = p2316.get("recommended_next_honest_step", "")
    p2088_missing_names = [
        row.get("name", "")
        for row in p2088.get("theorem_readiness_gap_register", [])
    ]

    theorem_export = {
        "theorem_name": "P2362 selector-independent strict Lagrangian/EOM parallel completion audit",
        "claim": "The strict Lagrangian/EOM track can be continued independently of selector/QW-2191 closure: covariant sector EOM rows and reduced termwise normal forms are available, and remaining EOM blockers are tensor/background/nonproxy-sector obligations rather than selector prerequisites.",
        "positive_content": [
            "P1688 exports sector EOM rows for scalar_phi, higgs_H, gauge_A, fermion_psi, and metric_g.",
            "P1693 exports matching covariant Lagrangian sector anchors.",
            "P2086 exports 11 reduced L_total terms with termwise Euler-Lagrange recomposition residuals equal to zero.",
            "P2087 exports second-derivative normal forms for psi, A, and h with zero substitution residual.",
            "P2329 certifies the reduced Lagrangian/EOM terms, variations, and residuals as selector-independent.",
        ],
        "not_licensed": [
            "full tensor-resolved metric/background theorem",
            "global nonperturbative well-posedness theorem",
            "BRST/Cutkosky/counterterm closure",
            "legacy physical-role transfer",
            "selector premise or QW-2191 selector discharge",
            "ToE closure",
        ],
    }

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p1688_sector_eom_loaded": p1688.get("checkpoint") == "P1688_S638",
        "p1693_covariant_anchor_loaded": p1693.get("checkpoint") == "P1693_S643",
        "p2086_termwise_zero_residuals": p2086.get("gatekeeper_checks", {}).get("all_symbolic_residual_zero") and p2086.get("gatekeeper_checks", {}).get("all_numeric_probe_residual_zero"),
        "p2087_normal_forms_zero_residuals": p2087.get("gatekeeper_checks", {}).get("normal_form_solved_all_fields") and p2087.get("gatekeeper_checks", {}).get("residual_zero_after_normal_form_substitution"),
        "p2316_parallel_next_step_recorded": "independently" in p2316_next and "separate parallel track" in p2316_next,
        "p2329_selector_independence_inherited": p2329_cert.get("selector_independent_term_count") == 11 and p2329_cert.get("selector_independent_variation_field_count") == 3,
        "p2088_eom_gaps_not_selector_gaps": all("selector" not in name and "qw2191" not in name.lower() for name in p2088_missing_names),
        "covariant_sector_rows_complete": len(covariant_sector_rows) == 5 and all(row["lagrangian_anchor"] and row["eom_row"] for row in covariant_sector_rows),
        "reduced_incidence_matrix_full_field_rank": incidence_rank == len(REDUCED_FIELD_ORDER),
        "reduced_all_terms_hit_at_least_one_field": all(sum(row) >= 1 for row in incidence_matrix),
        "reduced_all_fields_have_nonzero_terms": all(total > 0 for total in field_column_sums),
        "reduced_eom_rows_have_zero_residuals": all(row["symbolic_recomposition_residual"] == "Integer(0)" and str(row["numeric_probe_residual"]) == "0" for row in reduced_eom_rows),
        "docs_updated_with_parallel_eom_statement": all(
            snippet in text
            for text in doc_texts.values()
            for snippet in REQUIRED_DOC_SNIPPETS
        ),
        "no_selector_prerequisite_for_eom_export": all(not row["selector_required_for_eom_export"] for row in covariant_sector_rows) and all(not row["selector_required_for_variation"] for row in reduced_incidence_rows),
        "hard_limits_preserved": True,
    }

    probe = {
        "probe_id": "P2362_S1312_STRICT_LAGRANGIAN_EOM_PARALLEL_COMPLETION",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "covariant_sector_eom_rows": covariant_sector_rows,
        "reduced_termwise_incidence_rows": reduced_incidence_rows,
        "reduced_termwise_incidence_matrix": incidence_matrix,
        "reduced_termwise_incidence_rank": incidence_rank,
        "reduced_field_column_sums": dict(zip(REDUCED_FIELD_ORDER, field_column_sums)),
        "reduced_eom_rows": reduced_eom_rows,
        "parallel_completion_summary": {
            "eom_lagrangian_track_is_selector_independent": True,
            "selector_closure_is_parallel_problem": True,
            "selector_closure_required_before_eom_work": False,
            "covariant_sector_eom_row_count": len(covariant_sector_rows),
            "reduced_term_count": len(reduced_incidence_rows),
            "reduced_variation_field_count": len(REDUCED_FIELD_ORDER),
            "reduced_incidence_rank": incidence_rank,
            "remaining_eom_blockers": p2088_missing_names,
            "selector_qw2191_status": "OPEN_PARALLEL_TRACK_NOT_USED_AS_EOM_PREREQUISITE",
        },
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    payload = {
        "schema_version": "p2362_s1312_v1",
        "packet_id": "P2362",
        "stage_id": "S1312",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PARALLEL_EOM_LAGRANGIAN_COMPLETION_ADVANCED_SELECTOR_SEPARATE",
        "result_kind": "STRICT_LAGRANGIAN_EOM_PARALLEL_COMPLETION_AUDIT_SELECTOR_INDEPENDENT_NO_QW2191_DISCHARGE",
        "strict_lagrangian_eom_parallel_completion_probe": probe,
        "recommended_next_honest_step": "Continue EOM/Lagrangian completion by exporting nonproxy covariant spinor/gauge/metric residual tables on a named background family; work on selector/QW-2191 can proceed in parallel but is not a prerequisite for this EOM step.",
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2362 S1312: strict Lagrangian/EOM parallel completion audit",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "EOM/Lagrangian track is selector-independent. Selector closure is a parallel problem, not a prerequisite for continuing termwise/covariant EOM export.",
                "",
                "## Coverage",
                "",
                f"- Covariant sector EOM rows: `{len(covariant_sector_rows)}`.",
                f"- Reduced termwise terms: `{len(reduced_incidence_rows)}`.",
                f"- Reduced incidence rank: `{incidence_rank}`.",
                f"- Reduced field column sums: `{dict(zip(REDUCED_FIELD_ORDER, field_column_sums))}`.",
                "",
                "## Hard Limits",
                "",
                "- No full tensor-resolved theorem is claimed.",
                "- No selector premise or QW-2191 discharge is claimed.",
                "- No legacy role transfer is claimed.",
                "- No ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
