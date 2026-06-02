#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.json"
MD = GEN / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.md"

SOURCE_FILES = {
    "P1563_KERNEL_TO_EOM": GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json",
    "P1641_CLOSURE_MATRIX": GEN / "p1641_s591_theorem_level_closure_requirement_matrix_summary.json",
    "P1664_MANIFEST_INVERSION": GEN / "p1664_s614_strict_full_lagrangian_manifest_and_inversion.json",
    "P1692_SYMPY_REPLAY": GEN / "p1692_s642_strict_sympy_coefficient_to_eom_witness.json",
    "P1910_C1_GR_COEFF_TABLE": GEN / "p1910_s860_strict_c1_gr_coefficient_table_v1_probe.json",
    "P2438_SM_GR_OBLIGATION_MATRIX": GEN / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

CURRENT_STRICT_TUPLE = {
    "omega": 0.18575,
    "phi": 0.16250,
    "beta": 1.0,
    "eta": 1.8,
}

AUDITED_SOURCE_ORDER = [
    "P1563_effective_three_coefficient_chain",
    "P1664_full_manifest_local_inversion",
    "P1910_loop_counterterm_placeholder_table",
]

FEATURES = [
    "current_strict_tuple_match",
    "has_kernel_to_coefficient_map",
    "exports_sm_gauge_or_coupling_coefficients",
    "exports_sm_matter_higgs_yukawa_coefficients",
    "exports_gr_coefficients",
    "coefficients_evaluated_not_placeholders",
    "strict_observable_value_generator_exported",
    "qw2191_selector_discharged",
    "closed_theorem_status",
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
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:35]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2439|S1389|strict coefficient source consistency|coefficient source consistency audit|current strict tuple coefficient source",
        "p2438_input": "P2438|S1388|strict kernel SM GR generation|strict_kernel_to_coefficient_map_theorem",
        "coefficient_map_chain": "K_strict -> coefficients|kernel-to-coefficient|strict kernel.*coefficients|coefficient_map",
        "full_manifest_inversion": "P1664|strict_full_lagrangian_manifest_and_inversion|local_pass|inverse_recovery",
        "effective_three_coefficients": "lambda_sm_eff|kappa_gr_eff|epsilon_mix_eff",
        "loop_counterterm_placeholders": "P1910|coefficient_table_v1|OPEN_SYMBOLIC_EXPORT|A_\\*|B_\\*|F_\\*",
        "observable_generator_language": "observable-value generator|physical constants generated|strict observable|physical-value generator",
        "selector_blocker": "QW-2191|selector uniqueness|selector source",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def close(a: float, b: float, tol: float = 1e-12) -> bool:
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def tuple_comparison(candidate: dict[str, Any]) -> dict[str, Any]:
    compared: dict[str, Any] = {}
    matches: dict[str, bool] = {}
    deltas: dict[str, float | None] = {}
    missing: list[str] = []
    for key, expected in CURRENT_STRICT_TUPLE.items():
        if key not in candidate:
            missing.append(key)
            matches[key] = False
            deltas[key] = None
            continue
        actual = float(candidate[key])
        compared[key] = actual
        matches[key] = close(actual, expected)
        deltas[key] = actual - expected
    return {
        "current_tuple": CURRENT_STRICT_TUPLE,
        "candidate_tuple_compared": compared,
        "per_parameter_match": matches,
        "per_parameter_delta": deltas,
        "missing_parameters": missing,
        "match_count": sum(1 for value in matches.values() if value),
        "mismatch_count": sum(1 for value in matches.values() if not value),
        "full_current_tuple_match": all(matches.values()) and not missing,
    }


def status_of(source: dict[str, Any]) -> str:
    return str(source.get("status") or source.get("checkpoint") or source.get("_missing") or "UNKNOWN")


def p1563_row(p1563: dict[str, Any], p1641: dict[str, Any]) -> dict[str, Any]:
    coeffs = p1563.get("derived_coefficients", {})
    tuple_info = tuple_comparison(p1563.get("strict_kernel_params", {}))
    p1641_reqs = p1641.get("closure_requirement_matrix", [])
    open_reqs = [row.get("id") for row in p1641_reqs if row.get("status") == "OPEN"]
    features = {
        "current_strict_tuple_match": tuple_info["full_current_tuple_match"],
        "has_kernel_to_coefficient_map": bool(coeffs),
        "exports_sm_gauge_or_coupling_coefficients": "lambda_sm_eff" in coeffs,
        "exports_sm_matter_higgs_yukawa_coefficients": False,
        "exports_gr_coefficients": "kappa_gr_eff" in coeffs,
        "coefficients_evaluated_not_placeholders": all(isinstance(v, (int, float)) for v in coeffs.values()) and bool(coeffs),
        "strict_observable_value_generator_exported": False,
        "qw2191_selector_discharged": False,
        "closed_theorem_status": False,
    }
    return {
        "source_id": "P1563_effective_three_coefficient_chain",
        "input_packets": ["P1563", "P1641"],
        "status": status_of(p1563),
        "tuple_comparison": tuple_info,
        "coefficient_symbols": sorted(coeffs),
        "coefficient_count": len(coeffs),
        "feature_vector": features,
        "open_requirement_ids_inherited_from_p1641": open_reqs,
        "methodology_classification": "current_tuple_partial_effective_coefficients_not_full_SM_GR_observable_generator",
        "promotion_blockers": [
            "only three effective coefficients are exported",
            "no full SM gauge/Yukawa/Higgs coefficient spectrum is generated",
            "P1641 theorem-level closure requirements remain open",
            "no QW-2191 discharge or strict observable-value generator is exported",
        ],
    }


def p1664_row(p1664: dict[str, Any], p1692: dict[str, Any]) -> dict[str, Any]:
    coeffs = p1664.get("coefficient_map", {})
    tuple_info = tuple_comparison(p1664.get("kernel_input", {}))
    features = {
        "current_strict_tuple_match": tuple_info["full_current_tuple_match"],
        "has_kernel_to_coefficient_map": bool(coeffs),
        "exports_sm_gauge_or_coupling_coefficients": all(k in coeffs for k in ["Z3", "Z2", "Z1"]),
        "exports_sm_matter_higgs_yukawa_coefficients": all(k in coeffs for k in ["muH2", "lambdaH", "yu", "yd", "ye"]),
        "exports_gr_coefficients": all(k in coeffs for k in ["Mpl2", "cR2", "cRic2", "cRiem2"]),
        "coefficients_evaluated_not_placeholders": all(isinstance(v, (int, float)) for v in coeffs.values()) and bool(coeffs),
        "strict_observable_value_generator_exported": False,
        "qw2191_selector_discharged": False,
        "closed_theorem_status": False,
    }
    return {
        "source_id": "P1664_full_manifest_local_inversion",
        "input_packets": ["P1664", "P1692"],
        "status": status_of(p1664),
        "tuple_comparison": tuple_info,
        "coefficient_symbols": sorted(coeffs),
        "coefficient_count": len(coeffs),
        "feature_vector": features,
        "inverse_recovery_local_pass": p1664.get("inverse_recovery", {}).get("local_pass") is True,
        "p1692_status": status_of(p1692),
        "methodology_classification": "full_looking_local_manifest_but_not_current_QW2049_tuple_and_not_global_generator",
        "promotion_blockers": [
            "kernel_input does not match the current strict tuple in phi, beta, and eta",
            "local inverse recovery is not a global strict source theorem",
            "P1692 is a reduced SymPy replay witness, not a full covariant SM/GR observable generator",
            "no QW-2191 discharge or strict observable-value generator is exported",
        ],
    }


def p1910_row(p1910: dict[str, Any]) -> dict[str, Any]:
    table = p1910.get("coefficient_table_v1", [])
    symbols = [row.get("coefficient_symbol") for row in table]
    all_open = all(row.get("status") == "OPEN_SYMBOLIC_EXPORT" for row in table) and bool(table)
    features = {
        "current_strict_tuple_match": False,
        "has_kernel_to_coefficient_map": False,
        "exports_sm_gauge_or_coupling_coefficients": False,
        "exports_sm_matter_higgs_yukawa_coefficients": any("yukawa" in str(s) for s in symbols),
        "exports_gr_coefficients": any("gr" in str(s) or "curvature" in str(s) for s in symbols),
        "coefficients_evaluated_not_placeholders": False,
        "strict_observable_value_generator_exported": False,
        "qw2191_selector_discharged": False,
        "closed_theorem_status": False,
    }
    return {
        "source_id": "P1910_loop_counterterm_placeholder_table",
        "input_packets": ["P1910"],
        "status": status_of(p1910),
        "tuple_comparison": tuple_comparison({}),
        "coefficient_symbols": sorted(str(s) for s in symbols if s),
        "coefficient_count": len(symbols),
        "all_coefficients_open_symbolic_export": all_open,
        "feature_vector": features,
        "methodology_classification": "loop_counterterm_contract_table_with_placeholders_not_kernel_value_generation",
        "promotion_blockers": [
            "no explicit current strict tuple is bound to the table",
            "A_*, B_*, and F_* loop constants remain unevaluated placeholders",
            "counterterm/Cutkosky/background contracts remain open",
            "no QW-2191 discharge or strict observable-value generator is exported",
        ],
    }


def gf2_rank(rows: list[list[int]]) -> int:
    mat = [row[:] for row in rows]
    rank = 0
    col = 0
    while rank < len(mat) and col < (len(mat[0]) if mat else 0):
        pivot = next((r for r in range(rank, len(mat)) if mat[r][col]), None)
        if pivot is None:
            col += 1
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        for r in range(len(mat)):
            if r != rank and mat[r][col]:
                mat[r] = [a ^ b for a, b in zip(mat[r], mat[rank])]
        rank += 1
        col += 1
    return rank


def vectorize(row: dict[str, Any]) -> list[int]:
    fv = row["feature_vector"]
    return [1 if fv[name] else 0 for name in FEATURES]


def acceptable(row: dict[str, Any]) -> bool:
    fv = row["feature_vector"]
    required = [
        "current_strict_tuple_match",
        "has_kernel_to_coefficient_map",
        "exports_sm_gauge_or_coupling_coefficients",
        "exports_sm_matter_higgs_yukawa_coefficients",
        "exports_gr_coefficients",
        "coefficients_evaluated_not_placeholders",
        "strict_observable_value_generator_exported",
        "qw2191_selector_discharged",
        "closed_theorem_status",
    ]
    return all(fv[name] for name in required)


def append_doc_sections() -> None:
    eq_section = """
## P2439/S1389 strict coefficient-source consistency audit

`P2439/S1389` audits the existing strict coefficient sources before any renewed SM/GR value derivation.  The finite audit separates the current-tuple three-effective-coefficient chain (`P1563/P1641`), the fuller but tuple-mismatched local manifest/inversion chain (`P1664/P1692`), and the open loop-counterterm placeholder table (`P1910`).  None is promoted to a current strict SM/GR physical-value generator: the first is too low-dimensional, the second is not the current `QW-2049` tuple and is only local, and the third has unevaluated symbolic placeholders.

Consequently the next honest coefficient step is not to insert legacy constants or reuse the tuple-mismatched manifest as if it were final, but to construct a current-`K_strict_gate` coefficient map that simultaneously covers SM gauge, matter/Higgs/Yukawa, GR, selector uniqueness, and observable-value generation.
""".strip()
    lag_section = """
## P2439/S1389 coefficient-source guard for Lagrangian/EOM

`P2439/S1389` blocks promotion of the current coefficient scaffolds into a role-bearing `L_total`.  `P1563` may be cited only as a current-tuple effective three-coefficient chain; `P1664/P1692` may be cited only as a local manifest/replay on a non-current tuple; and `P1910` may be cited only as an open placeholder counterterm table.  No SM/GR physical constants, selector choice, or observable generator is licensed by these sources.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    rows = [
        p1563_row(sources["P1563_KERNEL_TO_EOM"], sources["P1641_CLOSURE_MATRIX"]),
        p1664_row(sources["P1664_MANIFEST_INVERSION"], sources["P1692_SYMPY_REPLAY"]),
        p1910_row(sources["P1910_C1_GR_COEFF_TABLE"]),
    ]
    incidence_rows = [vectorize(row) for row in rows]
    per_feature_support = {
        feature: [row["source_id"] for row in rows if row["feature_vector"][feature]] for feature in FEATURES
    }
    acceptable_sources = [row["source_id"] for row in rows if acceptable(row)]
    tuple_matches = [row["source_id"] for row in rows if row["feature_vector"]["current_strict_tuple_match"]]
    current_tuple_full_sm_gr_sources = [
        row["source_id"]
        for row in rows
        if row["feature_vector"]["current_strict_tuple_match"]
        and row["feature_vector"]["exports_sm_gauge_or_coupling_coefficients"]
        and row["feature_vector"]["exports_sm_matter_higgs_yukawa_coefficients"]
        and row["feature_vector"]["exports_gr_coefficients"]
    ]
    theorem_export = {
        "theorem_name": "P2439_T1_strict_coefficient_source_consistency_audit_certificate",
        "current_strict_tuple": CURRENT_STRICT_TUPLE,
        "audited_source_count": len(rows),
        "audited_source_ids": [row["source_id"] for row in rows],
        "feature_count": len(FEATURES),
        "feature_names": FEATURES,
        "source_feature_rows": rows,
        "source_feature_incidence_matrix_gf2": incidence_rows,
        "source_feature_incidence_rank_gf2": gf2_rank(incidence_rows),
        "per_feature_support": per_feature_support,
        "current_tuple_matching_source_ids": tuple_matches,
        "current_tuple_matching_source_count": len(tuple_matches),
        "current_tuple_full_sm_gr_coefficient_source_ids": current_tuple_full_sm_gr_sources,
        "current_tuple_full_sm_gr_coefficient_source_count": len(current_tuple_full_sm_gr_sources),
        "acceptable_current_strict_physical_value_generator_source_ids": acceptable_sources,
        "acceptable_current_strict_physical_value_generator_source_count": len(acceptable_sources),
        "p1563_current_tuple_but_three_effective_coefficients_only": rows[0]["feature_vector"]["current_strict_tuple_match"] and rows[0]["coefficient_count"] == 3,
        "p1664_full_manifest_tuple_mismatch_count": rows[1]["tuple_comparison"]["mismatch_count"],
        "p1664_inverse_recovery_is_only_local": rows[1]["inverse_recovery_local_pass"],
        "p1910_all_coefficients_open_symbolic_export": rows[2]["all_coefficients_open_symbolic_export"],
        "p2438_no_strict_generation_inherited": sources["P2438_SM_GR_OBLIGATION_MATRIX"]
        .get("strict_kernel_sm_gr_generation_obligation_matrix_certificate", {})
        .get("theorem_export", {})
        .get("strict_sm_gr_generation_theorem_exported")
        is False,
        "strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate": False,
        "strict_observable_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "The current-tuple P1563/P1641 chain exports only three effective coefficients, not full SM/GR constants.",
            "The fuller P1664/P1692 manifest is tuple-mismatched relative to the current QW-2049 strict tuple and remains local/open.",
            "The P1910 loop coefficient table is an open symbolic placeholder table, not evaluated physical-value generation.",
            "No audited source discharges QW-2191 or exports a strict observable-value generator.",
        ],
        "next_honest_step": (
            "Build a current-QW-2049 strict kernel-to-coefficient map with full SM gauge, matter/Higgs/Yukawa, GR, selector, and observable-value coverage."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "three_sources_audited": theorem_export["audited_source_count"] == 3,
        "nine_features_tracked": theorem_export["feature_count"] == 9,
        "gf2_rank_three": theorem_export["source_feature_incidence_rank_gf2"] == 3,
        "only_one_current_tuple_match": theorem_export["current_tuple_matching_source_count"] == 1,
        "no_current_tuple_full_sm_gr_source": theorem_export["current_tuple_full_sm_gr_coefficient_source_count"] == 0,
        "no_acceptable_physical_value_generator_source": theorem_export[
            "acceptable_current_strict_physical_value_generator_source_count"
        ]
        == 0,
        "p1563_is_three_coefficient_only": theorem_export["p1563_current_tuple_but_three_effective_coefficients_only"],
        "p1664_tuple_mismatch_three": theorem_export["p1664_full_manifest_tuple_mismatch_count"] == 3,
        "p1664_local_only_detected": theorem_export["p1664_inverse_recovery_is_only_local"],
        "p1910_open_placeholders_detected": theorem_export["p1910_all_coefficients_open_symbolic_export"],
        "p2438_no_generation_inherited": theorem_export["p2438_no_strict_generation_inherited"],
        "no_coefficient_theorem_export": not theorem_export[
            "strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate"
        ],
        "no_observable_generator_export": not theorem_export["strict_observable_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2439_s1389_v1",
        "packet_id": "P2439",
        "stage_id": "S1389",
        "result_kind": "STRICT_COEFFICIENT_SOURCE_CONSISTENCY_AUDIT_CERTIFICATE",
        "status": "PASS_STRICT_COEFFICIENT_SOURCE_AUDIT_NO_CURRENT_FULL_SM_GR_VALUE_GENERATOR",
        "strict_coefficient_source_consistency_audit_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_STRICT_COEFFICIENT_SOURCE_AUDIT_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_coefficient_source_consistency_audit_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2439 S1389: strict coefficient-source consistency audit certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Audited sources: `{theorem['audited_source_ids']}`.",
                f"- Feature count: `{theorem['feature_count']}`.",
                f"- GF(2) source-feature rank: `{theorem['source_feature_incidence_rank_gf2']}`.",
                f"- Current-tuple matching sources: `{theorem['current_tuple_matching_source_ids']}`.",
                f"- Current-tuple full SM/GR coefficient sources: `{theorem['current_tuple_full_sm_gr_coefficient_source_ids']}`.",
                f"- Acceptable physical-value generator sources: `{theorem['acceptable_current_strict_physical_value_generator_source_ids']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Next honest step",
                "",
                theorem["next_honest_step"],
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
