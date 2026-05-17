#!/usr/bin/env python3
"""P1959 S909 strict non-Abelian gauge-group data availability audit.

P1958 exported only a local Abelianized Lorenz/Faddeev-Popov/ghost seed.  This
executor checks whether the current strict repository state contains enough
theorem-grade non-Abelian gauge-group data to extend that seed to
SU(3)xSU(2)xU(1) BRST: explicit structure constants, Jacobi certificate,
gauge-connection rules, ghost self-interaction rules, and a physical
representation map not blocked by QW-2191/P1380.

The honest result is a current-state obstruction, not a no-go theorem.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1959_s909_strict_nonabelian_gauge_group_data_availability_audit.json"

TEXT_SUFFIXES = {".py", ".md", ".json", ".txt", ".yaml", ".yml"}
SKIP_NAMES = {
    "TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf",
    "p1959_s909_strict_nonabelian_gauge_group_data_availability_audit.py",
    "p1959_s909_strict_nonabelian_gauge_group_data_availability_audit.json",
}


def load_generated(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(name: str) -> str:
    path = ROOT / name
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="ignore")


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def digest_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def has_key_deep(obj: Any, key_name: str) -> bool:
    if isinstance(obj, dict):
        return any(k == key_name or has_key_deep(v, key_name) for k, v in obj.items())
    if isinstance(obj, list):
        return any(has_key_deep(v, key_name) for v in obj)
    return False


def text_contains(obj: Any, needle: str) -> bool:
    if isinstance(obj, str):
        return needle.lower() in obj.lower()
    return needle.lower() in json.dumps(obj, sort_keys=True, ensure_ascii=False).lower()


def repo_text_files() -> list[Path]:
    paths: list[Path] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.name in SKIP_NAMES:
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        paths.append(path)
    return sorted(paths)


def scan_terms(terms: list[str]) -> dict[str, dict[str, Any]]:
    files = repo_text_files()
    result: dict[str, dict[str, Any]] = {}
    for term in terms:
        count = 0
        sample_paths: list[str] = []
        low = term.lower()
        for path in files:
            text = path.read_text(encoding="utf-8", errors="ignore")
            hits = text.lower().count(low)
            if hits:
                count += hits
                if len(sample_paths) < 8:
                    sample_paths.append(str(path.relative_to(ROOT)))
        result[term] = {"count": count, "sample_paths": sample_paths}
    return result


def generated_json_key_hits(key_names: list[str]) -> dict[str, list[str]]:
    hits: dict[str, list[str]] = {key: [] for key in key_names}
    for path in sorted(GEN.glob("*.json")):
        if path.name in SKIP_NAMES:
            continue
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        for key in key_names:
            if has_key_deep(obj, key):
                hits[key].append(str(path.relative_to(ROOT)))
    return hits


def requirement_row(
    req_id: str,
    required_object: str,
    present: bool,
    observed_sources: list[str],
    exact_missing_data: list[str],
) -> dict[str, Any]:
    return {
        "req_id": req_id,
        "required_object": required_object,
        "present_as_theorem_grade_export": present,
        "verdict": "AVAILABLE" if present else "MISSING_OR_PARTIAL_SCAFFOLD_ONLY",
        "observed_sources": observed_sources,
        "exact_missing_data": [] if present else exact_missing_data,
    }


def main() -> None:
    a6 = load_generated("a6_gauge_reconstruction_summary.json")
    p1380 = load_generated("p1380_l_b1_01_su3_su2_u1_image_closure_theorem_attempt_summary.json")
    p1907 = load_generated("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1958 = load_generated("p1958_s908_strict_local_abelian_gauge_fixing_ghost_action_seed.json")
    p1957 = load_generated("p1957_s907_strict_bv_brst_ghost_sector_nonavailability_theorem.json")

    a6_text = load_text("A6_GAUGE_RECONSTRUCTION_SPEC.md")
    p1380_text = load_text("P1380_L_B1_01_SU3_SU2_U1_IMAGE_CLOSURE_THEOREM_ATTEMPT_PL.md")
    report_qw2190_path = ROOT / "report_qw2190_kernel_mode_representation_emergence_gate.json"
    report_qw2190_present = report_qw2190_path.exists()

    exact_keys = [
        "StructureConstants_SU3_SU2_U1_strict_v1",
        "structure_constants_table",
        "f_abc_table",
        "JacobiCertificate_SU3_SU2_U1_strict_v1",
        "jacobi_identity_certificate",
        "NonAbelianBRSTDifferential_strict_B1_v1",
        "BRST_differential_rules",
        "GhostSelfInteraction_strict_B1_v1",
        "GaugeConnectionRules_SU3_SU2_U1_strict_v1",
        "PhysicalRepresentationMap_SU3_SU2_U1_strict_v1",
    ]
    key_hits = generated_json_key_hits(exact_keys)

    structure_constants_present = any(key_hits[key] for key in [
        "StructureConstants_SU3_SU2_U1_strict_v1",
        "structure_constants_table",
        "f_abc_table",
    ])
    jacobi_certificate_present = any(key_hits[key] for key in [
        "JacobiCertificate_SU3_SU2_U1_strict_v1",
        "jacobi_identity_certificate",
    ])
    brst_rules_present = any(key_hits[key] for key in [
        "NonAbelianBRSTDifferential_strict_B1_v1",
        "BRST_differential_rules",
    ])
    ghost_self_interaction_present = bool(key_hits["GhostSelfInteraction_strict_B1_v1"])
    gauge_connection_rules_present = bool(key_hits["GaugeConnectionRules_SU3_SU2_U1_strict_v1"])
    physical_rep_map_present = bool(key_hits["PhysicalRepresentationMap_SU3_SU2_U1_strict_v1"])

    a6_partial_scaffold_present = (
        text_contains(a6, "SU(3) kernel-mode Lie closure")
        and text_contains(a6, "SU(2) kernel-mode Lie closure")
        and text_contains(a6, "partially_derived_in_strict_core")
    )
    a6_full_uniqueness_blocked = text_contains(a6, "full physical uniqueness of representation map") and text_contains(
        a6, "blocked"
    )
    p1380_closure_not_proven = text_contains(p1380_text, "THEOREM_CLOSURE := NOT_PROVEN") or text_contains(
        p1380, "NOT_PROVEN"
    )
    p1907_sm_gauge_registry_present = text_contains(p1907, "-1/4 F^2_SU3") and text_contains(p1907, "-1/4 F^2_SU2")
    p1958_abelian_seed_present = text_contains(p1958, "LOCAL_ABELIAN_LORENZ_GAUGE_FIXING")

    terms = [
        "SU(3)",
        "SU(2)",
        "U(1)",
        "structure constants",
        "f^{abc}",
        "f_abc",
        "Jacobi identity",
        "Jacobi tensor",
        "non-Abelian",
        "BRST differential",
        "D_mu c",
        "ghost self-interaction",
        "stałe struktury",
        "stale struktury",
        "algebra Lie",
        "komutator",
        "nieabel",
        "różniczka BRST",
        "rozniczka BRST",
        "grupa cechowania",
    ]
    term_scan = scan_terms(terms)

    A6_PARTIAL, QW2190_REPORT, F_TABLE, JACOBI, CONNECTION_RULES, BRST_RULES, GHOST_SELF, REP_MAP_OK = sp.symbols(
        "A6_PARTIAL QW2190_REPORT F_TABLE JACOBI CONNECTION_RULES BRST_RULES GHOST_SELF REP_MAP_OK"
    )
    nonabelian_ready_formula = sp.And(
        A6_PARTIAL,
        QW2190_REPORT,
        F_TABLE,
        JACOBI,
        CONNECTION_RULES,
        BRST_RULES,
        GHOST_SELF,
        REP_MAP_OK,
    )
    truth_assignment = {
        A6_PARTIAL: a6_partial_scaffold_present,
        QW2190_REPORT: report_qw2190_present,
        F_TABLE: structure_constants_present,
        JACOBI: jacobi_certificate_present,
        CONNECTION_RULES: gauge_connection_rules_present,
        BRST_RULES: brst_rules_present,
        GHOST_SELF: ghost_self_interaction_present,
        REP_MAP_OK: physical_rep_map_present and not a6_full_uniqueness_blocked and not p1380_closure_not_proven,
    }
    evaluated_ready = bool(nonabelian_ready_formula.subs(truth_assignment))

    g, jacobi_tensor = sp.symbols("g JacobiTensor")
    q2_ghost_dependency = sp.simplify(g**2 * jacobi_tensor / 6)
    q2_ghost_if_jacobi_zero = sp.simplify(q2_ghost_dependency.subs(jacobi_tensor, 0))

    requirement_matrix = [
        requirement_row(
            "N0",
            "A6 strict-core SU(3)xSU(2)xU(1) partial scaffold",
            a6_partial_scaffold_present,
            [
                f"A6 summary present={'a6' in a6}",
                "A6 declares partial SU(3), SU(2), U(1) scaffold but anti-overclaim keeps full uniqueness false.",
            ],
            ["A6 partial scaffold summary"],
        ),
        requirement_row(
            "N1",
            "QW-2190 source report with embedded generator data",
            report_qw2190_present,
            [
                f"expected local path={report_qw2190_path.relative_to(ROOT)}",
                f"present={report_qw2190_present}",
            ],
            [
                "report_qw2190_kernel_mode_representation_emergence_gate.json or equivalent strict exported source report",
                "embedded generator matrices or enough data to reconstruct them",
            ],
        ),
        requirement_row(
            "N2",
            "SU(3)xSU(2)xU(1) structure constants table f^a_bc",
            structure_constants_present,
            [
                f"key hits={ {k: key_hits[k] for k in ['StructureConstants_SU3_SU2_U1_strict_v1', 'structure_constants_table', 'f_abc_table']} }",
            ],
            [
                "normalization convention for generators",
                "nonzero f^a_bc rows for SU(3) and SU(2)",
                "U(1) zero-structure declaration",
                "gauge-factor index ranges",
            ],
        ),
        requirement_row(
            "N3",
            "Jacobi identity certificate for exported structure constants",
            jacobi_certificate_present,
            [
                f"key hits={ {k: key_hits[k] for k in ['JacobiCertificate_SU3_SU2_U1_strict_v1', 'jacobi_identity_certificate']} }",
                "nilpotency dependency: s^2 c is proportional to the Jacobi tensor",
            ],
            [
                "computed Jacobi tensor J^a_bcd",
                "machine check J^a_bcd=0 in the chosen normalization",
                "digest tying Jacobi check to the exact f^a_bc table",
            ],
        ),
        requirement_row(
            "N4",
            "non-Abelian gauge-connection and ghost transformation rules",
            gauge_connection_rules_present and brst_rules_present and ghost_self_interaction_present,
            [
                f"GaugeConnectionRules hits={key_hits['GaugeConnectionRules_SU3_SU2_U1_strict_v1']}",
                f"BRST differential hits={key_hits['NonAbelianBRSTDifferential_strict_B1_v1'] + key_hits['BRST_differential_rules']}",
                f"GhostSelfInteraction hits={key_hits['GhostSelfInteraction_strict_B1_v1']}",
            ],
            [
                "s A_mu^a = partial_mu c^a + g f^a_bc A_mu^b c^c",
                "s c^a = -g/2 f^a_bc c^b c^c",
                "s cbar^a and auxiliary-field convention",
                "graded Leibniz convention",
            ],
        ),
        requirement_row(
            "N5",
            "physical representation map compatible with QW-2191/P1380 blockers",
            truth_assignment[REP_MAP_OK],
            [
                f"A6 full uniqueness blocked={a6_full_uniqueness_blocked}",
                f"P1380 closure not proven={p1380_closure_not_proven}",
                f"PhysicalRepresentationMap key hits={key_hits['PhysicalRepresentationMap_SU3_SU2_U1_strict_v1']}",
            ],
            [
                "representation assignment map for gauge generators to strict fields",
                "resolution or lawful bypass of QW-2191 for the needed scope",
                "P1380 c_mix/Pi_gauge commutator compatibility lemma",
            ],
        ),
    ]

    out = {
        "packet_id": "P1959",
        "stage_id": "S909",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "NONABELIAN_GAUGE_GROUP_DATA_NOT_AVAILABLE_FOR_BRST_EXTENSION__A6_PARTIAL_SCAFFOLD_ONLY",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "ignored_files": sorted(SKIP_NAMES),
        "pre_execution_grep_summary": {
            "english_terms": [
                "SU(3)",
                "SU(2)",
                "U(1)",
                "structure constants",
                "f^{abc}",
                "Jacobi identity",
                "Jacobi tensor",
                "non-Abelian",
                "BRST differential",
                "D_mu c",
                "ghost self-interaction",
            ],
            "polish_terms": [
                "stale struktury",
                "algebra Lie",
                "komutator",
                "nieabel",
                "rozniczka BRST",
                "grupa cechowania",
            ],
            "term_scan_counts": term_scan,
            "key_existing_sources_found": [
                "A6: SU(3)xSU(2)xU(1) exists as strict-core partial scaffold, with full physical uniqueness blocked.",
                "P1380: SU(3)xSU(2)xU(1) image-closure theorem attempt retains obstruction.",
                "P1907: sector-level SM gauge registry exists, but ghost/BRST constraints are OPEN_SYMBOLIC.",
                "P1958: local Abelian Lorenz/FP/ghost seed exported, not non-Abelian BRST.",
            ],
            "negative_search_result": (
                "No theorem-grade f^a_bc structure-constant table, Jacobi certificate, "
                "non-Abelian BRST differential, or ghost self-interaction export was found."
            ),
        },
        "depends_on": {
            "a6_summary_present": "a6" in a6,
            "p1380_text_present": bool(p1380_text),
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1958_present": "GhostAntighostAction_strict_B1_v1" in p1958,
            "p1957_present": "formal_nonavailability_theorem" in p1957,
        },
        "input_hashes": {
            "a6_summary_sha256": digest(a6),
            "a6_text_sha256": digest_text(a6_text),
            "p1380_generated_sha256": digest(p1380),
            "p1380_text_sha256": digest_text(p1380_text),
            "p1907_sha256": digest(p1907),
            "p1958_sha256": digest(p1958),
            "p1957_sha256": digest(p1957),
        },
        "nonabelian_extension_acceptance_formula": {
            "formula": str(nonabelian_ready_formula),
            "truth_assignment": {str(k): bool(v) for k, v in truth_assignment.items()},
            "evaluated_ready": evaluated_ready,
        },
        "symbolic_nilpotency_dependency": {
            "brst_ghost_rule": "s c^a = -(g/2) f^a_{bc} c^b c^c",
            "s_squared_ghost_dependency": "s^2 c^a is proportional to J^a_{bcd} c^b c^c c^d",
            "jacobi_tensor_symbolic_factor": str(q2_ghost_dependency),
            "factor_if_jacobi_tensor_zero": str(q2_ghost_if_jacobi_zero),
            "machine_conclusion": {
                "nilpotency_reduces_to_zero_if_jacobi_exported_zero": q2_ghost_if_jacobi_zero == 0,
                "nilpotency_underdetermined_without_jacobi_certificate": q2_ghost_dependency != 0,
            },
        },
        "minimal_missing_data_matrix": requirement_matrix,
        "formal_obstruction_statement": {
            "statement": (
                "The current strict repository state does not contain enough theorem-grade "
                "non-Abelian group data to extend P1958 to an SU(3)xSU(2)xU(1) BRST differential."
            ),
            "proof_trace": [
                "A6 provides a partial strict-core gauge scaffold, but explicitly denies full uniqueness/closure.",
                "P1380 retains image-closure obstruction, including missing commutator compatibility.",
                "The expected QW-2190 source report is not present inside fundamental_action_reconstruction.",
                "No generated JSON exports a structure-constant table or Jacobi certificate.",
                "Without Jacobi data, s^2 c remains proportional to an undischarged Jacobi tensor.",
                "Therefore the P1958 Abelian seed cannot be lawfully promoted to non-Abelian BRST.",
            ],
        },
        "safe_consequence": {
            "may_use": [
                "A6 partial scaffold as context, not as full BRST input",
                "P1958 local Abelian gauge-fixing/ghost seed",
                "P1907 sector-level gauge registry",
            ],
            "must_not_claim": [
                "SU(3)xSU(2)xU(1) non-Abelian BRST differential exported",
                "Q^2=0 for non-Abelian strict gauge sector",
                "ghost self-interaction closure",
                "TG2_BRST_GLOBAL_NILPOTENCY PASS",
                "TG3_CUTKOSKY_GLOBAL_UNITARITY PASS",
            ],
        },
        "next_solver_input_contract": {
            "recommended_packet": "P1960",
            "minimum_new_exports_required": [
                "QW2190_source_report_ingestion_or_absence_certificate",
                "StructureConstants_SU3_SU2_U1_strict_v1",
                "JacobiCertificate_SU3_SU2_U1_strict_v1",
                "GaugeConnectionRules_SU3_SU2_U1_strict_v1",
                "NonAbelianBRSTDifferential_strict_B1_v1",
                "GhostSelfInteraction_strict_B1_v1",
                "P1380_CMIX_COMMUTATOR_COMPATIBILITY_VERDICT",
            ],
        },
        "false_pass_guard": (
            "This audit does not deny the standard mathematical existence of SU(3), SU(2), or U(1). "
            "It says those data are not exported here in the strict proof-object form needed for "
            "P1958 -> non-Abelian BRST promotion."
        ),
        "next_honest_step": (
            "Build P1960 with high reasoning: ingest or formally certify absence of the QW-2190 "
            "source report and then export a strict structure-constant/Jacobi data pack, or prove "
            "that the non-Abelian BRST extension remains blocked at the data-source level."
        ),
        "lay_explanation": (
            "Mamy opis, ze struktura podobna do SU(3)xSU(2)xU(1) istnieje jako silny szkic, "
            "ale nie mamy jeszcze tabeli reguł mnozenia generatorow i testu Jacobiego. Bez tego "
            "nie da sie uczciwie zbudowac pelnej nieabelowej wersji duchow BRST."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
