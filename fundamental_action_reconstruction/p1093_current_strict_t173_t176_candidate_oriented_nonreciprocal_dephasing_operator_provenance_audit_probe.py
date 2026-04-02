#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-29"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_V1 = ROOT / "V1_INFORMATIONAL_VISCOSITY_HYPOTHESIS_AUDIT.md"
IN_N457 = ROOT / "N457_CURRENT_FIRST_STRICT_PHASE_12_AUT_Z12_QUOTIENT_ORBIT_SPACE_TOPOLOGICAL_HOLONOMY_TRIVIALITY_BOUNDARY_THEOREM.md"
IN_NEURAL = REPO / "raport_qw540_544_neural.md"
IN_REPORT99 = REPO / "report_99_quick_win.md"
IN_REPORT100 = REPO / "report_100_quick_win.md"
IN_TEMP = REPO / "TEMP_FULL_REPORT_CONTEXT.md"
IN_P1092 = GENERATED / "p1092_current_strict_t173_t176_time_arrow_light_resonance_oriented_dephasing_competing_extension_audit_probe_summary.json"
IN_F967 = GENERATED / "f967_current_strict_t173_t176_time_arrow_light_resonance_oriented_dephasing_competing_extension_packet_summary.json"

OUT_JSON = GENERATED / "p1093_current_strict_t173_t176_candidate_oriented_nonreciprocal_dephasing_operator_provenance_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1093_current_strict_t173_t176_candidate_oriented_nonreciprocal_dephasing_operator_provenance_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_V1, IN_N457, IN_NEURAL, IN_REPORT99, IN_REPORT100, IN_TEMP, IN_P1092, IN_F967]
    missing = [rel(p) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1093",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    v1 = IN_V1.read_text(encoding="utf-8")
    n457 = IN_N457.read_text(encoding="utf-8")
    neural = IN_NEURAL.read_text(encoding="utf-8")
    report99 = IN_REPORT99.read_text(encoding="utf-8")
    report100 = IN_REPORT100.read_text(encoding="utf-8")
    temp = IN_TEMP.read_text(encoding="utf-8")
    p1092 = load_json(IN_P1092)
    f967 = load_json(IN_F967)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    add_check(
        "repo_roots_memory_retardation_motif",
        (
            "retard_phase" in v1
            and "selector-sector reduction" in v1
            and "O(2)" in v1
            and "competing extension hypothesis" in v1
        ),
        True,
        "V1 roots memory/retardation motifs but keeps them below solved selector status.",
    )
    add_check(
        "repo_roots_light_resonance_and_phase_modulation_motifs",
        (
            "modulacja fazowa" in report99
            and "Rezonans światło–materia" in report99
            and "Rabi frequency" in report99
            and "oscylacji fazy" in report99
        ),
        True,
        "Report 99 roots light resonance and phase modulation motifs.",
    )
    add_check(
        "repo_roots_chirality_motif",
        (
            "Light + Chirality + Mass" in report100
            and "Chirality measure" in report100
            and "helicity" in report99.lower()
        ),
        True,
        "Reports 99 and 100 root a chirality/helicity motif.",
    )
    add_check(
        "repo_roots_arrow_of_time_motif",
        (
            "Strzałka Czasu" in neural
            and "Entropia wzrosła" in neural
        ),
        True,
        "The old neural report roots an arrow-of-time motif tied to entropy growth.",
    )
    add_check(
        "repo_roots_decoherence_and_resonator_stability_motifs",
        (
            "Szybka dekoherencja" in temp
            and "rezonatorów" in temp
        ),
        True,
        "TEMP_FULL_REPORT_CONTEXT roots decoherence and resonator-stability motifs.",
    )
    add_check(
        "strict_phase_holonomy_alone_is_blocked",
        (
            "cannot support a nontrivial Berry holonomy" in n457
            and "additional typed" in n457
            and "new selector slot" in n457
        ),
        True,
        "N457 blocks pure phase/holonomy alone without extra selector structure.",
    )
    add_check(
        "current_lane_is_only_competing_extension",
        (
            p1092.get("hypothesis_grade") == "competing_extension_hypothesis_only"
            and f967.get("hypothesis_grade") == "competing_extension_hypothesis_only"
            and bool(f967.get("counts_as_lawful_supplier")) is False
            and bool(f967.get("counts_as_strict_physical_orientation_datum")) is False
        ),
        True,
        "P1092/F967 already classify the lane only as a competing extension hypothesis.",
    )

    discharged = len(blocking) == 0
    status = (
        "PASS_CURRENT_STRICT_T173_T176_CANDIDATE_ORIENTED_NONRECIPROCAL_DEPHASING_OPERATOR_PROVENANCE_AUDITED"
        if discharged
        else "P1093_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PROVENANCE_STATE"
    )

    repo_rooted_components = [
        {
            "component_id": "memory_retardation_motif",
            "grade": "repo_rooted_motif_only",
            "source": rel(IN_V1),
            "meaning": "damping / retard_phase / selector-lane competing extension support",
        },
        {
            "component_id": "light_resonance_motif",
            "grade": "repo_rooted_motif_only",
            "source": rel(IN_REPORT99),
            "meaning": "light-from-phase / Rabi / light-matter resonance exploratory support",
        },
        {
            "component_id": "phase_modulation_motif",
            "grade": "repo_rooted_motif_only",
            "source": rel(IN_REPORT99),
            "meaning": "phase modulation and feedback motif",
        },
        {
            "component_id": "chirality_motif",
            "grade": "repo_rooted_motif_only",
            "source": rel(IN_REPORT100),
            "meaning": "light+chirality+mass unification motif",
        },
        {
            "component_id": "arrow_of_time_entropy_motif",
            "grade": "repo_rooted_motif_only",
            "source": rel(IN_NEURAL),
            "meaning": "entropy-growth arrow-of-time motif",
        },
        {
            "component_id": "decoherence_resonator_stability_motif",
            "grade": "repo_rooted_motif_only",
            "source": rel(IN_TEMP),
            "meaning": "decoherence / resonator stability motif",
        },
        {
            "component_id": "strict_route_constraint_against_pure_phase_closure",
            "grade": "strict_boundary_only",
            "source": rel(IN_N457),
            "meaning": "pure phase/holonomy alone cannot close the route without extra selector structure",
        },
    ]

    new_import_components = [
        {
            "component_id": "density_matrix_core_rho_t_x",
            "grade": "genuine_new_import",
            "why_new": "no current strict export instantiates the candidate as a rho-evolution operator",
        },
        {
            "component_id": "spatiotemporal_memory_kernel_K_t_tau_x_xprime",
            "grade": "genuine_new_import",
            "why_new": "repo motifs mention kernel and memory, but no exported strict operator of this full form exists",
        },
        {
            "component_id": "orientation_drift_vector_v",
            "grade": "genuine_new_import",
            "why_new": "a privileged drift/orientation carrier is not exported on this lane",
        },
        {
            "component_id": "orientation_axis_n_and_jump_sector_L_n",
            "grade": "genuine_new_import",
            "why_new": "no current strict selector-safe source-side operator sector of this form is exported",
        },
        {
            "component_id": "nonreciprocal_rate_asymmetry_gamma_n_neq_gamma_minus_n",
            "grade": "genuine_new_import",
            "why_new": "this is the actual symmetry-breaking law, not an already exported repo object",
        },
        {
            "component_id": "oriented_commutator_retardation_term",
            "grade": "genuine_new_import",
            "why_new": "no strict current packet exports this history-dependent oriented commutator term",
        },
        {
            "component_id": "cp_cptp_proof_for_memory_generator",
            "grade": "genuine_new_import",
            "why_new": "physical well-posedness is not exported for this candidate lane",
        },
        {
            "component_id": "exact_reduction_to_active_missing_bridge",
            "grade": "genuine_new_import",
            "why_new": "no strict derivation yet ties this candidate directly to ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
        },
    ]

    artifact = {
        "stage": "P1093",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "V1": rel(IN_V1),
            "N457": rel(IN_N457),
            "raport_qw540_544_neural": rel(IN_NEURAL),
            "report_99_quick_win": rel(IN_REPORT99),
            "report_100_quick_win": rel(IN_REPORT100),
            "TEMP_FULL_REPORT_CONTEXT": rel(IN_TEMP),
            "P1092": rel(IN_P1092),
            "F967": rel(IN_F967),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "repo_rooted_components": repo_rooted_components,
        "new_import_components": new_import_components,
        "classification": {
            "complete_candidate_operator_exported": False,
            "repo_rooted_component_count": len(repo_rooted_components),
            "genuine_new_import_component_count": len(new_import_components),
            "minimal_new_import_boundary": "provenance_safe_orientation_anchor_plus_nonreciprocal_rate_asymmetry_plus_oriented_memory_rule",
            "active_missing_bridge": p1092.get("active_missing_bridge"),
            "counts_as_lawful_supplier": False,
            "counts_as_strict_physical_orientation_datum": False,
            "next_honest_move": "if pursued_further_freeze_only_the_irreducible_new_import_boundary_not_the_full_long_ansatz",
        },
        "hard_limits": [
            "No full operator export.",
            "No lawful supplier export.",
            "No strict physical orientation datum export.",
            "No T183 discharge.",
            "No T176 discharge.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "repo_rooted_component_count": artifact["classification"]["repo_rooted_component_count"],
        "genuine_new_import_component_count": artifact["classification"]["genuine_new_import_component_count"],
        "complete_candidate_operator_exported": artifact["classification"]["complete_candidate_operator_exported"],
        "minimal_new_import_boundary": artifact["classification"]["minimal_new_import_boundary"],
        "active_missing_bridge": artifact["classification"]["active_missing_bridge"],
        "counts_as_lawful_supplier": False,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
