#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

T159_SLOT_FREE_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
N489_THEOREM = (
    ROOT
    / "N489_CURRENT_FIRST_STRICT_T162_SLOT_FREE_SIGMA_INT_TO_THETA_PAIR_SOURCE_DISCHARGE_AND_T159_UPGRADE_THEOREM.md"
)

T165_THETA_FIX = GENERATED / "theta_fix_pair1_o2_cut_ord_reference_v1.json"
F454_SHANNON_ASSIGNMENT = (
    GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"
)
F454_SHANNON_ASSIGNMENT_SUMMARY = (
    GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1_summary.json"
)
N496_THEOREM = (
    ROOT
    / "N496_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR3_TO_PAIR5_O2_TO_Z2_UNIQUENESS_THEOREM.md"
)
F450_THETA_PAIR = GENERATED / "theta_pair_canonical_local_diagonal_strict_derived_v1.json"
P450_R1_POPULATION = (
    GENERATED
    / "r1_residual_orientation_datum_target_slot_population_strict_derived_from_canonical_local_diagonal_theta_pair_v1.json"
)
N487_THEOREM = (
    ROOT
    / "N487_CURRENT_FIRST_STRICT_QW2191_CONTINUOUS_O2_FAMILY_DISCHARGE_ON_ALL_FOURIER_DEGENERATE_PAIRS_VIA_CANONICAL_DIAGONAL_LOCAL_SECTOR_THEOREM.md"
)

F453_ASSIGNMENT = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
F453_ASSIGNMENT_SUMMARY = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1_summary.json"
N492_THEOREM = (
    ROOT
    / "N492_CURRENT_FIRST_ACTUAL_STRICT_CANONICAL_LOCAL_DIAGONAL_INTERNAL_ORIENTATION_DATUM_EXPORT_THEOREM.md"
)
P452_SUMMARY = (
    GENERATED / "p452_current_strict_qw2191_residual_z2_sign_flip_gauge_equivalence_audit_probe_summary.json"
)
N493_THEOREM = (
    ROOT / "N493_CURRENT_FIRST_STRICT_QW2191_RESIDUAL_Z2_SIGN_FLIP_GAUGE_EQUIVALENCE_THEOREM.md"
)
N494_THEOREM = (
    ROOT
    / "N494_CURRENT_FIRST_STRICT_QW2190_DIAGONAL_LOCAL_MODE_INDEX_CANONICALIZATION_UNIQUENESS_UP_TO_CONJUGATION_THEOREM.md"
)

P454_SUMMARY = GENERATED / "p454_current_strict_qw2191_o2_rotation_gauge_equivalence_audit_probe_summary.json"
N495_THEOREM = (
    ROOT
    / "N495_CURRENT_FIRST_STRICT_QW2191_O2_ROTATION_GAUGE_EQUIVALENCE_FOR_QW2190_EMBEDDING_AUDITS_THEOREM.md"
)

N490_THEOREM = (
    ROOT
    / "N490_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_DISCHARGE_THEOREM.md"
)

IOTA_SUPPORT = GENERATED / "iota_residual_datum_sigma_int_bridge_export_map_object_support_v1.json"

N491_THEOREM = (
    ROOT / "N491_CURRENT_FIRST_ACTUAL_STRICT_T2_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_DISCHARGE_THEOREM.md"
)

P438_SCRIPT = ROOT / "p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe.py"
P439_SCRIPT = ROOT / "p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe.py"

P438_SUMMARY = GENERATED / "p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe_summary.json"
P439_SUMMARY = GENERATED / "p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe_summary.json"

P455_SUMMARY = GENERATED / "p455_current_strict_mode_index_assignment_shannon_vs_diagonal_alignment_audit_probe_summary.json"

OUT_JSON = GENERATED / "p441_current_strict_global_closure_next_move_dashboard_probe.json"
OUT_SUMMARY = GENERATED / "p441_current_strict_global_closure_next_move_dashboard_probe_summary.json"

T170_GLOBAL_ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
N515_THEOREM = (
    ROOT / "N515_CURRENT_FIRST_STRICT_T170_GLOBAL_SELECTOR_ATLAS_AND_TRANSITION_OBJECT_DISCHARGE_THEOREM.md"
)
H39_GLOBAL_SELECTOR_STATE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
N516_THEOREM = (
    ROOT / "N516_CURRENT_FIRST_STRICT_H39_GLOBAL_PROJECTIVE_SELECTOR_STATE_OBJECT_EXPORT_THEOREM.md"
)
N517_THEOREM = (
    ROOT / "N517_CURRENT_FIRST_STRICT_H37_EVEN_REFERENCE_WEIGHTS_SIGN_DISTINCTION_OBSTRUCTION_THEOREM.md"
)
N518_THEOREM = (
    ROOT
    / "N518_CURRENT_FIRST_STRICT_H37_AUT_INVARIANT_REFERENCE_WEIGHTS_SIGN_DISTINCTION_OBSTRUCTION_THEOREM.md"
)
P472_SUMMARY = (
    GENERATED / "p472_current_strict_pair1_reflection_breaking_weight_repo_scan_probe_summary.json"
)
P473_SUMMARY = (
    GENERATED
    / "p473_current_strict_extension_lane_global_oriented_selector_state_projector_consistency_audit_probe_summary.json"
)
P474_SUMMARY = (
    GENERATED
    / "p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe_summary.json"
)
P475_SUMMARY = (
    GENERATED
    / "p475_current_strict_projective_only_continuation_decision_packet_summary.json"
)
P16_FACTORIZATION_SUMMARY = (
    GENERATED
    / "p16_existing_kernel_feedback_legacy_chart_reduced_operator_export_probe_summary.json"
)
P477_SUMMARY = (
    GENERATED
    / "p477_current_strict_r18_pair1_residual_zero_equations_value_instantiation_probe_summary.json"
)
P478_SUMMARY = (
    GENERATED
    / "p478_current_strict_t169_rordpow_sign_scan_for_r18_pair1_zero_equations_probe_summary.json"
)
N520_SUMMARY = (
    GENERATED
    / "n520_current_first_strict_r18_declared_pair1_residual_zero_equations_value_instance_obstruction_theorem_summary.json"
)
P61_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p61_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_source_eom_coefficient_defect_polynomial_packet_summary.json"
)
P46_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p46_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_action_coefficient_defect_polynomial_packet_summary.json"
)
P622_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p622_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_role_matching_packets_summary.json"
)
P623_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p623_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_common_plus3_assignment_slot_split_packets_summary.json"
)
P624_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p624_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi10_target_role_split_and_defect_polynomial_packets_summary.json"
)
P625_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p625_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi10_target_coherence_instance_summary.json"
)
P626_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p626_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_eom_coherence_instance_summary.json"
)
P627_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p627_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_slotwise_role_split_and_defect_polynomial_packets_summary.json"
)
P628_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p628_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi2_source_psi5_target_psi8_source_psi11_target_coherence_instances_summary.json"
)
P629_DIRECT_FORMAL_SUMMARY = (
    GENERATED
    / "p629_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_t169_constrained_lift_g4_g6_family_shift_defect_zero_witness_packet_summary.json"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "P441 global dashboard: summarize the strict closure next-move state across the two current QW-2191 attack lanes "
            "(diagonal/local accelerator vs Shannon symmetry-breaking audit)."
        )
    )
    p.add_argument(
        "--full-scan",
        action="store_true",
        help="Run the diagonal accelerator dashboard (P438) with --full-scan (slower but repo-wide, includes cache/vendor dirs).",
    )
    return p.parse_args()


def run_script(path: Path, extra_args: list[str] | None = None) -> dict[str, Any]:
    extra_args = extra_args or []
    proc = subprocess.run(
        ["python3", str(path), *extra_args],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "script": str(path.relative_to(REPO)),
        "args": extra_args,
        "returncode": proc.returncode,
        "stdout_tail": "\n".join(proc.stdout.splitlines()[-10:]),
        "stderr_tail": "\n".join(proc.stderr.splitlines()[-10:]),
    }


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _p16_value_instance_obstruction_note() -> str:
    notes: list[str] = []
    if N520_SUMMARY.exists():
        try:
            n520 = load_json(N520_SUMMARY)
            tr = n520.get("theorem_result")
            if isinstance(tr, dict) and bool(tr.get("discharged")):
                violated = tr.get("violated_equations")
                notes.append(
                    " Evidence: N520 packages a value-instance-only obstruction (from P477): under the current exported strict-derived "
                    f"value instance the R18 declared pair1 residual zero equations are violated (violated_equations={violated})."
                )
        except Exception:
            pass
    if P477_SUMMARY.exists():
        try:
            p477 = load_json(P477_SUMMARY)
            if p477.get("status") == "PASS_EVALUATION_ZERO_EQUATIONS_VIOLATED_UNDER_CURRENT_VALUE_INSTANCE":
                violated = p477.get("violated_equations")
                notes.append(
                    " Evidence: P477 evaluates the R18 declared pair1 residual zero equations on the current exported strict-derived "
                    f"value instance and reports violations (violated_equations={violated})."
                )
        except Exception:
            pass
    if P478_SUMMARY.exists():
        try:
            p478 = load_json(P478_SUMMARY)
            if (
                p478.get("status")
                == "PASS_SCAN_COMPLETE_NO_SIGN_VECTOR_SATISFIES_R18_PAIR1_ZERO_EQUATIONS_UNDER_N477_ON_RORDPOW_MAGNITUDES"
            ):
                min_abs = p478.get("min_abs_by_entry_over_scan")
                notes.append(
                    " Evidence: P478 exhaustively scans the full finite T169 r_ordpow sign space (fixing magnitudes and g4 as in F447) "
                    "under conditional N477 and reports no sign vector satisfies all three R18 declared pair1 zero equations; "
                    f"min_abs_by_entry_over_scan={min_abs}."
                )
        except Exception:
            pass
    return "".join(notes)


def _projective_only_continuation_selected() -> tuple[bool, dict[str, Any] | None]:
    if not P475_SUMMARY.exists():
        return False, None
    try:
        p475 = load_json(P475_SUMMARY)
    except Exception:
        return False, None
    selected = (p475.get("decision") == "PROJECTIVE_ONLY_CONTINUATION_SELECTED") or (
        p475.get("selected_continuation") == "projective_only"
    )
    return bool(selected), p475


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    args = parse_args()

    runs = [
        run_script(P438_SCRIPT, ["--full-scan"] if args.full_scan else []),
        run_script(P439_SCRIPT),
    ]

    missing_files: list[str] = []
    for p in (P438_SUMMARY, P439_SUMMARY):
        if not p.exists():
            missing_files.append(str(p.relative_to(REPO)))

    if missing_files:
        summary = {
            "stage": "P441",
            "status": "NOT_COMPUTABLE_MISSING_DEPENDENCY_SUMMARIES",
            "missing_summary_files": missing_files,
            "dependency_runs": runs,
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p438 = load_json(P438_SUMMARY)
    p439 = load_json(P439_SUMMARY)

    # Conservative global recommendation:
    # - the diagonal/local lane (T166) is the only currently exported strict candidate accelerator against QW-2191 on pair1;
    # - P438 already computes which strict target blocks it (T167 vs T168/T169);
    # - the Shannon lane exports strict selector ingredients (T165, and a full n=12 all-pairs mode-index assignment via F454),
    #   but it still does not decide the diagonal/local defect (T166) nor does it globally discharge QW-2191.
    recommended_next = str(p438.get("recommended_next_strict_target") or "T168")
    recommendation_reason = str(
        p438.get("recommendation_reason")
        or "P438 diagonal accelerator lane determines the next strict missing target beneath T166."
    )

    diagonal_note = "P438 summarizes the diagonal/local lane computability and next strict target."
    if N487_THEOREM.exists():
        diagonal_note = (
            "Diagonal/local lane packages a scoped discharge of the continuous QW-2191 O(2) family on all Fourier-degenerate "
            "pairs on n=12 via the canonical diagonal/local sector (N487). This is not ToE closure."
        )
    if F450_THETA_PAIR.exists():
        diagonal_note += (
            " The same lane also exports a strict-derived diagonal/local theta-pair source (F450) and a strict-derived "
            "R1 target-slot inhabitant instance constructed from it (P450), without promoting sigma-int corridor upgrade (T159)."
            if P450_R1_POPULATION.exists()
            else " The same lane also exports a strict-derived diagonal/local theta-pair source (F450)."
        )

    if F453_ASSIGNMENT.exists() or F453_ASSIGNMENT_SUMMARY.exists() or N492_THEOREM.exists():
        diagonal_note += (
            " Additionally, the lane exports an explicit strict-derived mode-index assignment basis object for all Fourier-degenerate "
            "pair planes (F453) and packages it as an internal orientation datum in the lane-scoped sense (N492), without implying global selector closure."
        )

    if P452_SUMMARY.exists() or N493_THEOREM.exists() or N494_THEOREM.exists():
        diagonal_note += (
            " Residual Z2 sign flips are audited to be conjugation-only for the QW-2190 SU(3)/SU(2) embedding audits (P452/N493); "
            "therefore the diagonal/local lane canonicalizes the embedding uniquely up to conjugation in its declared scope (N494)."
        )

    if P454_SUMMARY.exists() or N495_THEOREM.exists():
        diagonal_note += (
            " Moreover, full O(2) basis rotations are gauge-equivalent (conjugation-only) for the same QW-2190 embedding audits "
            "(N495; audited by P454), so the continuous O(2) family is not a changing observable at the embedding-audit level."
        )

    # If the sigma-int corridor strict-core theta supply is exported (T159 satisfied via T162),
    # the next honest frontier is no longer "T159", but the post-T148 object-support layer (T130/N302).
    if recommended_next == "T159" and (T159_SLOT_FREE_THETA_PAIR.exists() or N489_THEOREM.exists()):
        recommended_next = "T130"
        recommendation_reason = (
            "T159 strict-core theta supply is exported via a slot-free construction class (F451/N489). "
            "The next honest strict frontier is the post-T148 object-support target above the export-map object (T130), "
            "with N302 kept explicit."
        )
        diagonal_note += " Sigma-int corridor now exports slot-free theta supply (T159 via T162: F451/N489); next frontier shifts to T130."

    if recommended_next == "T130" and (IOTA_SUPPORT.exists() or N490_THEOREM.exists()):
        recommended_next = "T2"
        recommendation_reason = (
            "Post-T148 object-support above the exported sigma-int -> residual export-map object is now exported (F452/N490). "
            "The next honest strict frontier is theorem-level discharge of the conditional bridge theorem (T2) and/or "
            "continuation under explicit QW-2191 discipline (no implied selector closure)."
        )
        diagonal_note += " Sigma-int post-map object support is now exported (F452/N490); next frontier shifts to T2."

    if recommended_next == "T2" and N491_THEOREM.exists():
        recommended_next = "T170"
        recommendation_reason = (
            "T2 theorem-level bridge discharge is now exported (N491). The remaining strict frontier is explicit continuation "
            "under QW-2191 discipline (no implied selector closure): discharge the strict global selector atlas + transition/gluing "
            "object target (T170), keeping all residual sign handling explicit."
        )
        diagonal_note += " T2 theorem-level bridge discharge is now exported (N491); next frontier shifts to T170."

    if recommended_next == "T170" and (T170_GLOBAL_ATLAS.exists() or N515_THEOREM.exists()):
        recommended_next = "H39"
        recommendation_reason = (
            "T170 is now discharged at object level (F469/N515 export a global selector atlas + transition object on C_v1). "
            "The next honest strict frontier is the absence of a global physical selector object beyond chart locality (H39), "
            "and continuation under explicit QW-2191 discipline (no implied selector closure)."
        )
        diagonal_note += " T170 is now discharged by F469/N515; next frontier shifts to H39."

    if recommended_next == "H39" and (H39_GLOBAL_SELECTOR_STATE.exists() or N516_THEOREM.exists()):
        recommended_next = "H37"
        recommendation_reason = (
            "H39 object-existence layer is now resolved: a global projective selector state object is exported on C_v1 (F470/N516). "
            "The next honest strict frontier is to export a sign-sensitive/directed selector state datum or observable (H37/H36) "
            "and continue under explicit QW-2191 discipline (no implied selector closure)."
        )
        diagonal_note += " H39 object-existence layer is now resolved by F470/N516; next frontier shifts to H37."

    projective_selected, p475 = _projective_only_continuation_selected()

    if recommended_next == "H37":
        if N518_THEOREM.exists():
            note = (
                " Note: N518 strengthens N517: any direction-free Aut(Z_12)-invariant reference weight family cannot distinguish sign on the current exported pair1 sine axis "
                "via a scalar of the form Σ_x w(x) u_1(x); therefore H37 requires an explicit reflection-breaking/orientation source or a different observable class."
            )
            if "N518 strengthens N517" not in recommendation_reason:
                recommendation_reason += note
        elif N517_THEOREM.exists():
            note = (
                " Note: N517 shows even ord-reference weights (ord_Z12 / r_ord) cannot distinguish sign on the current exported pair1 sine axis "
                "via a scalar of the form Σ_x w(x) u_1(x); therefore H37 requires an explicit reflection-breaking/orientation source or a different observable class."
            )
            if "N517 shows even ord-reference weights" not in recommendation_reason:
                recommendation_reason += note
        if P472_SUMMARY.exists() and "Probe note: P472" not in recommendation_reason:
            try:
                p472 = load_json(P472_SUMMARY)
                outside = int(
                    p472.get(
                        "candidates_weight_like_strictish_reflection_breaking_and_dot_nonzero_outside_k_total_rows"
                    )
                    or 0
                )
                if outside == 0:
                    recommendation_reason += (
                        " Probe note: P472 scans exported generated artifacts and finds no strict(-derived) weight-like reflection-breaking per-site arrays "
                        "outside non-canonical marked-site K_total row vectors; therefore H37 still requires an explicit reflection-breaking/orientation source."
                    )
            except Exception:
                pass
        if P473_SUMMARY.exists() and "Extension-lane note: P473" not in recommendation_reason:
            try:
                p473 = load_json(P473_SUMMARY)
                if bool(p473.get("overall_pass")):
                    recommendation_reason += (
                        " Extension-lane note: P473 audits that the non-strict AX29 global oriented vector representative "
                        "is projector-consistent with the strict global projective selector state (F470), so the extension-lane sign-fix "
                        "changes only a sign-gauge representative and does not change strict core."
                    )
            except Exception:
                pass
        if P474_SUMMARY.exists() and "Strict note: P474" not in recommendation_reason:
            try:
                p474 = load_json(P474_SUMMARY)
                if bool(p474.get("overall_pass")):
                    recommendation_reason += (
                        " Strict note: P474 audits that the exported global projective selector state object is projector-level "
                        "glued/transported consistently by the exported global selector transition operators on {pair1..pair5}; "
                        "this is ray/projector-level only and does not lift residual sign."
                    )
            except Exception:
                pass

        # Professorial continuation decision:
        # - If projective-only continuation is explicitly selected (P475), we do not pursue a directed sign-lift in strict core (H37/T171).
        # - We keep H37 as an open frontier for the directed branch, but the recommended next move shifts to strict-only ToE closure tasks
        #   that do not require a sign-sensitive orientation datum (projective-only compatible).
        if projective_selected:
            if P16_FACTORIZATION_SUMMARY.exists():
                recommended_next = "P16"
                p16_note = ""
                try:
                    p16 = load_json(P16_FACTORIZATION_SUMMARY)
                    missing = p16.get("remaining_missing_upstream_objects")
                    req = p16.get("required_next_step")
                    p16_note = f" Current factorization frontier: P16. remaining_missing={missing}; required_next_step={req}."
                    p16_note += _p16_value_instance_obstruction_note()
                except Exception:
                    p16_note = " P16 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict ToE-closure bottleneck: continue the existing-kernel-feedback -> explicit H3 factorization route; "
                    "the current frontier is P16 (legacy chart-reduced operator export route)."
                    + p16_note
                )
            elif P629_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P629"
                p_note = ""
                try:
                    p629 = load_json(P629_DIRECT_FORMAL_SUMMARY)
                    missing = p629.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P629. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P629 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P629."
                    + p_note
                )
            elif P628_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P628"
                p_note = ""
                try:
                    p628 = load_json(P628_DIRECT_FORMAL_SUMMARY)
                    missing = p628.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P628. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P628 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P628."
                    + p_note
                )
            elif P627_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P627"
                p_note = ""
                try:
                    p627 = load_json(P627_DIRECT_FORMAL_SUMMARY)
                    missing = p627.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P627. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P627 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P627."
                    + p_note
                )
            elif P626_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P626"
                p_note = ""
                try:
                    p626 = load_json(P626_DIRECT_FORMAL_SUMMARY)
                    missing = p626.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P626. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P626 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P626."
                    + p_note
                )
            elif P625_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P625"
                p_note = ""
                try:
                    p625 = load_json(P625_DIRECT_FORMAL_SUMMARY)
                    missing = p625.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P625. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P625 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P625."
                    + p_note
                )
            elif P624_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P624"
                p_note = ""
                try:
                    p624 = load_json(P624_DIRECT_FORMAL_SUMMARY)
                    missing = p624.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P624. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P624 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P624."
                    + p_note
                )
            elif P623_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P623"
                p_note = ""
                try:
                    p623 = load_json(P623_DIRECT_FORMAL_SUMMARY)
                    missing = p623.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P623. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P623 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P623."
                    + p_note
                )
            elif P622_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P622"
                p_note = ""
                try:
                    p622 = load_json(P622_DIRECT_FORMAL_SUMMARY)
                    missing = p622.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P622. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P622 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P622."
                    + p_note
                )
            elif P61_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P61"
                p_note = ""
                try:
                    p61 = load_json(P61_DIRECT_FORMAL_SUMMARY)
                    missing = p61.get("remaining_missing_upstream_objects")
                    p_note = f" Current direct-formal frontier: P61. remaining_missing_upstream_objects={missing}."
                except Exception:
                    p_note = " P61 summary exists but could not be parsed."

                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P61."
                    + p_note
                )
            elif P46_DIRECT_FORMAL_SUMMARY.exists():
                recommended_next = "P46"
                p_note = ""
                try:
                    p46 = load_json(P46_DIRECT_FORMAL_SUMMARY)
                    missing = p46.get("remaining_missing_objects")
                    p_note = f" Current direct-formal frontier: P46. remaining_missing_objects={missing}."
                except Exception:
                    p_note = " P46 summary exists but could not be parsed."

                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict-only ToE-closure bottleneck: continue on the kernel-split-robust canonical-ontology-supported direct formal c1s1 family route "
                    "(F3 priority), currently tracked at P46."
                    + p_note
                )
            else:
                recommended_next = "P16" if P16_FACTORIZATION_SUMMARY.exists() else "P11"
                p16_note = ""
                if P16_FACTORIZATION_SUMMARY.exists():
                    try:
                        p16 = load_json(P16_FACTORIZATION_SUMMARY)
                        missing = p16.get("remaining_missing_upstream_objects")
                        req = p16.get("required_next_step")
                        p16_note = (
                            f" Current factorization frontier: P16. remaining_missing={missing}; required_next_step={req}."
                        )
                        p16_note += _p16_value_instance_obstruction_note()
                    except Exception:
                        p16_note = " P16 summary exists but could not be parsed."
                recommendation_reason = (
                    "Projective-only continuation is explicitly selected (P475): treat the exported global projective selector state as the strict physical state object "
                    "for the declared closure stack, keeping residual sign as a gauge/convention layer where proven irrelevant (N502, N519). "
                    "H37/T171 remain open for a future directed branch only. "
                    "Next strict ToE-closure bottleneck: continue the existing-kernel-feedback -> explicit H3 factorization route; "
                    "the current frontier is beyond P11 and is tracked by P16 (legacy chart-reduced operator export route)."
                    + p16_note
                )

    # Backward-compatible mapping: older P438 versions returned a packet label ("B3") here.
    if recommended_next == "B3" and N491_THEOREM.exists():
        recommended_next = "T170"
        recommendation_reason = (
            "T2 theorem-level bridge discharge is now exported (N491). The remaining strict frontier is explicit continuation "
            "under QW-2191 discipline (no implied selector closure): discharge the strict global selector atlas + transition/gluing "
            "object target (T170), proceeding via the B3 continuation while keeping residual sign handling explicit."
        )

    shannon_note = "P439 is probe-level only; Shannon-lane status is determined by strict theorem/packet exports."
    if F454_SHANNON_ASSIGNMENT.exists() or F454_SHANNON_ASSIGNMENT_SUMMARY.exists() or N496_THEOREM.exists():
        shannon_note = (
            "Strict Shannon element-order reference lane exports a full strict-core mode-index assignment basis object on n=12 (F454), "
            "cutting O(2) down to residual Z2 on all pair_m (m=1..5) via the cross-entropy objective "
            "(N480, N488, N496)."
        )
        if P455_SUMMARY.exists():
            shannon_note += " A cross-lane hygiene audit reports this axis choice aligns with the diagonal/local assignment up to residual sign (P455)."
        shannon_note += " This remains below T166, does not imply strict-core selector closure, and does not constitute a global QW-2191 discharge."
    elif T165_THETA_FIX.exists():
        shannon_note = (
            "T165 strict selector ingredient is exported (F446/N480): pair1 O(2)->Z2 cut with θ*=π/2 (mod π). "
            "This does not by itself discharge T166 or globally discharge QW-2191."
        )
    elif bool(p439.get("any_objective_unique_cluster_on_grid")):
        shannon_note = (
            "P439 found a unique near-minimum cluster on a theta-grid for at least one audited objective/reference, "
            "but without a theorem-level uniqueness proof this remains probe-only (no strict promotion)."
        )
    elif bool(p439.get("any_objective_z2_unique_on_grid")):
        shannon_note = (
            "P439 found a Z2-unique minimizer pattern on a dense theta-grid for at least one audited non-translation-invariant reference distribution r(x), "
            "but without a theorem-level uniqueness proof this remains probe-only (no strict promotion)."
        )

    artifact = {
        "stage": "P441",
        "goal": "global_dashboard_for_next_strict_move_on_QW2191_closure_attempts_across_diagonal_and_shannon_lanes",
        "inputs": {"p438_full_scan": bool(args.full_scan)},
        "dependency_runs": runs,
        "dependency_summaries": {
            "P438": str(P438_SUMMARY.relative_to(REPO)),
            "P439": str(P439_SUMMARY.relative_to(REPO)),
            "P455": (str(P455_SUMMARY.relative_to(REPO)) if P455_SUMMARY.exists() else None),
        },
        "status_snapshot": {"P438": p438, "P439": p439},
        "result": {
            "recommended_next_strict_target": recommended_next,
            "recommendation_reason": recommendation_reason,
            "diagonal_lane_note": diagonal_note,
            "shannon_lane_note": shannon_note,
            "projective_only_continuation": {
                "selected": bool(projective_selected),
                "p475_summary": (str(P475_SUMMARY.relative_to(REPO)) if p475 is not None else None),
                "note": (
                    "If selected, H37/T171 remain open for directed continuation only; the recommended next move shifts to ToE closure tasks "
                    "that depend only on projectors/spans."
                ),
            },
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P441",
        "status": "PASS_DASHBOARD_READY",
        "recommended_next_strict_target": recommended_next,
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
