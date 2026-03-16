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

F446_REF = GENERATED / "r_ord_z12_v1_reference_distribution.json"
T165_THETA_FIX = GENERATED / "theta_fix_pair1_o2_cut_ord_reference_v1.json"

P434_INPUT = GENERATED / "p434_input_sigma_opposite_pair_sum_values_candidate.json"

P432_SCRIPT = ROOT / "p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe.py"
P434_SCRIPT = ROOT / "p434_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_evaluation_probe.py"
P444_SCRIPT = ROOT / "p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe.py"
P437_SCRIPT = (
    ROOT
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.py"
)
P449_SCRIPT = ROOT / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe.py"

P432_SUMMARY = (
    GENERATED
    / "p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe_summary.json"
)
P434_SUMMARY = (
    GENERATED
    / "p434_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_evaluation_probe_summary.json"
)
P444_SUMMARY = GENERATED / "p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe_summary.json"
P437_SUMMARY = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json"
)
P449_SUMMARY = (
    GENERATED / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe_summary.json"
)

# If present, this theorem file indicates the T166 value-instantiation decision object is already exported.
N482_THEOREM = ROOT / "N482_CURRENT_FIRST_STRICT_T166_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_NONZERO_DECISION_THEOREM.md"
N487_THEOREM = (
    ROOT
    / "N487_CURRENT_FIRST_STRICT_QW2191_CONTINUOUS_O2_FAMILY_DISCHARGE_ON_ALL_FOURIER_DEGENERATE_PAIRS_VIA_CANONICAL_DIAGONAL_LOCAL_SECTOR_THEOREM.md"
)

N489_THEOREM = (
    ROOT
    / "N489_CURRENT_FIRST_STRICT_T162_SLOT_FREE_SIGMA_INT_TO_THETA_PAIR_SOURCE_DISCHARGE_AND_T159_UPGRADE_THEOREM.md"
)

N490_THEOREM = (
    ROOT
    / "N490_CURRENT_FIRST_ACTUAL_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_DISCHARGE_THEOREM.md"
)

IOTA_SUPPORT = GENERATED / "iota_residual_datum_sigma_int_bridge_export_map_object_support_v1.json"

N491_THEOREM = (
    ROOT / "N491_CURRENT_FIRST_ACTUAL_STRICT_T2_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_DISCHARGE_THEOREM.md"
)

OUT_JSON = (
    GENERATED / "p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe.json"
)
OUT_SUMMARY = (
    GENERATED / "p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe_summary.json"
)

T170_GLOBAL_ATLAS = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
N515_THEOREM = (
    ROOT / "N515_CURRENT_FIRST_STRICT_T170_GLOBAL_SELECTOR_ATLAS_AND_TRANSITION_OBJECT_DISCHARGE_THEOREM.md"
)
H39_GLOBAL_SELECTOR_STATE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
N516_THEOREM = (
    ROOT / "N516_CURRENT_FIRST_STRICT_H39_GLOBAL_PROJECTIVE_SELECTOR_STATE_OBJECT_EXPORT_THEOREM.md"
)
H37_DIRECTED_STATE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"
N524_THEOREM = (
    ROOT
    / "N524_CURRENT_FIRST_STRICT_T171_GLOBAL_DIRECTED_SELECTOR_STATE_DATUM_DISCHARGE_THEOREM.md"
)
P632_SUMMARY = (
    GENERATED / "p632_current_strict_directed_continuation_decision_packet_summary.json"
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
P633_SCRIPT = ROOT / "p633_current_strict_source_seed_route_selection_decision_packet.py"
P633_SUMMARY = (
    GENERATED / "p633_current_strict_source_seed_route_selection_decision_packet_summary.json"
)
P119_SUMMARY = (
    GENERATED / "p119_first_source_seed_construction_target_probe_summary.json"
)
P392_SUMMARY = (
    GENERATED
    / "p392_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_probe_summary.json"
)
P634_SUMMARY = (
    GENERATED
    / "p634_genuinely_new_strict_core_source_object_clause_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
N525_SUMMARY = (
    GENERATED / "n525_current_first_clause_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
F635_SUMMARY = (
    GENERATED
    / "f635_future_genuinely_new_source_object_lift_bind_target_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P635_SUMMARY = (
    GENERATED
    / "p635_future_genuinely_new_source_object_lift_bind_target_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
F636_SUMMARY = (
    GENERATED
    / "f636_first_future_genuinely_new_source_object_lift_bind_attempt_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P636_SUMMARY = (
    GENERATED
    / "p636_first_future_genuinely_new_source_object_lift_bind_attempt_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
F637_SUMMARY = (
    GENERATED
    / "f637_first_future_constructed_source_object_realization_target_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P637_SUMMARY = (
    GENERATED
    / "p637_first_future_constructed_source_object_realization_target_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
N528_SUMMARY = (
    GENERATED
    / "n528_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_target_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
F638_SUMMARY = (
    GENERATED
    / "f638_first_future_constructed_source_object_realization_attempt_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P638_SUMMARY = (
    GENERATED
    / "p638_first_future_constructed_source_object_realization_attempt_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
N529_SUMMARY = (
    GENERATED
    / "n529_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_attempt_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
F639_SUMMARY = (
    GENERATED
    / "f639_first_future_constructed_source_object_realization_verdict_target_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P639_SUMMARY = (
    GENERATED
    / "p639_first_future_constructed_source_object_realization_verdict_target_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
N530_SUMMARY = (
    GENERATED
    / "n530_next_constructive_move_reduced_to_one_first_future_constructed_source_object_realization_verdict_target_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
F640_SUMMARY = (
    GENERATED
    / "f640_first_future_constructed_source_object_realization_verdict_branch_refinement_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P640_SUMMARY = (
    GENERATED
    / "p640_first_future_constructed_source_object_realization_verdict_branch_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
N531_SUMMARY = (
    GENERATED
    / "n531_next_constructive_move_reduced_to_one_explicit_success_failure_branch_split_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
F641_SUMMARY = (
    GENERATED
    / "f641_first_conservative_realization_failure_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P641_SUMMARY = (
    GENERATED
    / "p641_current_failure_verdict_discharge_probe_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
)
N532_SUMMARY = (
    GENERATED
    / "n532_current_failure_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
)
F642_SUMMARY = (
    GENERATED
    / "f642_remaining_realization_success_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P642_SUMMARY = (
    GENERATED
    / "p642_current_success_verdict_discharge_probe_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
)
N533_SUMMARY = (
    GENERATED
    / "n533_current_success_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_realization_attempt_summary.json"
)
F643_SUMMARY = (
    GENERATED
    / "f643_first_post_verdict_admissibility_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P643_SUMMARY = (
    GENERATED
    / "p643_current_admissibility_branch_discharge_probe_for_future_constructed_source_object_for_s_sel_int_candidate_seed_v1_summary.json"
)
N534_SUMMARY = (
    GENERATED
    / "n534_current_admissibility_branch_obstruction_theorem_for_future_constructed_source_object_for_s_sel_int_candidate_seed_v1_summary.json"
)
F644_SUMMARY = (
    GENERATED
    / "f644_first_post_admissibility_orientation_export_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P644_SUMMARY = (
    GENERATED
    / "p644_current_orientation_export_branch_discharge_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
N535_SUMMARY = (
    GENERATED
    / "n535_current_orientation_export_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
F645_SUMMARY = (
    GENERATED
    / "f645_last_downstream_completion_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)
P645_SUMMARY = (
    GENERATED
    / "p645_current_downstream_completion_branch_discharge_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)
N536_SUMMARY = (
    GENERATED
    / "n536_current_downstream_completion_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
N537_SUMMARY = (
    GENERATED
    / "n537_current_post_verdict_lower_branch_full_negative_closure_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
)
F646_SUMMARY = (
    GENERATED
    / "f646_strict_witness_provider_contract_packet_for_seed_v1_realization_attempt_summary.json"
)
P646_SUMMARY = (
    GENERATED
    / "p646_current_strict_witness_provider_scan_probe_for_seed_v1_realization_attempt_summary.json"
)
N538_SUMMARY = (
    GENERATED
    / "n538_current_strict_witness_provider_absence_theorem_for_seed_v1_realization_attempt_summary.json"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "P438 dashboard: run P432/P444/P437/P434 and summarize whether T166 is computable, plus the next strict target."
        )
    )
    p.add_argument(
        "--full-scan",
        action="store_true",
        help="Run P432 and P444 with --full-scan (include vendor/cache dirs). Defaults to the faster hygiene scan.",
    )
    return p.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def has_explicit_strict_derived_marker(obj: dict[str, Any]) -> bool:
    status = str(obj.get("status") or "")
    classification = str(obj.get("classification") or "")
    return ("STRICT_DERIVED" in status.upper()) or ("strict_derived" in status.lower()) or ("strict_derived" in classification.lower())


def sigma_six_keys() -> list[str]:
    return [f"Sigma_psi{i}_psi{i+6}" for i in range(6)]


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


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    args = parse_args()

    runs = [
        run_script(P432_SCRIPT, ["--full-scan"] if args.full_scan else []),
        run_script(P444_SCRIPT, ["--full-scan"] if args.full_scan else []),
        run_script(P437_SCRIPT),
        run_script(P449_SCRIPT),
        run_script(P434_SCRIPT),
        run_script(P633_SCRIPT),
    ]

    missing_files: list[str] = []
    for p in (P432_SUMMARY, P444_SUMMARY, P437_SUMMARY, P449_SUMMARY, P434_SUMMARY, P633_SUMMARY):
        if not p.exists():
            missing_files.append(str(p.relative_to(REPO)))

    if missing_files:
        summary = {
            "stage": "P438",
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

    p432 = load_json(P432_SUMMARY)
    p444 = load_json(P444_SUMMARY)
    p437 = load_json(P437_SUMMARY)
    p449 = load_json(P449_SUMMARY)
    p434 = load_json(P434_SUMMARY)

    decision_ready = bool(p432.get("decision_ready_from_repo_values"))
    t165_selector_ingredient_exported = bool(F446_REF.exists() and T165_THETA_FIX.exists())
    t166_decision_theorem_exported = bool(N482_THEOREM.exists())
    q2191_scoped_discharge_theorem_exported = bool(N487_THEOREM.exists())

    t167_strict_sigma_input_present = False
    if P434_INPUT.exists():
        try:
            p434_in = load_json(P434_INPUT)
            sigma_keys = sigma_six_keys()
            sigma_numeric = all(isinstance(p434_in.get(k), (int, float)) for k in sigma_keys)
            sigma_marked_strict = has_explicit_strict_derived_marker(p434_in)
            provider_marked_strict = False
            src = p434_in.get("source_provider")
            if isinstance(src, str):
                src_path = Path(src)
                if not src_path.is_absolute():
                    src_path = REPO / src_path
                if src_path.exists():
                    provider_marked_strict = has_explicit_strict_derived_marker(load_json(src_path))
            t167_strict_sigma_input_present = bool(sigma_numeric and sigma_marked_strict and provider_marked_strict)
        except Exception:
            t167_strict_sigma_input_present = False

    recommended_next_target = "T167"
    recommendation_reason = "P434 sigma-six-sum input is not exported as an explicitly strict-derived value instantiation today; discharge T167."
    if p437.get("status") == "NOT_COMPUTABLE_MISSING_INPUT_VALUES":
        recommended_next_target = "T168"
        recommendation_reason = "P437 is blocked on missing upstream values (vpsi,g4,g6); discharge strict vacuum/self-coupling value-provider target"
        if t165_selector_ingredient_exported:
            recommended_next_target = "T169"
            recommendation_reason = (
                "F446/N480 export a direction-free strict selector ingredient (T165) on Z_12; "
                "the remaining strict bottleneck under the T166 diagonal lane is now best stated as the constrained lift "
                "from strict scalar closure + strict selector ingredient into a T168-consumable per-site provider (T169)."
            )
    elif p434.get("status") == "PASS_EVALUATED" and t167_strict_sigma_input_present:
        recommended_next_target = "T166"
        recommendation_reason = (
            "P434 evaluated F2(d) from an explicit strict-derived sigma-six-sum instantiation; discharge the T166 decision "
            "object/theorem (F2(d)=0 vs F2(d)≠0) and export the induced pair1 axis datum."
        )
        if t166_decision_theorem_exported:
            recommended_next_target = "T159"
            recommendation_reason = (
                "T166 value-instantiation decision theorem is already exported (N482). The next strict frontier is now beyond "
                "the diagonal/local lane: proceed to the strict sigma-int → theta strict-core upgrade target (T159)."
            )
            if q2191_scoped_discharge_theorem_exported and bool(p449.get("all_pairs_cut")):
                recommendation_reason = (
                    "Diagonal/local lane now packages a scoped discharge of the continuous QW-2191 O(2) family on all "
                    "Fourier-degenerate pairs on n=12 (N487, supported by P449/N485). The next strict frontier is the "
                    "sigma-int → theta strict-core upgrade target (T159) and/or upgrading the diagonal decision from "
                    "value-instantiation to coefficient-class discharge (N472/P431)."
                )

    if recommended_next_target == "T159" and N489_THEOREM.exists():
        recommended_next_target = "T130"
        recommendation_reason = (
            "T159 strict-core theta supply is now exported via a slot-free construction class (N489). "
            "The next honest strict frontier is the post-T148 object-support target above the export-map object (T130), "
            "with N302 kept explicit."
        )

    if recommended_next_target == "T130" and (IOTA_SUPPORT.exists() or N490_THEOREM.exists()):
        recommended_next_target = "T2"
        recommendation_reason = (
            "Post-T148 object-support above the exported sigma-int -> residual export-map object is now exported (F452/N490). "
            "The next honest strict frontier is theorem-level discharge of the conditional bridge theorem (T2) and/or "
            "continuation under explicit QW-2191 discipline (no implied selector closure)."
        )

    if recommended_next_target == "T2" and N491_THEOREM.exists():
        recommended_next_target = "T170"
        recommendation_reason = (
            "T2 theorem-level bridge discharge is now exported (N491). The remaining strict frontier is explicit continuation "
            "under QW-2191 discipline (no implied selector closure): discharge the strict global selector atlas + transition/gluing "
            "object target (T170), proceeding via the B3 topological selector bridge continuation while keeping the residual "
            "Z2 sign handling explicit (e.g. sign gauge-irrelevance where proven, or sign frozen as a tracked convention layer)."
        )

    if recommended_next_target == "T170" and (T170_GLOBAL_ATLAS.exists() or N515_THEOREM.exists()):
        recommended_next_target = "H39"
        recommendation_reason = (
            "T170 is now discharged at object level (F469/N515 export a global selector atlas + transition object on C_v1). "
            "The next honest strict frontier is the absence of a global physical selector object beyond chart locality (H39), "
            "and continuation under explicit QW-2191 discipline (no implied selector closure)."
        )

    if recommended_next_target == "H39" and (H39_GLOBAL_SELECTOR_STATE.exists() or N516_THEOREM.exists()):
        recommended_next_target = "H37"
        recommendation_reason = (
            "H39 object-existence layer is now resolved: a global projective selector state object is exported on C_v1 (F470/N516). "
            "The next honest strict frontier is to export a sign-sensitive/directed selector state datum or observable (H37/H36) "
            "and continue under explicit QW-2191 discipline (no implied selector closure)."
        )

    # Post-projective directed frontier resolved: if the directed selector state datum is now exported (T171 discharged),
    # do not keep recommending H37.
    if recommended_next_target == "H37" and (H37_DIRECTED_STATE.exists() or N524_THEOREM.exists()):
        # If both main strict-only ToE closure continuations are already frozen negative (P480 + P631),
        # shift the next move to the strict-core source-seed construction frontier (P119), but only if
        # the explicit professorial decision packet is present and selected (P633).
        if P633_SUMMARY.exists():
            try:
                p633 = load_json(P633_SUMMARY)
                if str(p633.get("decision") or "") == "STRICT_CORE_SOURCE_SEED_ROUTE_SELECTED":
                    recommended_next_target = "P119"
                    recommendation_reason = (
                        "Post-projective directed frontier is now resolved (T171 discharged), and both previously pursued strict-only ToE closure continuations "
                        "are explicitly frozen negative on the current strict branch (P480: freeze P16; P631: freeze direct-formal residual-cancellation on T166 nonzero). "
                        "Therefore the next honest strict bottleneck shifts to the genuinely-new strict-core internal selector source seed construction frontier "
                        "for S_sel_int, tracked by the first source-seed construction target probe (P119). Professorial routing decision: P633."
                    )
                else:
                    recommended_next_target = "P11"
            except Exception:
                recommended_next_target = "P11"
        else:
            recommended_next_target = "P11"
        p_note = ""
        if P632_SUMMARY.exists():
            try:
                p632 = load_json(P632_SUMMARY)
                if str(p632.get("decision") or "") == "DIRECTED_CONTINUATION_SELECTED":
                    p_note = " Professorial note: directed continuation is explicitly selected (P632)."
            except Exception:
                p_note = " Professorial note: P632 summary exists but could not be parsed."
        if recommended_next_target == "P11":
            recommendation_reason = (
                "Post-projective directed frontier is now resolved: a premise-based strict fixing datum (T164) and a strict global directed selector state datum (T171) are exported "
                "(F473/N523 and F474/N524). Therefore H37 is discharged in the declared directed scope and is no longer the next strict blocker."
                " Next strict bottleneck shifts back to strict-only ToE closure tasks that can now treat the selector state as directed where needed (P11)."
                + p_note
            )

    # Source-seed route continuation: advance beyond the initial target probe (P119) when downstream strict artifacts exist.
    if recommended_next_target == "P119" and P119_SUMMARY.exists():
        try:
            p119 = load_json(P119_SUMMARY)
            reduced = bool(
                p119.get("target_state", {}).get(
                    "last_positive_branch_reduced_to_one_first_source_seed_target"
                )
            )
            if reduced:
                if not P392_SUMMARY.exists():
                    recommended_next_target = "P392"
                    recommendation_reason = (
                        "P119 fixes S_sel_int as the next strict-core source-seed construction target. "
                        "The next required strict artifact is the strict sigma-int upgraded seed candidate instance S_sel_int_candidate_seed_v1 (F318), "
                        "probed by P392."
                    )
                elif not P634_SUMMARY.exists():
                    recommended_next_target = "P634"
                    recommendation_reason = (
                        "P392 confirms the strict sigma-int upgraded seed candidate instance S_sel_int_candidate_seed_v1 is exported. "
                        "The next honest move is to rerun the first admissibility clause "
                        "(genuinely_new_strict_core_source_object_required) on seed v1 (P634), "
                        "without implying admissible S_sel_int nor strict-core selector closure."
                    )
                else:
                    p634 = load_json(P634_SUMMARY)
                    first_ok = bool(p634.get("clause_test_result", {}).get("currently_satisfied"))
                    if not first_ok:
                        if not F635_SUMMARY.exists():
                            recommended_next_target = "F635"
                            recommendation_reason = (
                                "P634 keeps the first clause blocked for the seed-v1 attempt. "
                                "The next honest constructive move is to freeze one explicit future genuinely-new source-object lift/bind target "
                                "using only the seed-v1 materials (F635), without implying admissibility."
                            )
                        elif not P635_SUMMARY.exists():
                            recommended_next_target = "P635"
                            recommendation_reason = (
                                "F635 exports an explicit future lift/bind construction target for a genuinely-new strict-core source object "
                                "above the seed-v1 materials. Next move: run P635 to confirm the constructive move is reduced to that one target."
                            )
                        else:
                            if not F636_SUMMARY.exists():
                                recommended_next_target = "F636"
                                recommendation_reason = (
                                    "P635 reduces the next constructive move to one explicit future target S_sel_int_new_object_target_v1. "
                                    "Next move: export the first strict-core lift/bind attempt instance on that target (F636)."
                                )
                            elif not P636_SUMMARY.exists():
                                recommended_next_target = "P636"
                                recommendation_reason = (
                                    "F636 exports the first lift/bind attempt instance S_sel_int_new_object_lift_bind_attempt_v1. "
                                    "Next move: run P636 to confirm the constructive move is reduced to that one attempt instance."
                                )
                            else:
                                if not F637_SUMMARY.exists():
                                    recommended_next_target = "F637"
                                    recommendation_reason = (
                                        "P636 reduces the next constructive move to the first lift/bind attempt instance S_sel_int_new_object_lift_bind_attempt_v1. "
                                        "Next move: freeze the first realization target above that attempt (F637)."
                                    )
                                elif not P637_SUMMARY.exists():
                                    recommended_next_target = "P637"
                                    recommendation_reason = (
                                        "F637 exports the seed-v1 constructed-source-object realization target. Next move: run P637 to confirm reduction to that target."
                                    )
                                elif not N528_SUMMARY.exists():
                                    recommended_next_target = "N528"
                                    recommendation_reason = (
                                        "P637 confirms the next move is reduced to the seed-v1 realization target; next move: package the reduction theorem (N528)."
                                    )
                                elif not F638_SUMMARY.exists():
                                    recommended_next_target = "F638"
                                    recommendation_reason = (
                                        "N528 fixes the seed-v1 realization target. Next move: export the first realization attempt instance (F638)."
                                    )
                                elif not P638_SUMMARY.exists():
                                    recommended_next_target = "P638"
                                    recommendation_reason = (
                                        "F638 exports the seed-v1 realization attempt instance. Next move: run P638 to confirm reduction to that attempt."
                                    )
                                elif not N529_SUMMARY.exists():
                                    recommended_next_target = "N529"
                                    recommendation_reason = (
                                        "P638 confirms the seed-v1 realization attempt reduction; next move: package the reduction theorem (N529)."
                                    )
                                elif not F639_SUMMARY.exists():
                                    recommended_next_target = "F639"
                                    recommendation_reason = (
                                        "N529 fixes the seed-v1 realization attempt instance. Next move: export the first verdict target (F639)."
                                    )
                                elif not P639_SUMMARY.exists():
                                    recommended_next_target = "P639"
                                    recommendation_reason = (
                                        "F639 exports the seed-v1 realization verdict target. Next move: run P639 to confirm reduction to that verdict target."
                                    )
                                elif not N530_SUMMARY.exists():
                                    recommended_next_target = "N530"
                                    recommendation_reason = (
                                        "P639 confirms the seed-v1 verdict target reduction; next move: package the reduction theorem (N530)."
                                    )
                                elif not F640_SUMMARY.exists():
                                    recommended_next_target = "F640"
                                    recommendation_reason = (
                                        "N530 fixes the seed-v1 verdict target. Next move: refine the two explicit success/failure branches (F640)."
                                    )
                                elif not P640_SUMMARY.exists():
                                    recommended_next_target = "P640"
                                    recommendation_reason = (
                                        "F640 exports explicit success/failure branch names (v1). Next move: run P640 to confirm reduction to the binary split."
                                    )
                                elif not N531_SUMMARY.exists():
                                    recommended_next_target = "N531"
                                    recommendation_reason = (
                                        "P640 confirms the explicit success/failure branch split (v1). Next move: package the split theorem (N531)."
                                    )
                                else:
                                    if not F641_SUMMARY.exists():
                                        recommended_next_target = "F641"
                                        recommendation_reason = (
                                            "N531 packages the explicit success/failure branch split for the seed-v1 realization. "
                                            "Next move (conservative): freeze the failure-first branch ordering (F641)."
                                        )
                                    elif not P641_SUMMARY.exists():
                                        recommended_next_target = "P641"
                                        recommendation_reason = (
                                            "F641 freezes the failure-first branch ordering (v1). "
                                            "Next move: run P641 to test whether an explicit failure verdict is already exported for the fixed realization attempt."
                                        )
                                    elif not N532_SUMMARY.exists():
                                        recommended_next_target = "N532"
                                        recommendation_reason = (
                                            "P641 shows no explicit failure verdict export on current repo state (v1). "
                                            "Next move: package the failure-branch obstruction theorem (N532)."
                                        )
                                    elif not F642_SUMMARY.exists():
                                        recommended_next_target = "F642"
                                        recommendation_reason = (
                                            "N532 packages the current failure-side obstruction (v1). "
                                            "Next move: freeze the remaining success branch as the only remaining branch to attack (F642)."
                                        )
                                    elif not P642_SUMMARY.exists():
                                        recommended_next_target = "P642"
                                        recommendation_reason = (
                                            "F642 freezes the remaining success branch ordering (v1). "
                                            "Next move: run P642 to test whether an explicit success verdict is already exported for the fixed realization attempt."
                                        )
                                    elif not N533_SUMMARY.exists():
                                        recommended_next_target = "N533"
                                        recommendation_reason = (
                                            "P642 shows no explicit success verdict export on current repo state (v1). "
                                            "Next move: package the success-branch obstruction theorem (N533)."
                                        )
                                    else:
                                        recommended_next_target = "EXPORT_EXPLICIT_SUCCESS_OR_FAILURE_VERDICT_FOR_S_SEL_INT_NEW_OBJECT_CONSTRUCTED_REALIZATION_ATTEMPT_V1"
                                        recommendation_reason = (
                                            "N533 packages the current success-branch obstruction (v1) after the failure-side obstruction (N532). "
                                            "Next move: proceed to the post-verdict lower branches (admissibility -> E_orient -> downstream) under no-false-pass ordering for seed v1."
                                        )
                                        if not F643_SUMMARY.exists():
                                            recommended_next_target = "F643"
                                            recommendation_reason = (
                                                "N533 packages the current success-branch obstruction (v1) after the failure-side obstruction (N532). "
                                                "Next move: freeze the admissibility branch as the first lower branch below the exhausted verdict layer (F643)."
                                            )
                                        elif not P643_SUMMARY.exists():
                                            recommended_next_target = "P643"
                                            recommendation_reason = (
                                                "F643 freezes admissibility as the first post-verdict lower branch (seed v1). "
                                                "Next move: run P643 to test whether an explicit admissibility-branch discharge is already exported."
                                            )
                                        elif not N534_SUMMARY.exists():
                                            recommended_next_target = "N534"
                                            recommendation_reason = (
                                                "P643 shows no explicit admissibility-branch discharge export on current repo state (seed v1). "
                                                "Next move: package the admissibility-branch obstruction theorem (N534)."
                                            )
                                        elif not F644_SUMMARY.exists():
                                            recommended_next_target = "F644"
                                            recommendation_reason = (
                                                "N534 packages the current admissibility-branch obstruction (seed v1). "
                                                "Next move: freeze the orientation-export branch as the next lower branch (F644)."
                                            )
                                        elif not P644_SUMMARY.exists():
                                            recommended_next_target = "P644"
                                            recommendation_reason = (
                                                "F644 freezes the orientation-export branch ordering (seed v1). "
                                                "Next move: run P644 to test whether an explicit orientation-export discharge is already exported."
                                            )
                                        elif not N535_SUMMARY.exists():
                                            recommended_next_target = "N535"
                                            recommendation_reason = (
                                                "P644 shows no explicit orientation-export discharge export on current repo state (seed v1). "
                                                "Next move: package the orientation-export obstruction theorem (N535)."
                                            )
                                        elif not F645_SUMMARY.exists():
                                            recommended_next_target = "F645"
                                            recommendation_reason = (
                                                "N535 packages the current orientation-export obstruction (seed v1). "
                                                "Next move: freeze the downstream-completion branch as the last lower branch (F645)."
                                            )
                                        elif not P645_SUMMARY.exists():
                                            recommended_next_target = "P645"
                                            recommendation_reason = (
                                                "F645 freezes the downstream-completion branch ordering (seed v1). "
                                                "Next move: run P645 to test whether an explicit downstream-completion discharge is already exported."
                                            )
                                        elif not N536_SUMMARY.exists():
                                            recommended_next_target = "N536"
                                            recommendation_reason = (
                                                "P645 shows no explicit downstream-completion discharge export on current repo state (seed v1). "
                                                "Next move: package the downstream-completion obstruction theorem (N536)."
                                            )
                                        elif not N537_SUMMARY.exists():
                                            recommended_next_target = "N537"
                                            recommendation_reason = (
                                                "N536 packages the downstream-completion obstruction (seed v1). "
                                                "Next move: package the full post-verdict lower-branch negative closure theorem (N537)."
                                            )
                                        else:
                                            if not F646_SUMMARY.exists():
                                                recommended_next_target = "F646"
                                                recommendation_reason = (
                                                    "N537 closes the post-verdict lower-branch frontier negatively on current repo state (seed v1). "
                                                    "Next move: freeze the minimal strict witness-provider contract + scan signature (F646)."
                                                )
                                            elif not P646_SUMMARY.exists():
                                                recommended_next_target = "P646"
                                                recommendation_reason = (
                                                    "F646 freezes the strict witness-provider contract + scan signature. "
                                                    "Next move: run P646 to mechanically scan whether any exported artifact already matches it."
                                                )
                                            elif not N538_SUMMARY.exists():
                                                recommended_next_target = "N538"
                                                recommendation_reason = (
                                                    "P646 finds no strict witness provider on current repo state. "
                                                    "Next move: package the absence theorem (N538)."
                                                )
                                            else:
                                                recommended_next_target = "IMPLEMENT_STRICT_WITNESS_PROVIDER_EXPORT_PACKET_FOR_SEED_V1_REALIZATION_ATTEMPT"
                                                recommendation_reason = (
                                                    "N538 confirms absence of any strict witness provider matching the F646 signature on current repo state. "
                                                    "Next move: implement one strict witness provider export packet meeting F646 (constructed source object export + admissibility interface), "
                                                    "or prove a scope-level non-bridge impossibility, without implying selector closure / QW-2191 discharge."
                                                )
        except Exception:
            pass

    if recommended_next_target == "H37":
        if N518_THEOREM.exists():
            recommendation_reason += (
                " Note: N518 strengthens N517: any direction-free Aut(Z_12)-invariant reference weight family cannot distinguish sign on the current exported pair1 sine axis "
                "via a scalar of the form Σ_x w(x) u_1(x); therefore H37 requires an explicit reflection-breaking/orientation source or a different observable class."
            )
        elif N517_THEOREM.exists():
            recommendation_reason += (
                " Note: N517 shows even ord-reference weights (ord_Z12 / r_ord) cannot distinguish sign on the current exported pair1 sine axis "
                "via a scalar of the form Σ_x w(x) u_1(x); therefore H37 requires an explicit reflection-breaking/orientation source or a different observable class."
            )
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
                # Probe note is optional hygiene only; do not fail the dashboard if it cannot be parsed.
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
                # Probe note is optional hygiene only; do not fail the dashboard if it cannot be parsed.
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
                # Probe note is optional hygiene only; do not fail the dashboard if it cannot be parsed.
                pass

    artifact = {
        "stage": "P438",
        "goal": "single_dashboard_summary_of_the_T166_diagonal_local_pair1_O2_cut_lane_computability_and_missing_inputs",
        "inputs": {"p432_full_scan": bool(args.full_scan)},
        "dependency_runs": runs,
        "dependency_summaries": {
            "P432": str(P432_SUMMARY.relative_to(REPO)),
            "P444": str(P444_SUMMARY.relative_to(REPO)),
            "P437": str(P437_SUMMARY.relative_to(REPO)),
            "P449": str(P449_SUMMARY.relative_to(REPO)),
            "P434": str(P434_SUMMARY.relative_to(REPO)),
            "P633": (str(P633_SUMMARY.relative_to(REPO)) if P633_SUMMARY.exists() else None),
        },
        "status_snapshot": {
            "P432": p432,
            "P444": p444,
            "P437": p437,
            "P449": p449,
            "P434": p434,
            "P633": (load_json(P633_SUMMARY) if P633_SUMMARY.exists() else None),
        },
        "result": {
            "decision_ready_from_repo_values": decision_ready,
            "any_strict_derived_t168_provider_found": bool(p444.get("any_explicit_strict_derived_t168_provider_found")),
            "P437_computable": str(p437.get("status") or "").startswith("PASS_COMPUTED"),
            "P434_computable": p434.get("status") == "PASS_EVALUATED",
            "P449_all_pairs_cut": bool(p449.get("all_pairs_cut")),
            "t165_selector_ingredient_exported": t165_selector_ingredient_exported,
            "t167_strict_sigma_input_present": t167_strict_sigma_input_present,
            "t166_decision_theorem_exported": t166_decision_theorem_exported,
            "q2191_scoped_discharge_theorem_exported": q2191_scoped_discharge_theorem_exported,
            "recommended_next_strict_target": recommended_next_target,
            "recommendation_reason": recommendation_reason,
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P438",
        "status": "PASS_DASHBOARD_READY",
        "decision_ready_from_repo_values": decision_ready,
        "recommended_next_strict_target": recommended_next_target,
        "recommendation_reason": recommendation_reason,
        "t166_decision_theorem_exported": t166_decision_theorem_exported,
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
