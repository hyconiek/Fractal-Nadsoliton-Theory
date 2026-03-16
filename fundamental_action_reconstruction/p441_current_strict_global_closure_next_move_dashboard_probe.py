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
