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
    ]

    missing_files: list[str] = []
    for p in (P432_SUMMARY, P444_SUMMARY, P437_SUMMARY, P449_SUMMARY, P434_SUMMARY):
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
        recommended_next_target = "B3"
        recommendation_reason = (
            "T2 theorem-level bridge discharge is now exported (N491). The remaining strict frontier is explicit continuation "
            "under QW-2191 discipline (no implied selector closure): proceed via the B3 topological selector bridge packet "
            "continuation (residual sign lift or explicit sign gauge-irrelevance where required), or explicitly freeze residual "
            "sign as a tracked convention layer while continuing strict-only closure per S2."
        )

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
        },
        "status_snapshot": {
            "P432": p432,
            "P444": p444,
            "P437": p437,
            "P449": p449,
            "P434": p434,
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
