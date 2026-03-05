#!/usr/bin/env python3
"""
QW-2114: Remaining strict closure roadmap gate.

Purpose:
- formalize closure contract for the current strict-unresolved package IDs,
- stay explicit about what is resolved vs unresolved,
- avoid false closure promotion.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2114_remaining_strict_closure_roadmap.json"
OUT_MD = ROOT / "RAPORT_QW2114_REMAINING_STRICT_CLOSURE_ROADMAP.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def find_entry(entries: List[Dict], pid: str) -> Dict:
    for e in entries:
        if str(e.get("id")) == pid:
            return e
    raise KeyError(f"Missing entry in QW-2069 package: {pid}")


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2071 = load_json("report_qw2071_sm_gr_full_precision_closure_gate.json")
    r2097 = load_json("report_qw2097_ckm_cp_target_refinement_gate.json")
    r2064 = load_json("report_qw2064_micro_derived_renormalization_constants_gate.json")

    entries = r2069["entries"]
    unresolved_ids = list(r2069["coverage"]["strict_unresolved_ids"])

    ckm = find_entry(entries, "delta_cp_ckm")
    grav = find_entry(entries, "gravity_hierarchy_beta20")
    ckm_sel = r2097.get("selected_candidate", {})

    gravity_resolved = bool(
        grav.get("strict_level") == "strict_internal_gate"
        and str(grav.get("status", "")).startswith("derived")
        and bool(grav.get("within_tolerance"))
    )

    roadmap: List[Dict] = []
    if "delta_cp_ckm" in unresolved_ids:
        roadmap.append(
            {
                "id": "delta_cp_ckm",
                "priority_tier": "T1",
                "current_status": ckm.get("status"),
                "strict_level": ckm.get("strict_level"),
                "predicted_value_rad": ckm.get("predicted_value"),
                "reference_value_rad": ckm.get("reference_value"),
                "rel_err_pct": ckm.get("rel_err_pct"),
                "tolerance_rel_pct": ckm.get("tolerance_rel_pct"),
                "selected_candidate": {
                    "matrix_mean_rel_pct": ckm_sel.get("matrix_mean_rel_pct"),
                    "delta_best_rel_err_pct": ckm_sel.get("delta_best_rel_err_pct"),
                    "best_branch": ckm_sel.get("best_branch"),
                },
                "hard_closure_condition": (
                    "Deterministic no-scan/no-retune CKM complex-phase construction must satisfy "
                    "matrix-fidelity and CP-tolerance simultaneously."
                ),
                "required_next_artifacts": [
                    "Deterministic CKM complex-phase law extension pre-registered before execution.",
                    "Independent CP-sensitive observable constraints integrated without kernel retune.",
                    "Strict gate with explicit permutation/convention audit preserving unitarity and matrix fidelity.",
                ],
                "next_gate": "QW_2097_CKM_CP_TARGET_REFINEMENT_GATE",
                "notes": "Current deterministic candidate preserves CKM matrix shape but misses CP tolerance.",
            }
        )

    if "gravity_hierarchy_beta20" in unresolved_ids:
        roadmap.append(
            {
                "id": "gravity_hierarchy_beta20",
                "priority_tier": "T2",
                "current_status": grav.get("status"),
                "strict_level": grav.get("strict_level"),
                "predicted_value": grav.get("predicted_value"),
                "reference_value": grav.get("reference_value"),
                "rel_err_pct": grav.get("rel_err_pct"),
                "tolerance_rel_pct": grav.get("tolerance_rel_pct"),
                "beta_uv_canonical": r2064.get("uv_canonical_constants", {}).get("beta_uv"),
                "hard_closure_condition": (
                    "Promote from model-formula to strict gate only if beta hierarchy is derived through "
                    "a non-circular micro-supported chain with reproducible provenance."
                ),
                "required_next_artifacts": [
                    "Dedicated gravity-hierarchy strict gate with anti-circularity checks.",
                    "Micro-to-bridge derivation trace linked to independently supported observables.",
                ],
                "next_gate": "QW_2115_GRAVITY_HIERARCHY_STRICT_BRIDGE_GATE",
                "notes": "If unresolved here, remains model-formula-level only.",
            }
        )

    flags = {
        "strict_internal_strengthened_pass": bool(r2069.get("strict_internal_strengthened_pass", False)),
        "q2069_q2071_unresolved_consistent": bool(
            sorted(unresolved_ids) == sorted(r2071.get("strict_unresolved_parameters", []))
        ),
        "ckm_target_miss_confirmed_if_unresolved": bool(
            ("delta_cp_ckm" not in unresolved_ids)
            or (str(r2097.get("verdict", "")) == "CKM_CP_TARGET_REFINEMENT_GATE_TARGET_MISS")
        ),
        "gravity_resolved_or_listed_unresolved": bool(
            gravity_resolved or ("gravity_hierarchy_beta20" in unresolved_ids)
        ),
        "roadmap_matches_unresolved_set": bool(
            sorted([str(r.get("id")) for r in roadmap]) == sorted(unresolved_ids)
        ),
        "unresolved_count_nonnegative": bool(len(unresolved_ids) >= 0),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "REMAINING_STRICT_CLOSURE_ROADMAP_READY"
        if all(flags.values())
        else "REMAINING_STRICT_CLOSURE_ROADMAP_INCONSISTENT_INPUTS"
    )

    next_step = (
        "NONE_STRICT_UNRESOLVED_SET_EMPTY"
        if len(unresolved_ids) == 0
        else "EXECUTE_T1_CKM_PHASE_EXTENSION"
        if unresolved_ids == ["delta_cp_ckm"]
        else "EXECUTE_LISTED_ROADMAP_GATES"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2069": "report_qw2069_full_sm_gr_derivation_package.json",
            "q2071": "report_qw2071_sm_gr_full_precision_closure_gate.json",
            "q2097": "report_qw2097_ckm_cp_target_refinement_gate.json",
            "q2064": "report_qw2064_micro_derived_renormalization_constants_gate.json",
        },
        "strict_unresolved_ids": unresolved_ids,
        "gravity_hierarchy_resolution_status": (
            "RESOLVED_STRICT_INTERNAL" if gravity_resolved else "UNRESOLVED"
        ),
        "roadmap": roadmap,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": next_step if verdict.endswith("READY") else "REPAIR_INPUT_CHAIN_AND_RERUN_QW2114",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    title = (
        "REMAINING-0 STRICT CLOSURE ROADMAP"
        if len(unresolved_ids) == 0
        else "REMAINING-1 STRICT CLOSURE ROADMAP"
        if len(unresolved_ids) == 1
        else "REMAINING STRICT CLOSURE ROADMAP"
    )
    lines = [
        f"# RAPORT QW-2114: {title}",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Current Strict-Unresolved Set",
    ]
    if unresolved_ids:
        for pid in unresolved_ids:
            lines.append(f"- `{pid}`")
    else:
        lines.append("- (empty)")

    lines.extend(["", "## Priority Roadmap"])
    if roadmap:
        for row in roadmap:
            lines.extend(
                [
                    f"- `{row['id']}` ({row['priority_tier']}): {row['hard_closure_condition']}",
                    f"  next_gate: `{row['next_gate']}`",
                ]
            )
    else:
        lines.append("- No unresolved IDs to roadmap.")

    lines.extend(["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2114] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2114] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2114] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

