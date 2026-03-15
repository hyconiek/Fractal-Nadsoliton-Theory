#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

QW2122 = REPO / "report_qw2122_psi_potential_diagonal_floor_gate.json"
QW2123 = REPO / "report_qw2123_vacuum_branch_selection_strict_gate.json"
QW2124 = REPO / "report_qw2124_scalar_vacuum_closure_branch_resolved_gate.json"

P437_IN = GENERATED / "p437_input_vpsi_g4_g6_candidate.json"

OUT_JSON = GENERATED / "p442_current_strict_t168_vacuum_and_self_coupling_provider_frontier_audit_probe.json"
OUT_SUMMARY = (
    GENERATED / "p442_current_strict_t168_vacuum_and_self_coupling_provider_frontier_audit_probe_summary.json"
)


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing_files: list[str] = []
    for p in (QW2122, QW2123, QW2124, P437_IN):
        if not p.exists():
            missing_files.append(str(p.relative_to(REPO)) if p.is_absolute() else str(p))

    if missing_files:
        summary = {
            "stage": "P442",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing_files,
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    r2122 = load_json(QW2122)
    r2123 = load_json(QW2123)
    r2124 = load_json(QW2124)
    p437_in = load_json(P437_IN)

    # Strict scalar-side values (available in strict chain).
    lambda_psi_strict = (r2122.get("inputs") or {}).get("lambda_psi_strict")
    rho_star_sq = ((r2122.get("branch_results") or {}).get("broken_branch_strict") or {}).get("rho_star_sq")
    diag_floor = ((r2122.get("branch_results") or {}).get("broken_branch_strict") or {}).get("diag_floor")

    branch_selection_verdict = str(r2123.get("verdict", ""))
    branch_selection_result = str((r2123.get("selection_rule") or {}).get("result", ""))
    scalar_closure_verdict = str(r2124.get("verdict", ""))

    scalar_sources_ok = {
        "qw2122_lambda_psi_strict_present": is_number(lambda_psi_strict),
        "qw2122_broken_branch_rho_star_sq_present": is_number(rho_star_sq) and float(rho_star_sq) > 0.0,
        "qw2122_broken_branch_diag_floor_present": is_number(diag_floor),
        "qw2123_branch_selection_pass_broken_branch_required": branch_selection_verdict.startswith(
            "VACUUM_BRANCH_SELECTION_STRICT_GATE_PASS"
        )
        and branch_selection_result == "broken_branch_required",
        "qw2124_scalar_vacuum_closure_branch_resolved_strict_pass": scalar_closure_verdict.startswith(
            "SCALAR_VACUUM_CLOSURE_BRANCH_RESOLVED_STRICT_PASS"
        ),
    }

    # T168 still needs a strict-derived canonical provider for vpsi[0..11], g4[0..11], g6[0..11].
    # We only audit the designated harness input for now (values are intentionally null today).
    def has_numeric_len12(name: str) -> bool:
        arr = p437_in.get(name)
        return bool(
            isinstance(arr, list)
            and len(arr) == 12
            and all(is_number(x) for x in arr)
            and all(float(x) != 0.0 for x in arr)  # N477 division premise
        )

    have_vpsi = has_numeric_len12("vpsi")
    have_g4 = has_numeric_len12("g4")
    have_g6 = has_numeric_len12("g6")

    missing: list[str] = []
    if not have_vpsi:
        missing.append("strict-derived numeric vpsi[0..11] provider (length-12, nonzero)")
    if not have_g4:
        missing.append("strict-derived numeric g4[0..11] provider (length-12)")
    if not have_g6:
        missing.append("strict-derived numeric g6[0..11] provider (length-12)")

    # Even with strict scalar-side values, the lift to per-site arrays is not exported today without mapping premises.
    mapping_gap = True
    recommended_next = "T168"
    recommendation_reason = (
        "P437 requires strict-derived (vpsi,g4,g6) arrays; repo currently exports only strict scalar-side psi potential "
        "values (QW-2122/2123/2124) and no strict lift/mapping to the canonical per-site families."
    )

    artifact = {
        "stage": "P442",
        "goal": "audit_frontier_for_T168_vacuum_and_self_coupling_value_provider_needed_by_P437_on_the_T166_diagonal_lane",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "strict_scalar_sources": {
            "QW-2122": str(QW2122.relative_to(REPO)),
            "QW-2123": str(QW2123.relative_to(REPO)),
            "QW-2124": str(QW2124.relative_to(REPO)),
        },
        "designated_harness_input": str(P437_IN.relative_to(REPO)),
        "scalar_values": {
            "lambda_psi_strict": float(lambda_psi_strict) if is_number(lambda_psi_strict) else None,
            "broken_branch_rho_star_sq": float(rho_star_sq) if is_number(rho_star_sq) else None,
            "broken_branch_diag_floor": float(diag_floor) if is_number(diag_floor) else None,
            "qw2123_verdict": branch_selection_verdict,
            "qw2123_selection_rule_result": branch_selection_result,
            "qw2124_verdict": scalar_closure_verdict,
        },
        "audit": {
            "scalar_sources_ok": scalar_sources_ok,
            "p437_input_has_numeric_arrays": {"vpsi": have_vpsi, "g4": have_g4, "g6": have_g6},
            "mapping_gap_remains": mapping_gap,
            "missing": missing,
        },
        "result": {
            "t168_discharged": False,
            "recommended_next_strict_target": recommended_next,
            "recommendation_reason": recommendation_reason,
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P442",
        "status": "PASS_AUDIT_READY_T168_NOT_DISCHARGED",
        "scalar_sources_ok_all": bool(all(scalar_sources_ok.values())),
        "mapping_gap_remains": mapping_gap,
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

