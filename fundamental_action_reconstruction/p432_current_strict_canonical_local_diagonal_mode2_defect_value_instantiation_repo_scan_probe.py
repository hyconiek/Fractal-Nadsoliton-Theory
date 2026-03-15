#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe_summary.json"
)


def run_rg(pattern: str, max_count: int = 50) -> dict[str, Any]:
    cmd = [
        "rg",
        "--no-heading",
        "-n",
        "--max-count",
        str(max_count),
        "--glob",
        "!fundamental_action_reconstruction/generated/p432_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_repo_scan_probe*.json",
        pattern,
        str(REPO),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    lines = [line for line in proc.stdout.splitlines() if line.strip()]
    return {
        "pattern": pattern,
        "returncode": proc.returncode,
        "match_count_capped": len(lines),
        "matches_sample": lines,
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    number = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

    # We intentionally scan for *numeric value assignments* (JSON or text form).
    # Pure symbolic occurrences (e.g. in expressions) are not counted as instantiations.
    searches: list[dict[str, str]] = [
        {
            "id": "json_numeric_m2_psi",
            "pattern": rf"\"m2_psi\d+\"\s*:\s*{number}",
        },
        {
            "id": "json_numeric_g4_psi",
            "pattern": rf"\"g4_psi\d+\"\s*:\s*{number}",
        },
        {
            "id": "json_numeric_g6_psi",
            "pattern": rf"\"g6_psi\d+\"\s*:\s*{number}",
        },
        {
            "id": "json_numeric_gY",
            "pattern": rf"\"gY\d+\"\s*:\s*{number}",
        },
        {
            "id": "json_numeric_vpsi",
            "pattern": rf"\"vpsi\d+\"\s*:\s*{number}",
        },
        {
            "id": "json_numeric_vphi",
            "pattern": rf"\"vphi\"\s*:\s*{number}",
        },
        {
            "id": "text_numeric_m2_psi",
            "pattern": rf"\bm2_psi\d+\b\s*=\s*{number}",
        },
        {
            "id": "text_numeric_g4_psi",
            "pattern": rf"\bg4_psi\d+\b\s*=\s*{number}",
        },
        {
            "id": "text_numeric_g6_psi",
            "pattern": rf"\bg6_psi\d+\b\s*=\s*{number}",
        },
        {
            "id": "text_numeric_gY",
            "pattern": rf"\bgY\d+\b\s*=\s*{number}",
        },
        {
            "id": "text_numeric_vpsi",
            "pattern": rf"\bvpsi\d+\b\s*=\s*{number}",
        },
        {
            "id": "text_numeric_vphi",
            "pattern": rf"\bvphi\b\s*=\s*{number}",
        },
        {
            "id": "text_numeric_sigma_opposite_pair_sums",
            "pattern": rf"\bSigma_psi\d+_psi\d+\b\s*[:=]\s*{number}",
        },
        {
            "id": "json_numeric_sigma_opposite_pair_sums",
            "pattern": rf"\"Sigma_psi\d+_psi\d+\"\s*:\s*{number}",
        },
    ]

    results: list[dict[str, Any]] = []
    for row in searches:
        out = run_rg(row["pattern"])
        out["id"] = row["id"]
        results.append(out)

    hits_by_id = {row["id"]: row["match_count_capped"] > 0 for row in results}
    any_numeric_assignment_hits = any(hits_by_id.values())

    # A strict decision of F2(d) from exported values is realistically possible only if the repo already exports
    # numeric instantiations of the six opposite-pair sums (or an equivalent direct diagonal profile). This scan
    # checks for those opposite-pair sums explicitly in both text and JSON forms.
    sigma_sums_numeric_hits = bool(
        hits_by_id["text_numeric_sigma_opposite_pair_sums"] or hits_by_id["json_numeric_sigma_opposite_pair_sums"]
    )
    decision_ready_from_repo_values = sigma_sums_numeric_hits

    artifact = {
        "stage": "P432",
        "goal": "repo_state_scan_for_any_exported_numeric_value_instantiation_sufficient_to_evaluate_the_canonical_FIN_D_local_residual_mode2_defect_F2(d)",
        "searches": results,
        "result": {
            "any_numeric_assignment_hits_found": any_numeric_assignment_hits,
            "sigma_opposite_pair_sum_numeric_instantiation_hits_found": sigma_sums_numeric_hits,
            "decision_ready_from_repo_values": decision_ready_from_repo_values,
            "strict_decision_of_F2_for_canonical_D_local_residual_obtained": False,
        },
        "frontier": "T166_open",
        "no_false_pass": True,
    }

    summary = {
        "stage": "P432",
        "status": "PASS_REPO_SCAN_NO_DECISION_READY_CANONICAL_DIAGONAL_VALUE_EXPORT_FOUND"
        if not decision_ready_from_repo_values
        else "PASS_REPO_SCAN_DECISION_READY_VALUE_EXPORT_HITS_FOUND_REQUIRES_MANUAL_REVIEW",
        "any_numeric_assignment_hits_found": any_numeric_assignment_hits,
        "sigma_opposite_pair_sum_numeric_instantiation_hits_found": sigma_sums_numeric_hits,
        "decision_ready_from_repo_values": decision_ready_from_repo_values,
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
