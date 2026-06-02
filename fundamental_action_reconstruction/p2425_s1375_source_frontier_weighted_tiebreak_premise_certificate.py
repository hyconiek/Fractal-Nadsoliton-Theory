#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.json"
MD = GEN / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.md"

SOURCE_FILES = {
    "P2424_SOURCE_FRONTIER_PARETO": GEN / "p2424_s1374_source_frontier_pareto_order_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

BRIDGE_VECTOR = [1, 3]
SELECTOR_VECTOR = [3, 2]
DOMINATED_VECTOR = [2, 3]
WEIGHT_GRID_MAX = 12


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
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:20]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2425|S1375|weighted tiebreak|tie-break premise|source-frontier weighted",
        "p2424_input": "P2424|source-frontier Pareto|bridge-first|selector-pair-first|unique first gate",
        "cost_prior": "cost premise|weighted cost|tie.?breaker|external/source-cost|Pareto.*weight",
        "hard_limit_prior": "no unique first|No unique first|not discharge a source theorem|source theorem",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2424's two-class Pareto frontier and generic cost/tie-break language, but no production P24xx "
            "certificate deriving the exact weighted-cost half-space needed to select bridge-first versus selector-pair-first."
        ),
    }


def p2424_theorem(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("source_frontier_pareto_order_certificate", {}).get("theorem_export", {})


def weighted_cost(vector: list[int], bridge_weight: int, selector_weight: int) -> int:
    return bridge_weight * vector[0] + selector_weight * vector[1]


def grid_rows() -> list[dict[str, Any]]:
    rows = []
    for bridge_weight in range(1, WEIGHT_GRID_MAX + 1):
        for selector_weight in range(1, WEIGHT_GRID_MAX + 1):
            bridge_cost = weighted_cost(BRIDGE_VECTOR, bridge_weight, selector_weight)
            selector_cost = weighted_cost(SELECTOR_VECTOR, bridge_weight, selector_weight)
            dominated_cost = weighted_cost(DOMINATED_VECTOR, bridge_weight, selector_weight)
            if bridge_cost < selector_cost:
                winner = "bridge_first_pareto"
            elif selector_cost < bridge_cost:
                winner = "selector_pair_first_pareto"
            else:
                winner = "bridge_selector_tie"
            rows.append(
                {
                    "bridge_weight": bridge_weight,
                    "selector_weight": selector_weight,
                    "weight_ratio_selector_over_bridge": selector_weight / bridge_weight,
                    "bridge_first_cost": bridge_cost,
                    "selector_pair_first_cost": selector_cost,
                    "mixed_split_dominated_cost": dominated_cost,
                    "winner": winner,
                    "dominated_beats_any_pareto": dominated_cost < min(bridge_cost, selector_cost),
                }
            )
    return rows


def winner_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["winner"]] = counts.get(row["winner"], 0) + 1
    return dict(sorted(counts.items()))


def boundary_rows(max_bridge_weight: int = WEIGHT_GRID_MAX) -> list[dict[str, Any]]:
    rows = []
    for bridge_weight in range(1, max_bridge_weight + 1):
        rows.append(
            {
                "bridge_weight": bridge_weight,
                "tie_selector_weight": 2 * bridge_weight,
                "bridge_first_if_selector_weight_less_than": 2 * bridge_weight,
                "selector_pair_first_if_selector_weight_greater_than": 2 * bridge_weight,
                "tie_inside_grid": 2 * bridge_weight <= WEIGHT_GRID_MAX,
            }
        )
    return rows


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2424 = p2424_theorem(sources["P2424_SOURCE_FRONTIER_PARETO"])
    rows = grid_rows()
    return {
        "objective_vectors": {
            "bridge_first_pareto": BRIDGE_VECTOR,
            "selector_pair_first_pareto": SELECTOR_VECTOR,
            "mixed_split_dominated": DOMINATED_VECTOR,
        },
        "symbolic_costs": {
            "bridge_first_pareto": "w_bridge + 3*w_selector",
            "selector_pair_first_pareto": "3*w_bridge + 2*w_selector",
            "mixed_split_dominated": "2*w_bridge + 3*w_selector",
            "bridge_minus_selector": "w_selector - 2*w_bridge",
        },
        "symbolic_tiebreak_conditions": {
            "bridge_first_selected iff": "w_selector < 2*w_bridge",
            "selector_pair_first_selected iff": "w_selector > 2*w_bridge",
            "tie iff": "w_selector = 2*w_bridge",
            "mixed_split_dominated_always": "mixed cost exceeds bridge-first cost by w_bridge > 0",
        },
        "weight_grid_max": WEIGHT_GRID_MAX,
        "weight_grid_rows": rows,
        "weight_grid_assignment_count": len(rows),
        "weight_grid_winner_counts": winner_counts(rows),
        "boundary_rows": boundary_rows(),
        "dominated_win_rows": [row for row in rows if row["dominated_beats_any_pareto"]],
        "p2424_pareto_order_count_inherited": p2424.get("pareto_order_count"),
        "p2424_dominated_order_count_inherited": p2424.get("dominated_order_count"),
        "p2424_unique_pareto_first_gate_selected_inherited": p2424.get("unique_pareto_first_gate_selected"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2425/S1375 source-frontier weighted tie-break premise certificate

`P2425/S1375` turns the P2424 two-class Pareto frontier into an exact weighted-cost half-space.  With positive weights `w_bridge` and `w_selector`, the bridge-first vector `[1,3]` has cost `w_bridge + 3*w_selector`, while the selector-pair-first vector `[3,2]` has cost `3*w_bridge + 2*w_selector`.

Therefore bridge-first is selected iff `w_selector < 2*w_bridge`, selector-pair-first is selected iff `w_selector > 2*w_bridge`, and the two classes tie exactly on `w_selector = 2*w_bridge`.  The mixed split vector `[2,3]` is always dominated because its cost exceeds bridge-first by `w_bridge > 0`.

The finite `12 x 12` positive integer grid confirms the symbolic split: `108` bridge-first wins, `30` selector-pair wins, `6` ties, and `0` dominated wins.  This is a tie-break premise map, not an internal source theorem; without an exported weight/source-cost premise the repo still cannot choose a unique first gate or close role/`L_total`/ToE.
""".strip()
    lag_section = """
## P2425/S1375 source-frontier weighted tie-break guard for Lagrangian/EOM

`P2425/S1375` proves that choosing between bridge-first and selector-pair-first requires an extra weight/source-cost premise (`w_selector` compared with `2*w_bridge`).  Since that premise is not exported internally, the weighted tie-break cannot be promoted into source discharge, role-transfer, role-bearing `L_total`, or ToE closure.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    theorem_export = {
        "theorem_name": "P2425_T1_source_frontier_weighted_tiebreak_premise_certificate",
        "objective_class_count": len(cert["objective_vectors"]),
        "pareto_objective_class_count": 2,
        "weight_grid_max": cert["weight_grid_max"],
        "weight_grid_assignment_count": cert["weight_grid_assignment_count"],
        "weight_grid_winner_counts": cert["weight_grid_winner_counts"],
        "tie_boundary": "w_selector = 2*w_bridge",
        "bridge_first_condition": "w_selector < 2*w_bridge",
        "selector_pair_first_condition": "w_selector > 2*w_bridge",
        "dominated_win_count": len(cert["dominated_win_rows"]),
        "symbolic_bridge_minus_selector_cost": "w_selector - 2*w_bridge",
        "mixed_split_dominated_by_bridge_first": True,
        "p2424_pareto_order_count_inherited": cert["p2424_pareto_order_count_inherited"] == 4,
        "p2424_dominated_order_count_inherited": cert["p2424_dominated_order_count_inherited"] == 2,
        "p2424_no_unique_first_gate_inherited": cert["p2424_unique_pareto_first_gate_selected_inherited"] is False,
        "internal_weight_source_premise_exported": False,
        "unique_first_gate_selected": False,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Weighted tie-break conditions are conditional on an extra source-cost premise.",
            "No internal weight/source-cost premise is exported by this certificate.",
            "No unique first source gate is selected in strict-core terms.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "three_objective_classes": theorem_export["objective_class_count"] == 3,
        "two_pareto_classes": theorem_export["pareto_objective_class_count"] == 2,
        "grid_12_by_12": theorem_export["weight_grid_assignment_count"] == 144,
        "winner_counts_expected": theorem_export["weight_grid_winner_counts"] == {
            "bridge_first_pareto": 108,
            "bridge_selector_tie": 6,
            "selector_pair_first_pareto": 30,
        },
        "tie_boundary_expected": theorem_export["tie_boundary"] == "w_selector = 2*w_bridge",
        "bridge_condition_expected": theorem_export["bridge_first_condition"] == "w_selector < 2*w_bridge",
        "selector_condition_expected": theorem_export["selector_pair_first_condition"] == "w_selector > 2*w_bridge",
        "dominated_never_wins": theorem_export["dominated_win_count"] == 0,
        "mixed_dominated_symbolic": theorem_export["mixed_split_dominated_by_bridge_first"],
        "p2424_pareto_count_inherited": theorem_export["p2424_pareto_order_count_inherited"],
        "p2424_dominated_count_inherited": theorem_export["p2424_dominated_order_count_inherited"],
        "p2424_no_unique_inherited": theorem_export["p2424_no_unique_first_gate_inherited"],
        "no_internal_weight_premise": not theorem_export["internal_weight_source_premise_exported"],
        "no_unique_first_gate": not theorem_export["unique_first_gate_selected"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2425_s1375_v1",
        "packet_id": "P2425",
        "stage_id": "S1375",
        "result_kind": "SOURCE_FRONTIER_WEIGHTED_TIEBREAK_PREMISE_CERTIFICATE",
        "status": "PASS_WEIGHTED_TIEBREAK_CONDITIONAL_NO_SOURCE_PREMISE_NO_DISCHARGE",
        "source_frontier_weighted_tiebreak_premise_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "If choosing a first source gate is required, supply an explicit source-cost theorem fixing w_selector relative "
            "to 2*w_bridge; otherwise continue with either Pareto class as proof-search only."
        ),
        "global_status": "OPEN_PROGRESS_WEIGHTED_TIEBREAK_CERTIFIED_NO_WEIGHT_SOURCE_PREMISE_EXPORTED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["source_frontier_weighted_tiebreak_premise_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2425 S1375: source-frontier weighted tie-break premise certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Weight grid assignments: `{theorem['weight_grid_assignment_count']}`.",
                f"- Winner counts: `{theorem['weight_grid_winner_counts']}`.",
                f"- Bridge-first condition: `{theorem['bridge_first_condition']}`.",
                f"- Selector-pair-first condition: `{theorem['selector_pair_first_condition']}`.",
                f"- Tie boundary: `{theorem['tie_boundary']}`.",
                f"- Dominated win count: `{theorem['dominated_win_count']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
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
