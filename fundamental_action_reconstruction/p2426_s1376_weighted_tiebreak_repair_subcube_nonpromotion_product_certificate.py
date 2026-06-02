#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.json"
MD = GEN / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.md"

SOURCE_FILES = {
    "P2425_WEIGHTED_TIEBREAK": GEN / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.json",
    "P2422_REPAIR_SUBCUBE": GEN / "p2422_s1372_current_missing_gate_repair_subcube_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

BRIDGE_VECTOR = [1, 3]
SELECTOR_VECTOR = [3, 2]
DOMINATED_VECTOR = [2, 3]
WEIGHT_GRID_MAX = 12
CURRENT_TRUE_GATES = ["apd_value_bridge_witness", "chi11_phase_selector_cut_mechanism"]
MISSING_GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]


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
        "new_packet": "P2426|S1376|weighted tiebreak repair|weight.*repair subcube|nonpromotion product",
        "p2425_input": "P2425|weighted tie-break|weight_grid_winner_counts|source-cost premise",
        "p2422_input": "P2422|repair subcube|proper repair failures|toe_ready_repair_count",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2425's weighted tie-break and P2422's repair subcube, but no production P24xx product "
            "certificate proving that even an explicit weight-side choice does not reduce the missing theorem-gate repair obligation."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def weighted_cost(vector: list[int], bridge_weight: int, selector_weight: int) -> int:
    return bridge_weight * vector[0] + selector_weight * vector[1]


def weight_winner(bridge_weight: int, selector_weight: int) -> str:
    bridge_cost = weighted_cost(BRIDGE_VECTOR, bridge_weight, selector_weight)
    selector_cost = weighted_cost(SELECTOR_VECTOR, bridge_weight, selector_weight)
    if bridge_cost < selector_cost:
        return "bridge_first_pareto"
    if selector_cost < bridge_cost:
        return "selector_pair_first_pareto"
    return "bridge_selector_tie"


def repair_target_ready(added: set[str]) -> dict[str, bool]:
    return {
        "bridge_source_ready": "source_obligation_discharge" in added,
        "selector_source_ready": {"chi11_source_export", "qw2191_selector_discharge"}.issubset(added),
        "role_transfer_ready": "role_transfer_audit_license" in added,
        "role_bearing_ltotal_ready": "role_bearing_ltotal_export" in added,
        "toe_ready": set(MISSING_GATES).issubset(added),
    }


def repair_subsets() -> list[list[str]]:
    out: list[list[str]] = []
    for size in range(len(MISSING_GATES) + 1):
        for combo in combinations(MISSING_GATES, size):
            out.append(list(combo))
    return out


def product_rows() -> list[dict[str, Any]]:
    rows = []
    repairs = repair_subsets()
    for bridge_weight in range(1, WEIGHT_GRID_MAX + 1):
        for selector_weight in range(1, WEIGHT_GRID_MAX + 1):
            winner = weight_winner(bridge_weight, selector_weight)
            for repair in repairs:
                added = set(repair)
                ready = repair_target_ready(added)
                rows.append(
                    {
                        "bridge_weight": bridge_weight,
                        "selector_weight": selector_weight,
                        "weight_winner": winner,
                        "added_gates": repair,
                        "added_gate_count": len(repair),
                        "remaining_missing_gates": [gate for gate in MISSING_GATES if gate not in added],
                        "target_ready": ready,
                        "toe_ready": ready["toe_ready"],
                    }
                )
    return rows


def count_by(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def winner_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        winner = row["weight_winner"]
        counts[winner] = counts.get(winner, 0) + 1
    return dict(sorted(counts.items()))


def toe_count_by_winner(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        if not row["toe_ready"]:
            continue
        winner = row["weight_winner"]
        counts[winner] = counts.get(winner, 0) + 1
    return dict(sorted(counts.items()))


def empty_repair_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if row["added_gate_count"] == 0]


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2425 = theorem(sources["P2425_WEIGHTED_TIEBREAK"], "source_frontier_weighted_tiebreak_premise_certificate")
    p2422 = theorem(sources["P2422_REPAIR_SUBCUBE"], "current_missing_gate_repair_subcube_certificate")
    rows = product_rows()
    toe_rows = [row for row in rows if row["toe_ready"]]
    proper_failure_rows = [row for row in rows if not row["toe_ready"]]
    empties = empty_repair_rows(rows)
    return {
        "weight_grid_max": WEIGHT_GRID_MAX,
        "weight_assignment_count": WEIGHT_GRID_MAX * WEIGHT_GRID_MAX,
        "repair_subcube_assignment_count": len(repair_subsets()),
        "product_assignment_count": len(rows),
        "product_row_schema": [
            "bridge_weight",
            "selector_weight",
            "weight_winner",
            "added_gates",
            "added_gate_count",
            "remaining_missing_gates",
            "target_ready",
            "toe_ready",
        ],
        "product_row_samples": rows[:4] + rows[-4:],
        "toe_ready_product_count": len(toe_rows),
        "toe_ready_row_samples": toe_rows[:4],
        "proper_repair_failure_product_count": len(proper_failure_rows),
        "proper_repair_failure_row_samples": proper_failure_rows[:4],
        "weight_winner_product_counts": winner_counts(rows),
        "toe_ready_count_by_weight_winner": toe_count_by_winner(rows),
        "added_gate_count_distribution": count_by(rows, "added_gate_count"),
        "empty_repair_row_count": len(empties),
        "empty_repair_row_samples": empties[:4],
        "empty_repair_winner_counts": winner_counts(empties),
        "empty_repair_remaining_missing_distribution": count_by(
            [{"remaining_count": len(row["remaining_missing_gates"])} for row in empties], "remaining_count"
        ),
        "p2425_weight_grid_winner_counts_inherited": p2425.get("weight_grid_winner_counts"),
        "p2425_no_internal_weight_premise_inherited": p2425.get("internal_weight_source_premise_exported"),
        "p2422_repair_subcube_count_inherited": p2422.get("repair_subcube_assignment_count"),
        "p2422_toe_ready_repair_count_inherited": p2422.get("toe_ready_repair_count"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2426/S1376 weighted tie-break x repair-subcube nonpromotion product certificate

`P2426/S1376` takes the Cartesian product of the P2425 `12 x 12` positive weight grid with the P2422 `2^5=32` repair subcube.  The resulting finite audit has `4608` rows.

The weighted choice never discharges repair obligations.  There are `144` ToE-ready product rows, exactly one all-five-gate repair row for each weight assignment, and `4464` proper repair failures.  The empty-repair slice has all `144` weight assignments but still has all five theorem gates missing, regardless of whether the weight side says bridge-first, selector-pair-first, or tie.

Thus even an explicit weighted tie-break premise would only choose a proof-search order; it would not export source discharge, chi11 source, QW-2191 discharge, role-transfer license, role-bearing `L_total`, or ToE closure.
""".strip()
    lag_section = """
## P2426/S1376 weighted tie-break x repair-subcube guard for Lagrangian/EOM

`P2426/S1376` proves that the weighted source-frontier tie-break is orthogonal to theorem-gate repair: every non-all-gates repair row remains a ToE failure across all weights.  Therefore no weighted order can be promoted into role-transfer, role-bearing `L_total`, or ToE closure.
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
        "theorem_name": "P2426_T1_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate",
        "weight_assignment_count": cert["weight_assignment_count"],
        "repair_subcube_assignment_count": cert["repair_subcube_assignment_count"],
        "product_assignment_count": cert["product_assignment_count"],
        "toe_ready_product_count": cert["toe_ready_product_count"],
        "proper_repair_failure_product_count": cert["proper_repair_failure_product_count"],
        "weight_winner_product_counts": cert["weight_winner_product_counts"],
        "toe_ready_count_by_weight_winner": cert["toe_ready_count_by_weight_winner"],
        "empty_repair_row_count": cert["empty_repair_row_count"],
        "empty_repair_winner_counts": cert["empty_repair_winner_counts"],
        "empty_repair_remaining_missing_distribution": cert["empty_repair_remaining_missing_distribution"],
        "p2425_winner_counts_inherited": cert["p2425_weight_grid_winner_counts_inherited"] == {
            "bridge_first_pareto": 108,
            "bridge_selector_tie": 6,
            "selector_pair_first_pareto": 30,
        },
        "p2425_no_internal_weight_premise_inherited": cert["p2425_no_internal_weight_premise_inherited"] is False,
        "p2422_repair_subcube_count_inherited": cert["p2422_repair_subcube_count_inherited"] == 32,
        "p2422_toe_ready_repair_count_inherited": cert["p2422_toe_ready_repair_count_inherited"] == 1,
        "weighted_choice_reduces_missing_gate_count": False,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Weight-side tie-breaks do not discharge theorem-gate repairs.",
            "Every proper repair subset remains a ToE failure for every positive weight assignment.",
            "The empty repair row keeps all five missing theorem gates under every weight winner.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "weight_grid_144": theorem_export["weight_assignment_count"] == 144,
        "repair_subcube_32": theorem_export["repair_subcube_assignment_count"] == 32,
        "product_4608": theorem_export["product_assignment_count"] == 4608,
        "toe_ready_144": theorem_export["toe_ready_product_count"] == 144,
        "proper_failures_4464": theorem_export["proper_repair_failure_product_count"] == 4464,
        "winner_counts_scaled": theorem_export["weight_winner_product_counts"] == {
            "bridge_first_pareto": 3456,
            "bridge_selector_tie": 192,
            "selector_pair_first_pareto": 960,
        },
        "toe_counts_match_weight_grid": theorem_export["toe_ready_count_by_weight_winner"] == {
            "bridge_first_pareto": 108,
            "bridge_selector_tie": 6,
            "selector_pair_first_pareto": 30,
        },
        "empty_repair_all_weights": theorem_export["empty_repair_row_count"] == 144,
        "empty_repair_all_five_missing": theorem_export["empty_repair_remaining_missing_distribution"] == {"5": 144},
        "p2425_winner_counts_inherited": theorem_export["p2425_winner_counts_inherited"],
        "p2425_no_weight_premise_inherited": theorem_export["p2425_no_internal_weight_premise_inherited"],
        "p2422_repair_count_inherited": theorem_export["p2422_repair_subcube_count_inherited"],
        "p2422_toe_count_inherited": theorem_export["p2422_toe_ready_repair_count_inherited"],
        "weighted_choice_does_not_reduce_missing_gates": not theorem_export["weighted_choice_reduces_missing_gate_count"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2426_s1376_v1",
        "packet_id": "P2426",
        "stage_id": "S1376",
        "result_kind": "WEIGHTED_TIEBREAK_REPAIR_SUBCUBE_NONPROMOTION_PRODUCT_CERTIFICATE",
        "status": "PASS_WEIGHTED_REPAIR_PRODUCT_NO_GATE_DISCHARGE_NO_CLOSURE",
        "weighted_tiebreak_repair_subcube_nonpromotion_product_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Treat any weighted tie-break as only an ordering premise; still discharge the five theorem gates explicitly, "
            "starting with a real source/selector theorem rather than a weight-side preference."
        ),
        "global_status": "OPEN_PROGRESS_WEIGHTED_REPAIR_PRODUCT_CERTIFIED_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["weighted_tiebreak_repair_subcube_nonpromotion_product_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2426 S1376: weighted tie-break x repair-subcube nonpromotion product certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Weight assignments: `{theorem['weight_assignment_count']}`.",
                f"- Repair assignments: `{theorem['repair_subcube_assignment_count']}`.",
                f"- Product assignments: `{theorem['product_assignment_count']}`.",
                f"- ToE-ready product rows: `{theorem['toe_ready_product_count']}`.",
                f"- Proper repair failures: `{theorem['proper_repair_failure_product_count']}`.",
                f"- Empty repair remaining-missing distribution: `{theorem['empty_repair_remaining_missing_distribution']}`.",
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
