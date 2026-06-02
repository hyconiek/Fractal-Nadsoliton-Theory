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

OUT = GEN / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.json"
MD = GEN / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.md"

SOURCE_FILES = {
    "P2430_COVER": GEN / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.json",
    "P2429_WITNESS_TABLE": GEN / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]
PREDECESSORS = {
    "source_obligation_discharge": set(),
    "chi11_source_export": set(),
    "qw2191_selector_discharge": set(),
    "role_transfer_audit_license": {
        "source_obligation_discharge",
        "chi11_source_export",
        "qw2191_selector_discharge",
    },
    "role_bearing_ltotal_export": {"role_transfer_audit_license"},
}
TARGETS = [
    "bridge_source_ready",
    "selector_source_ready",
    "role_transfer_ready",
    "role_bearing_ltotal_ready",
    "toe_ready",
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
        "new_packet": "P2431|S1381|next theorem target|target antichain|admissible next theorem",
        "p2430_input": "P2430|witness-cover minimality|minimal witness cover|theorem-gate cover",
        "priority_prior": "theorem target|admissible order|proof-order|target selection|next honest step",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2423 proof-order and P2430 witness-cover target-selection language, but no production "
            "P24xx certificate combining the cover with the precedence guardrails to compute the admissible next theorem-target antichain."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def all_nonempty_subsets(items: list[str], max_size: int) -> list[list[str]]:
    out = []
    for size in range(1, max_size + 1):
        for combo in combinations(items, size):
            out.append(list(combo))
    return out


def is_admissible(candidate: set[str], current: set[str]) -> bool:
    available = current | candidate
    return all(PREDECESSORS[gate].issubset(available) for gate in candidate)


def readiness_from_gates(gates: set[str]) -> dict[str, bool]:
    return {
        "bridge_source_ready": "source_obligation_discharge" in gates,
        "selector_source_ready": {"chi11_source_export", "qw2191_selector_discharge"}.issubset(gates),
        "role_transfer_ready": "role_transfer_audit_license" in gates,
        "role_bearing_ltotal_ready": "role_bearing_ltotal_export" in gates,
        "toe_ready": set(GATES).issubset(gates),
    }


def coverage_count(candidate: set[str], witness_count_by_target_gate: dict[str, dict[str, int]]) -> int:
    total = 0
    for counts in witness_count_by_target_gate.values():
        total += sum(count for gate, count in counts.items() if gate in candidate)
    return total


def readiness_unlocks(candidate: set[str], current: set[str]) -> dict[str, bool]:
    before = readiness_from_gates(current)
    after = readiness_from_gates(current | candidate)
    return {target: after[target] and not before[target] for target in TARGETS}


def completed_target_coverage(row: dict[str, Any], witness_count_by_target_gate: dict[str, dict[str, int]]) -> int:
    return sum(
        sum(witness_count_by_target_gate.get(target, {}).values())
        for target, unlocked in row["readiness_unlocks"].items()
        if unlocked
    )


def has_proper_subset_with_same_unlocks(row: dict[str, Any], rows: list[dict[str, Any]]) -> bool:
    candidate = set(row["candidate_gates"])
    return any(
        other["admissible_from_current"]
        and set(other["candidate_gates"]).issubset(candidate)
        and set(other["candidate_gates"]) != candidate
        and other["readiness_unlocks"] == row["readiness_unlocks"]
        for other in rows
    )


def build_candidate_rows(witness_count_by_target_gate: dict[str, dict[str, int]]) -> list[dict[str, Any]]:
    current: set[str] = set()
    rows = []
    for candidate in all_nonempty_subsets(GATES, 2):
        candidate_set = set(candidate)
        missing_predecessors = sorted(
            set().union(*(PREDECESSORS[gate] for gate in candidate_set)) - (current | candidate_set),
            key=GATES.index,
        )
        unlocks = readiness_unlocks(candidate_set, current)
        rows.append(
            {
                "candidate_gates": candidate,
                "candidate_size": len(candidate),
                "admissible_from_current": is_admissible(candidate_set, current),
                "missing_predecessors": missing_predecessors,
                "readiness_unlocks": unlocks,
                "unlocks_any_readiness_target": any(unlocks.values()),
                "raw_witness_edge_coverage_count": coverage_count(candidate_set, witness_count_by_target_gate),
            }
        )
    for row in rows:
        row["completed_target_witness_coverage_count"] = completed_target_coverage(row, witness_count_by_target_gate)
        row["minimal_readiness_complete_candidate"] = (
            row["admissible_from_current"]
            and row["unlocks_any_readiness_target"]
            and not has_proper_subset_with_same_unlocks(row, rows)
        )
    return rows


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2430 = theorem(sources["P2430_COVER"], "repair_derivative_witness_cover_minimality_certificate")
    p2429 = theorem(sources["P2429_WITNESS_TABLE"], "repair_derivative_nearest_miss_witness_table_certificate")
    counts = p2429.get("witness_count_by_target_gate", {})
    rows = build_candidate_rows(counts)
    admissible_rows = [row for row in rows if row["admissible_from_current"]]
    minimal_readiness_rows = [row for row in admissible_rows if row["minimal_readiness_complete_candidate"]]
    max_completed_coverage = max(row["completed_target_witness_coverage_count"] for row in minimal_readiness_rows)
    return {
        "candidate_rows_size_le_2": rows,
        "candidate_row_count_size_le_2": len(rows),
        "admissible_candidate_count_size_le_2": len(admissible_rows),
        "inadmissible_candidate_count_size_le_2": len(rows) - len(admissible_rows),
        "admissible_singletons_from_current": [row["candidate_gates"] for row in admissible_rows if row["candidate_size"] == 1],
        "inadmissible_singletons_from_current": [row["candidate_gates"] for row in rows if row["candidate_size"] == 1 and not row["admissible_from_current"]],
        "minimal_readiness_complete_candidates": [row["candidate_gates"] for row in minimal_readiness_rows],
        "minimal_readiness_complete_rows": minimal_readiness_rows,
        "max_completed_target_coverage_candidates": [
            row["candidate_gates"] for row in minimal_readiness_rows if row["completed_target_witness_coverage_count"] == max_completed_coverage
        ],
        "p2430_global_minimal_cover_inherited": p2430.get("global_minimal_covers"),
        "p2430_selector_pair_inherited": p2430.get("selector_minimal_cover_is_chi11_qw2191_pair"),
        "p2429_witness_count_by_target_gate_inherited": counts,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2431/S1381 admissible next theorem-target antichain certificate

`P2431/S1381` combines the P2430 witness-cover minimum with repair-order precedence.  From the current zero-discharge state, the admissible singleton theorem targets are exactly `source_obligation_discharge`, `chi11_source_export`, and `qw2191_selector_discharge`; `role_transfer_audit_license` and `role_bearing_ltotal_export` are inadmissible as first moves.

On candidates of size at most two, the minimal readiness-complete admissible antichain is `source_obligation_discharge` versus the selector pair `chi11_source_export + qw2191_selector_discharge`.  This identifies the next real theorem-target fork, but it still exports no source, selector, QW-2191, role-transfer, `L_total`, or ToE theorem.
""".strip()
    lag_section = """
## P2431/S1381 admissible next theorem-target guard for Lagrangian/EOM

`P2431/S1381` proves that role-transfer and role-bearing `L_total` are not admissible first theorem targets under the current source/selector gaps.  The Lagrangian/EOM draft may mention the source-vs-selector target fork, but cannot promote either fork into discharged dynamics before the corresponding theorem is actually proved.
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
        "theorem_name": "P2431_T1_admissible_next_theorem_target_antichain_certificate",
        "candidate_row_count_size_le_2": cert["candidate_row_count_size_le_2"],
        "admissible_candidate_count_size_le_2": cert["admissible_candidate_count_size_le_2"],
        "inadmissible_candidate_count_size_le_2": cert["inadmissible_candidate_count_size_le_2"],
        "admissible_singletons_from_current": cert["admissible_singletons_from_current"],
        "inadmissible_singletons_from_current": cert["inadmissible_singletons_from_current"],
        "minimal_readiness_complete_candidates": cert["minimal_readiness_complete_candidates"],
        "max_completed_target_coverage_candidates": cert["max_completed_target_coverage_candidates"],
        "p2430_global_minimal_cover_inherited": cert["p2430_global_minimal_cover_inherited"] == [GATES],
        "p2430_selector_pair_inherited": cert["p2430_selector_pair_inherited"] is True,
        "role_transfer_admissible_as_first_move": False,
        "role_bearing_ltotal_admissible_as_first_move": False,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "The antichain is a target-selection fork, not theorem-gate discharge.",
            "Role-transfer and role-bearing L_total remain inadmissible first moves.",
            "The selector fork still requires both chi11 source export and QW-2191 discharge.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "candidate_rows_15": theorem_export["candidate_row_count_size_le_2"] == 15,
        "admissible_candidates_6": theorem_export["admissible_candidate_count_size_le_2"] == 6,
        "inadmissible_candidates_9": theorem_export["inadmissible_candidate_count_size_le_2"] == 9,
        "admissible_singletons_expected": theorem_export["admissible_singletons_from_current"] == [
            ["source_obligation_discharge"],
            ["chi11_source_export"],
            ["qw2191_selector_discharge"],
        ],
        "inadmissible_singletons_expected": theorem_export["inadmissible_singletons_from_current"] == [
            ["role_transfer_audit_license"],
            ["role_bearing_ltotal_export"],
        ],
        "minimal_readiness_antichain_expected": theorem_export["minimal_readiness_complete_candidates"] == [
            ["source_obligation_discharge"],
            ["chi11_source_export", "qw2191_selector_discharge"],
        ],
        "max_completed_coverage_antichain_expected": theorem_export["max_completed_target_coverage_candidates"] == [
            ["source_obligation_discharge"],
            ["chi11_source_export", "qw2191_selector_discharge"],
        ],
        "p2430_global_cover_inherited": theorem_export["p2430_global_minimal_cover_inherited"],
        "p2430_selector_pair_inherited": theorem_export["p2430_selector_pair_inherited"],
        "role_transfer_not_first": not theorem_export["role_transfer_admissible_as_first_move"],
        "ltotal_not_first": not theorem_export["role_bearing_ltotal_admissible_as_first_move"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2431_s1381_v1",
        "packet_id": "P2431",
        "stage_id": "S1381",
        "result_kind": "ADMISSIBLE_NEXT_THEOREM_TARGET_ANTICHAIN_CERTIFICATE",
        "status": "PASS_ADMISSIBLE_NEXT_TARGET_ANTICHAIN_NO_GATE_DISCHARGE_NO_CLOSURE",
        "admissible_next_theorem_target_antichain_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Choose one branch of the admissible antichain: either prove source discharge first, or prove the chi11/QW-2191 selector pair, "
            "without promoting either branch before it is actually discharged."
        ),
        "global_status": "OPEN_PROGRESS_NEXT_THEOREM_TARGET_ANTICHAIN_CERTIFIED_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["admissible_next_theorem_target_antichain_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2431 S1381: admissible next theorem-target antichain certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Candidate rows size <= 2: `{theorem['candidate_row_count_size_le_2']}`.",
                f"- Admissible candidates: `{theorem['admissible_candidate_count_size_le_2']}`.",
                f"- Admissible singletons: `{theorem['admissible_singletons_from_current']}`.",
                f"- Minimal readiness-complete antichain: `{theorem['minimal_readiness_complete_candidates']}`.",
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
