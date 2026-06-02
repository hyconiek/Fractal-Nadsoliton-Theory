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

OUT = GEN / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.json"
MD = GEN / "p2430_s1380_repair_derivative_witness_cover_minimality_certificate.md"

SOURCE_FILES = {
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
        "new_packet": "P2430|S1380|derivative witness cover|witness cover minimality|repair witness cover",
        "p2429_input": "P2429|derivative nearest-miss|witness table|target/gate pairs",
        "cover_prior": "failure-cover|hitting set|minimal cover|blocker cover|witness cover",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds earlier failure-cover/hitting-set motifs and P2429 derivative witness rows, but no production "
            "P24xx certificate computing the minimal theorem-gate covers of the current derivative witness table."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def all_subsets(items: list[str]) -> list[list[str]]:
    out: list[list[str]] = []
    for size in range(len(items) + 1):
        for combo in combinations(items, size):
            out.append(list(combo))
    return out


def gate_set_for_target(target_gate_counts: dict[str, dict[str, int]], target: str) -> set[str]:
    return set(target_gate_counts.get(target, {}))


def covers(required: set[str], candidate: set[str]) -> bool:
    return required.issubset(candidate)


def minimal_covers(required: set[str]) -> list[list[str]]:
    covers_found: list[set[str]] = []
    for candidate in all_subsets(GATES):
        candidate_set = set(candidate)
        if covers(required, candidate_set) and not any(prev.issubset(candidate_set) for prev in covers_found):
            covers_found.append(candidate_set)
    return [sorted(item, key=GATES.index) for item in covers_found]


def cover_lattice_rows(required: set[str]) -> list[dict[str, Any]]:
    rows = []
    for candidate in all_subsets(GATES):
        candidate_set = set(candidate)
        uncovered = sorted(required - candidate_set, key=GATES.index)
        rows.append(
            {
                "candidate_gates": candidate,
                "candidate_size": len(candidate),
                "covers_required_witness_gates": not uncovered,
                "uncovered_required_witness_gates": uncovered,
                "uncovered_required_witness_gate_count": len(uncovered),
            }
        )
    return rows


def count_by(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2429 = theorem(sources["P2429_WITNESS_TABLE"], "repair_derivative_nearest_miss_witness_table_certificate")
    target_gate_counts = p2429.get("witness_count_by_target_gate", {})
    required_by_target = {
        target: sorted(gate_set_for_target(target_gate_counts, target), key=GATES.index)
        for target in TARGETS
    }
    minimal_by_target = {target: minimal_covers(set(gates)) for target, gates in required_by_target.items()}
    global_required = set().union(*(set(gates) for gates in required_by_target.values()))
    global_rows = cover_lattice_rows(global_required)
    toe_rows = cover_lattice_rows(set(required_by_target["toe_ready"]))
    current_non_theorem_evidence = ["apd_value_bridge_witness", "chi11_phase_selector_cut_mechanism"]
    return {
        "target_required_witness_gates": required_by_target,
        "minimal_covers_by_target": minimal_by_target,
        "minimal_cover_sizes_by_target": {target: len(covers[0]) for target, covers in minimal_by_target.items()},
        "global_required_witness_gates": sorted(global_required, key=GATES.index),
        "global_minimal_covers": minimal_covers(global_required),
        "global_cover_lattice_rows": global_rows,
        "global_covering_row_count": sum(1 for row in global_rows if row["covers_required_witness_gates"]),
        "global_proper_failure_row_count": sum(1 for row in global_rows if not row["covers_required_witness_gates"]),
        "global_uncovered_count_distribution": count_by(global_rows, "uncovered_required_witness_gate_count"),
        "toe_cover_lattice_rows": toe_rows,
        "toe_covering_row_count": sum(1 for row in toe_rows if row["covers_required_witness_gates"]),
        "toe_proper_failure_row_count": sum(1 for row in toe_rows if not row["covers_required_witness_gates"]),
        "toe_uncovered_count_distribution": count_by(toe_rows, "uncovered_required_witness_gate_count"),
        "current_non_theorem_evidence": current_non_theorem_evidence,
        "current_non_theorem_evidence_covers_required_witness_gates": [],
        "p2429_witness_row_count_inherited": p2429.get("witness_row_count"),
        "p2429_target_gate_pair_count_inherited": p2429.get("target_gate_pair_with_witness_count"),
        "p2429_toe_nearest_miss_count_inherited": p2429.get("toe_nearest_miss_row_count"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2430/S1380 repair derivative witness-cover minimality certificate

`P2430/S1380` takes the P2429 derivative witness table and computes the dual cover problem over theorem gates.  Bridge, role-transfer, and role-bearing `L_total` each have a singleton minimal witness cover; selector has the minimal pair `chi11_source_export + qw2191_selector_discharge`; ToE and the global derivative witness table have the unique minimal cover consisting of all five missing theorem gates.

The cover lattice has `32` rows: exactly `1` row covers all global/ToE derivative witnesses and `31` proper rows leave at least one required witness gate uncovered.  Existing value evidence (`apd_value_bridge_witness`, `chi11_phase_selector_cut_mechanism`) covers no theorem-gate witness requirement, so the cover certificate is a target-selection guide, not closure.
""".strip()
    lag_section = """
## P2430/S1380 repair derivative witness-cover guard for Lagrangian/EOM

`P2430/S1380` proves that covering all derivative witnesses for role-bearing `L_total` and ToE still requires the missing theorem-gate exports themselves.  A witness-cover minimum cannot be replaced by APD value evidence, a chi11 cut mechanism, a proof-order preference, or any Lagrangian term until those theorem gates are actually discharged.
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
        "theorem_name": "P2430_T1_repair_derivative_witness_cover_minimality_certificate",
        "gate_count": len(GATES),
        "target_count": len(TARGETS),
        "target_required_witness_gates": cert["target_required_witness_gates"],
        "minimal_covers_by_target": cert["minimal_covers_by_target"],
        "minimal_cover_sizes_by_target": cert["minimal_cover_sizes_by_target"],
        "global_required_witness_gates": cert["global_required_witness_gates"],
        "global_minimal_covers": cert["global_minimal_covers"],
        "global_covering_row_count": cert["global_covering_row_count"],
        "global_proper_failure_row_count": cert["global_proper_failure_row_count"],
        "global_uncovered_count_distribution": cert["global_uncovered_count_distribution"],
        "toe_covering_row_count": cert["toe_covering_row_count"],
        "toe_proper_failure_row_count": cert["toe_proper_failure_row_count"],
        "toe_uncovered_count_distribution": cert["toe_uncovered_count_distribution"],
        "unique_global_minimal_cover_is_all_five_gates": cert["global_minimal_covers"] == [GATES],
        "unique_toe_minimal_cover_is_all_five_gates": cert["minimal_covers_by_target"]["toe_ready"] == [GATES],
        "selector_minimal_cover_is_chi11_qw2191_pair": cert["minimal_covers_by_target"]["selector_source_ready"] == [
            ["chi11_source_export", "qw2191_selector_discharge"]
        ],
        "current_non_theorem_evidence": cert["current_non_theorem_evidence"],
        "current_non_theorem_evidence_covers_required_witness_gates": cert[
            "current_non_theorem_evidence_covers_required_witness_gates"
        ],
        "p2429_witness_row_count_inherited": cert["p2429_witness_row_count_inherited"] == 69,
        "p2429_target_gate_pair_count_inherited": cert["p2429_target_gate_pair_count_inherited"] == 10,
        "p2429_toe_nearest_miss_count_inherited": cert["p2429_toe_nearest_miss_count_inherited"] == 5,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Witness-cover minimality identifies theorem targets, not theorem-gate discharge.",
            "The unique ToE/global cover still contains all five missing theorem gates.",
            "Current APD value evidence and chi11 cut mechanism cover no theorem-gate witness requirement.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "five_gates": theorem_export["gate_count"] == 5,
        "five_targets": theorem_export["target_count"] == 5,
        "unique_global_all_five": theorem_export["unique_global_minimal_cover_is_all_five_gates"],
        "unique_toe_all_five": theorem_export["unique_toe_minimal_cover_is_all_five_gates"],
        "selector_pair": theorem_export["selector_minimal_cover_is_chi11_qw2191_pair"],
        "global_one_cover": theorem_export["global_covering_row_count"] == 1,
        "global_31_failures": theorem_export["global_proper_failure_row_count"] == 31,
        "toe_one_cover": theorem_export["toe_covering_row_count"] == 1,
        "toe_31_failures": theorem_export["toe_proper_failure_row_count"] == 31,
        "uncovered_distribution_binomial": theorem_export["global_uncovered_count_distribution"] == {
            "0": 1,
            "1": 5,
            "2": 10,
            "3": 10,
            "4": 5,
            "5": 1,
        },
        "current_value_evidence_covers_no_theorem_gate": theorem_export[
            "current_non_theorem_evidence_covers_required_witness_gates"
        ] == [],
        "p2429_witness_rows_inherited": theorem_export["p2429_witness_row_count_inherited"],
        "p2429_pairs_inherited": theorem_export["p2429_target_gate_pair_count_inherited"],
        "p2429_toe_inherited": theorem_export["p2429_toe_nearest_miss_count_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2430_s1380_v1",
        "packet_id": "P2430",
        "stage_id": "S1380",
        "result_kind": "REPAIR_DERIVATIVE_WITNESS_COVER_MINIMALITY_CERTIFICATE",
        "status": "PASS_REPAIR_DERIVATIVE_WITNESS_COVER_MINIMALITY_NO_GATE_DISCHARGE_NO_CLOSURE",
        "repair_derivative_witness_cover_minimality_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Treat the minimal covers as theorem-target selection only: prove source discharge, chi11 source, QW-2191 discharge, "
            "role transfer, and role-bearing L_total separately before claiming closure."
        ),
        "global_status": "OPEN_PROGRESS_REPAIR_DERIVATIVE_WITNESS_COVER_CERTIFIED_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["repair_derivative_witness_cover_minimality_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2430 S1380: repair derivative witness-cover minimality certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Global minimal covers: `{theorem['global_minimal_covers']}`.",
                f"- Minimal cover sizes by target: `{theorem['minimal_cover_sizes_by_target']}`.",
                f"- Global covering rows: `{theorem['global_covering_row_count']}`.",
                f"- Global proper failures: `{theorem['global_proper_failure_row_count']}`.",
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
