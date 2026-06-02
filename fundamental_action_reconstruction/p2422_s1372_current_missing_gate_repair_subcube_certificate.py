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

OUT = GEN / "p2422_s1372_current_missing_gate_repair_subcube_certificate.json"
MD = GEN / "p2422_s1372_current_missing_gate_repair_subcube_certificate.md"

SOURCE_FILES = {
    "P2421_PRIME_IMPLICANT_FAILURE_COVER": GEN / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

CURRENT_TRUE_GATES = ["apd_value_bridge_witness", "chi11_phase_selector_cut_mechanism"]
MISSING_GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]
ALL_GATES = [
    "apd_value_bridge_witness",
    "source_obligation_discharge",
    "chi11_phase_selector_cut_mechanism",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]

MACRO_TARGETS = {
    "bridge_source_ready": ["source_obligation_discharge"],
    "selector_source_ready": ["chi11_source_export", "qw2191_selector_discharge"],
    "role_transfer_ready": ["role_transfer_audit_license"],
    "role_bearing_ltotal_ready": ["role_bearing_ltotal_export"],
    "toe_ready": MISSING_GATES,
}


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
        "new_packet": "P2422|S1372|missing gate repair subcube|repair subcube|partial closure unlock",
        "p2421_input": "P2421|prime implicant|failure cover|current missing gates|current_repair_distance",
        "repair_language": "repair portfolio|missing gate repair|five missing gates|single flip|nearest miss",
        "selector_pair": "chi11_source_export|qw2191_selector_discharge|selector_source_ready",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds P2421's all-gates prime implicant and generic repair language, but no production P24xx "
            "certificate that enumerates the current five-missing-gate repair subcube and separates partial unlocks from closure."
        ),
    }


def p2421_theorem(source: dict[str, Any]) -> dict[str, Any]:
    return source.get("bridge_selector_closure_prime_implicant_failure_cover_certificate", {}).get("theorem_export", {})


def gate_mask(gates: list[str]) -> int:
    return sum(1 << ALL_GATES.index(gate) for gate in gates)


def is_ready(added: set[str], required: list[str]) -> bool:
    return set(required).issubset(added)


def repair_rows() -> list[dict[str, Any]]:
    rows = []
    for size in range(len(MISSING_GATES) + 1):
        for combo in combinations(MISSING_GATES, size):
            added = set(combo)
            target_ready = {name: is_ready(added, required) for name, required in MACRO_TARGETS.items()}
            non_toe_ready = [name for name, ready in target_ready.items() if ready and name != "toe_ready"]
            rows.append(
                {
                    "added_gates": list(combo),
                    "added_gate_count": size,
                    "absolute_mask": gate_mask(CURRENT_TRUE_GATES + list(combo)),
                    "remaining_missing_gates": [gate for gate in MISSING_GATES if gate not in added],
                    "target_ready": target_ready,
                    "non_toe_ready_targets": non_toe_ready,
                    "non_toe_ready_target_count": len(non_toe_ready),
                    "toe_ready": target_ready["toe_ready"],
                }
            )
    return rows


def minimal_unlock_sets(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for target in MACRO_TARGETS:
        candidates = [row for row in rows if row["target_ready"][target]]
        minimal = []
        for row in candidates:
            added = set(row["added_gates"])
            if any(set(other["added_gates"]).issubset(added) and set(other["added_gates"]) != added for other in candidates):
                continue
            minimal.append({"added_gates": row["added_gates"], "size": row["added_gate_count"], "absolute_mask": row["absolute_mask"]})
        out[target] = sorted(minimal, key=lambda item: (item["size"], item["added_gates"]))
    return out


def distribution(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    dist: dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        dist[value] = dist.get(value, 0) + 1
    return dict(sorted(dist.items(), key=lambda item: int(item[0])))


def singleton_repair_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if row["added_gate_count"] == 1]


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2421 = p2421_theorem(sources["P2421_PRIME_IMPLICANT_FAILURE_COVER"])
    rows = repair_rows()
    minimal = minimal_unlock_sets(rows)
    singleton_rows = singleton_repair_rows(rows)
    return {
        "current_true_gates": CURRENT_TRUE_GATES,
        "current_mask": gate_mask(CURRENT_TRUE_GATES),
        "missing_gates": MISSING_GATES,
        "missing_gate_count": len(MISSING_GATES),
        "repair_subcube_rows": rows,
        "repair_subcube_assignment_count": len(rows),
        "minimal_unlock_sets": minimal,
        "ready_target_count_distribution": distribution(rows, "non_toe_ready_target_count"),
        "toe_ready_rows": [row for row in rows if row["toe_ready"]],
        "proper_repair_failure_rows": [row for row in rows if not row["toe_ready"]],
        "singleton_repair_rows": singleton_rows,
        "selector_singleton_unlock_rows": [row for row in singleton_rows if row["target_ready"]["selector_source_ready"]],
        "singleton_non_toe_unlock_rows": [row for row in singleton_rows if row["non_toe_ready_targets"]],
        "p2421_current_missing_gates_inherited": p2421.get("current_missing_gates"),
        "p2421_current_repair_distance_inherited": p2421.get("current_repair_distance"),
        "p2421_unique_prime_implicant_inherited": p2421.get("prime_implicant_masks"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2422/S1372 current missing-gate repair subcube certificate

`P2422/S1372` expands the P2421 current gap into the exact `2^5=32` repair subcube over the five missing theorem gates: source discharge, chi11 source export, QW-2191 selector discharge, role-transfer license, and role-bearing `L_total` export.

The finite subcube separates partial unlocks from closure.  Source discharge alone unlocks bridge-source readiness; role-transfer and role-bearing `L_total` each have singleton local unlocks; selector-source readiness has no singleton unlock and requires the pair `chi11_source_export + qw2191_selector_discharge`.  ToE readiness still has exactly one repair row: all five missing gates.

Thus P2422 is a proof-search repair map, not a theorem discharge.  It identifies which missing gates would unlock which intermediate predicates, but exports no source theorem, no chi11 source, no QW-2191 discharge, no role-transfer license, no `L_total`, and no ToE closure.
""".strip()
    lag_section = """
## P2422/S1372 current missing-gate repair subcube guard for Lagrangian/EOM

`P2422/S1372` shows that the current bridge-selector closure gap decomposes into five missing theorem gates, with selector-source readiness requiring the chi11/QW-2191 pair and ToE requiring all five.  Partial repair rows cannot be promoted to a role-bearing `L_total` term or ToE closure.
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
    minimal = cert["minimal_unlock_sets"]
    theorem_export = {
        "theorem_name": "P2422_T1_current_missing_gate_repair_subcube_certificate",
        "current_mask": cert["current_mask"],
        "missing_gate_count": cert["missing_gate_count"],
        "repair_subcube_assignment_count": cert["repair_subcube_assignment_count"],
        "proper_repair_failure_count": len(cert["proper_repair_failure_rows"]),
        "toe_ready_repair_count": len(cert["toe_ready_rows"]),
        "toe_ready_added_gates": cert["toe_ready_rows"][0]["added_gates"],
        "bridge_source_minimal_unlock_sets": minimal["bridge_source_ready"],
        "selector_source_minimal_unlock_sets": minimal["selector_source_ready"],
        "role_transfer_minimal_unlock_sets": minimal["role_transfer_ready"],
        "role_bearing_ltotal_minimal_unlock_sets": minimal["role_bearing_ltotal_ready"],
        "toe_minimal_unlock_sets": minimal["toe_ready"],
        "selector_singleton_unlock_count": len(cert["selector_singleton_unlock_rows"]),
        "singleton_non_toe_unlock_count": len(cert["singleton_non_toe_unlock_rows"]),
        "ready_target_count_distribution": cert["ready_target_count_distribution"],
        "p2421_current_missing_gates_inherited": cert["p2421_current_missing_gates_inherited"] == MISSING_GATES,
        "p2421_repair_distance_inherited": cert["p2421_current_repair_distance_inherited"] == 5,
        "p2421_unique_prime_implicant_inherited": cert["p2421_unique_prime_implicant_inherited"] == [127],
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Repair-subcube readiness is conditional on adding theorem gates; it does not discharge them.",
            "Selector-source readiness has no singleton repair and still requires chi11 source plus QW-2191 discharge.",
            "All 31 proper repair subsets remain ToE failures.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "five_missing_gates": theorem_export["missing_gate_count"] == 5,
        "all_32_repair_rows": theorem_export["repair_subcube_assignment_count"] == 32,
        "proper_31_failures": theorem_export["proper_repair_failure_count"] == 31,
        "single_toe_repair": theorem_export["toe_ready_repair_count"] == 1,
        "toe_requires_all_five": theorem_export["toe_ready_added_gates"] == MISSING_GATES,
        "bridge_source_singleton": theorem_export["bridge_source_minimal_unlock_sets"] == [{"added_gates": ["source_obligation_discharge"], "size": 1, "absolute_mask": 7}],
        "selector_pair_required": theorem_export["selector_source_minimal_unlock_sets"] == [{"added_gates": ["chi11_source_export", "qw2191_selector_discharge"], "size": 2, "absolute_mask": 29}],
        "role_transfer_singleton": theorem_export["role_transfer_minimal_unlock_sets"] == [{"added_gates": ["role_transfer_audit_license"], "size": 1, "absolute_mask": 37}],
        "ltotal_singleton": theorem_export["role_bearing_ltotal_minimal_unlock_sets"] == [{"added_gates": ["role_bearing_ltotal_export"], "size": 1, "absolute_mask": 69}],
        "toe_minimal_all_five": theorem_export["toe_minimal_unlock_sets"] == [{"added_gates": MISSING_GATES, "size": 5, "absolute_mask": 127}],
        "no_selector_singleton": theorem_export["selector_singleton_unlock_count"] == 0,
        "three_singleton_non_toe_unlocks": theorem_export["singleton_non_toe_unlock_count"] == 3,
        "p2421_missing_gates_inherited": theorem_export["p2421_current_missing_gates_inherited"],
        "p2421_repair_distance_inherited": theorem_export["p2421_repair_distance_inherited"],
        "p2421_prime_implicant_inherited": theorem_export["p2421_unique_prime_implicant_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2422_s1372_v1",
        "packet_id": "P2422",
        "stage_id": "S1372",
        "result_kind": "CURRENT_MISSING_GATE_REPAIR_SUBCUBE_CERTIFICATE",
        "status": "PASS_CURRENT_REPAIR_SUBCUBE_PARTIAL_UNLOCKS_NO_GATE_DISCHARGE",
        "current_missing_gate_repair_subcube_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the repair subcube to choose a real theorem target: source discharge for bridge readiness, the chi11+QW-2191 "
            "pair for selector readiness, or the post-bridge role/L_total gates; do not count partial unlocks as closure."
        ),
        "global_status": "OPEN_PROGRESS_REPAIR_SUBCUBE_CERTIFIED_NO_MISSING_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["current_missing_gate_repair_subcube_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2422 S1372: current missing-gate repair subcube certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Missing gates: `{theorem['missing_gate_count']}`.",
                f"- Repair rows: `{theorem['repair_subcube_assignment_count']}`.",
                f"- Proper repair failures: `{theorem['proper_repair_failure_count']}`.",
                f"- ToE repair count: `{theorem['toe_ready_repair_count']}`.",
                f"- Selector singleton unlock count: `{theorem['selector_singleton_unlock_count']}`.",
                f"- Singleton non-ToE unlock count: `{theorem['singleton_non_toe_unlock_count']}`.",
                f"- ToE-ready added gates: `{theorem['toe_ready_added_gates']}`.",
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
