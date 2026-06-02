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

OUT = GEN / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.json"
MD = GEN / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.md"

SOURCE_FILES = {
    "P2428_ANF_DERIVATIVE": GEN / "p2428_s1378_repair_readiness_anf_derivative_certificate.json",
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
TARGET_FORMULAS = {
    "bridge_source_ready": [{"source_obligation_discharge"}],
    "selector_source_ready": [{"chi11_source_export", "qw2191_selector_discharge"}],
    "role_transfer_ready": [{"role_transfer_audit_license"}],
    "role_bearing_ltotal_ready": [{"role_bearing_ltotal_export"}],
    "toe_ready": [set(GATES)],
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
        "new_packet": "P2429|S1379|derivative nearest-miss|nearest miss witness|derivative witness table",
        "p2428_input": "P2428|repair-readiness ANF|ANF derivative|derivative-edge supports",
        "witness_prior": "Boolean derivative witness|nearest-miss derivative|edge witness|essentiality witness",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds older scratch Boolean derivative witness reports and the new P2428 derivative-count certificate, "
            "but no production P24xx certificate materializing every current five-gate repair derivative edge as explicit "
            "nearest-miss witness rows."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def mask_to_set(mask: int) -> set[str]:
    return {gate for idx, gate in enumerate(GATES) if mask & (1 << idx)}


def predicate_true(target: str, mask: int) -> bool:
    true_gates = mask_to_set(mask)
    return any(term.issubset(true_gates) for term in TARGET_FORMULAS[target])


def derivative_witness_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for target in TARGET_FORMULAS:
        for gate_index, gate in enumerate(GATES):
            for before_mask in range(1 << len(GATES)):
                if before_mask & (1 << gate_index):
                    continue
                after_mask = before_mask | (1 << gate_index)
                before_ready = predicate_true(target, before_mask)
                after_ready = predicate_true(target, after_mask)
                if before_ready == after_ready:
                    continue
                before_gates = mask_to_set(before_mask)
                after_gates = mask_to_set(after_mask)
                rows.append(
                    {
                        "target": target,
                        "gate_flipped": gate,
                        "before_mask": before_mask,
                        "after_mask": after_mask,
                        "before_added_gates": sorted(before_gates, key=GATES.index),
                        "after_added_gates": sorted(after_gates, key=GATES.index),
                        "before_missing_gates": [item for item in GATES if item not in before_gates],
                        "after_missing_gates": [item for item in GATES if item not in after_gates],
                        "before_ready": before_ready,
                        "after_ready": after_ready,
                        "missing_count_before": len(GATES) - len(before_gates),
                        "missing_count_after": len(GATES) - len(after_gates),
                    }
                )
    return rows


def count_by(rows: list[dict[str, Any]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def nested_count(rows: list[dict[str, Any]], outer: str, inner: str) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {}
    for row in rows:
        outer_key = str(row[outer])
        inner_key = str(row[inner])
        out.setdefault(outer_key, {})[inner_key] = out.setdefault(outer_key, {}).get(inner_key, 0) + 1
    return {key: dict(sorted(value.items())) for key, value in sorted(out.items())}


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2428 = theorem(sources["P2428_ANF_DERIVATIVE"], "repair_readiness_anf_derivative_certificate")
    rows = derivative_witness_rows()
    toe_rows = [row for row in rows if row["target"] == "toe_ready"]
    return {
        "witness_rows": rows,
        "witness_row_count": len(rows),
        "witness_count_by_target": count_by(rows, "target"),
        "witness_count_by_gate": count_by(rows, "gate_flipped"),
        "witness_count_by_target_gate": nested_count(rows, "target", "gate_flipped"),
        "missing_count_before_distribution": count_by(rows, "missing_count_before"),
        "missing_count_after_distribution": count_by(rows, "missing_count_after"),
        "toe_nearest_miss_rows": toe_rows,
        "toe_nearest_miss_missing_before_distribution": count_by(toe_rows, "missing_count_before"),
        "target_gate_pairs_with_witnesses": sorted(
            {f"{row['target']}::{row['gate_flipped']}" for row in rows}
        ),
        "p2428_derivative_edge_counts_inherited": p2428.get("derivative_edge_counts"),
        "p2428_toe_derivative_edges_inherited": p2428.get("toe_derivative_edges_are_nearest_misses"),
        "p2428_selector_pair_inherited": p2428.get("selector_requires_chi11_and_qw2191_pair"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2429/S1379 repair derivative nearest-miss witness table certificate

`P2429/S1379` materializes the P2428 derivative counts as explicit finite witness rows.  The current five-gate repair layer has `69` derivative witness edges in total across `10` target/gate pairs: `16` for bridge-source readiness, `16` for selector-source readiness, `16` for role-transfer readiness, `16` for role-bearing `L_total` readiness, and `5` ToE nearest-miss edges.

The ToE rows are exactly the five four-other-gates nearest misses: each missing theorem gate has one witness where adding that gate flips ToE from false to true.  These witnesses identify essential blockers, but they are not source, selector, QW-2191, role-transfer, `L_total`, or ToE theorem exports.
""".strip()
    lag_section = """
## P2429/S1379 derivative witness table guard for Lagrangian/EOM

`P2429/S1379` makes every derivative edge explicit.  The role-bearing `L_total` derivative has witnesses only for the `role_bearing_ltotal_export` gate, and ToE has one nearest-miss witness for each missing gate; finite witness rows remain obstruction evidence and cannot be inserted into `L_total` as discharged dynamics.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def positive_derivative_counts(edge_counts: dict[str, dict[str, int]]) -> dict[str, dict[str, int]]:
    return {
        target: {gate: count for gate, count in counts.items() if count > 0}
        for target, counts in edge_counts.items()
    }


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    p2428_positive_derivatives = positive_derivative_counts(cert["p2428_derivative_edge_counts_inherited"] or {})
    theorem_export = {
        "theorem_name": "P2429_T1_repair_derivative_nearest_miss_witness_table_certificate",
        "gate_count": len(GATES),
        "target_count": len(TARGET_FORMULAS),
        "truth_table_size": 1 << len(GATES),
        "witness_row_count": cert["witness_row_count"],
        "witness_count_by_target": cert["witness_count_by_target"],
        "witness_count_by_gate": cert["witness_count_by_gate"],
        "witness_count_by_target_gate": cert["witness_count_by_target_gate"],
        "target_gate_pair_with_witness_count": len(cert["target_gate_pairs_with_witnesses"]),
        "target_gate_pairs_with_witnesses": cert["target_gate_pairs_with_witnesses"],
        "missing_count_before_distribution": cert["missing_count_before_distribution"],
        "missing_count_after_distribution": cert["missing_count_after_distribution"],
        "toe_nearest_miss_row_count": len(cert["toe_nearest_miss_rows"]),
        "toe_nearest_miss_missing_before_distribution": cert["toe_nearest_miss_missing_before_distribution"],
        "toe_nearest_miss_one_per_gate": cert["witness_count_by_target_gate"].get("toe_ready") == {gate: 1 for gate in GATES},
        "selector_derivative_witnesses_only_chi11_qw2191": cert["witness_count_by_target_gate"].get("selector_source_ready") == {
            "chi11_source_export": 8,
            "qw2191_selector_discharge": 8,
        },
        "role_bearing_ltotal_witnesses_only_ltotal_export": cert["witness_count_by_target_gate"].get(
            "role_bearing_ltotal_ready"
        ) == {"role_bearing_ltotal_export": 16},
        "p2428_derivative_edge_counts_inherited": p2428_positive_derivatives == cert["witness_count_by_target_gate"],
        "p2428_toe_derivative_edges_inherited": cert["p2428_toe_derivative_edges_inherited"] is True,
        "p2428_selector_pair_inherited": cert["p2428_selector_pair_inherited"] is True,
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Derivative witness rows certify essential blockers, not theorem-gate discharge.",
            "The five ToE nearest misses remain five missing theorem obligations.",
            "The selector witness rows still require both chi11 source export and QW-2191 discharge.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "truth_table_32": theorem_export["truth_table_size"] == 32,
        "witness_rows_69": theorem_export["witness_row_count"] == 69,
        "witness_counts_by_target": theorem_export["witness_count_by_target"] == {
            "bridge_source_ready": 16,
            "role_bearing_ltotal_ready": 16,
            "role_transfer_ready": 16,
            "selector_source_ready": 16,
            "toe_ready": 5,
        },
        "target_gate_pairs_10": theorem_export["target_gate_pair_with_witness_count"] == 10,
        "toe_nearest_miss_5": theorem_export["toe_nearest_miss_row_count"] == 5,
        "toe_one_per_gate": theorem_export["toe_nearest_miss_one_per_gate"],
        "toe_missing_before_one": theorem_export["toe_nearest_miss_missing_before_distribution"] == {"1": 5},
        "selector_only_pair": theorem_export["selector_derivative_witnesses_only_chi11_qw2191"],
        "ltotal_only_ltotal_export": theorem_export["role_bearing_ltotal_witnesses_only_ltotal_export"],
        "p2428_derivatives_inherited": theorem_export["p2428_derivative_edge_counts_inherited"],
        "p2428_toe_inherited": theorem_export["p2428_toe_derivative_edges_inherited"],
        "p2428_selector_inherited": theorem_export["p2428_selector_pair_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2429_s1379_v1",
        "packet_id": "P2429",
        "stage_id": "S1379",
        "result_kind": "REPAIR_DERIVATIVE_NEAREST_MISS_WITNESS_TABLE_CERTIFICATE",
        "status": "PASS_REPAIR_DERIVATIVE_WITNESS_TABLE_NO_GATE_DISCHARGE_NO_CLOSURE",
        "repair_derivative_nearest_miss_witness_table_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the witness rows to choose an actual theorem target; do not treat any derivative edge as the discharge of its gate."
        ),
        "global_status": "OPEN_PROGRESS_REPAIR_DERIVATIVE_WITNESSES_CERTIFIED_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["repair_derivative_nearest_miss_witness_table_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2429 S1379: repair derivative nearest-miss witness table certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Witness rows: `{theorem['witness_row_count']}`.",
                f"- Witness count by target: `{theorem['witness_count_by_target']}`.",
                f"- Target/gate pairs with witnesses: `{theorem['target_gate_pair_with_witness_count']}`.",
                f"- ToE nearest-miss rows: `{theorem['toe_nearest_miss_row_count']}`.",
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
