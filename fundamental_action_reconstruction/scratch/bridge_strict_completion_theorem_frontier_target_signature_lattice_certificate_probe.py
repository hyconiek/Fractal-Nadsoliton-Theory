#!/usr/bin/env python3
"""Scratch probe: theorem-frontier target-signature lattice certificate.

The truth-table report enumerates open theorem-atom assignments; the atom-
influence report ranks atom criticality.  This probe projects the same 128
assignments onto the four closure targets and audits the reachable target
signatures.  It proves finite implications such as role-transfer closure
implies selector closure, and ToE closure implies all targets.

This is a target-lattice readiness certificate only.  It exports no theorem
atom and does not close bridge, role-transfer, selector/QW-2191, or ToE.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.md"

SOURCE_REPORTS = {
    "theorem_frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "theorem_frontier_atom_influence": HERE / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.json",
    "theorem_frontier_cut": HERE / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json",
}

TARGETS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def signature_key(values: tuple[bool, ...]) -> str:
    return "".join("1" if value else "0" for value in values)


def signature_from_row(row: dict[str, Any]) -> tuple[bool, ...]:
    return tuple(bool(row[target]) for target in TARGETS)


def true_targets(values: tuple[bool, ...]) -> list[str]:
    return [target for target, value in zip(TARGETS, values) if value]


def assignment_weight(row: dict[str, Any]) -> int:
    return int(row["true_atom_count"])


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    truth_table = loaded["theorem_frontier_truth_table"]
    atom_influence = loaded["theorem_frontier_atom_influence"]
    frontier_cut = loaded["theorem_frontier_cut"]

    all_signature_values = list(itertools.product([False, True], repeat=len(TARGETS)))
    truth_rows = truth_table["truth_table_rows"]
    rows_by_signature: dict[str, list[dict[str, Any]]] = {signature_key(sig): [] for sig in all_signature_values}
    for row in truth_rows:
        rows_by_signature[signature_key(signature_from_row(row))].append(row)

    signature_rows = []
    for sig in all_signature_values:
        key = signature_key(sig)
        rows = rows_by_signature[key]
        reachable = bool(rows)
        min_weight = min((assignment_weight(row) for row in rows), default=None)
        min_sets = [row["true_atoms"] for row in rows if min_weight is not None and assignment_weight(row) == min_weight]
        signature_rows.append(
            {
                "signature": key,
                "target_values": dict(zip(TARGETS, sig)),
                "true_targets": true_targets(sig),
                "reachable": reachable,
                "assignment_count": len(rows),
                "minimal_true_atom_count": min_weight,
                "minimal_true_atom_sets": min_sets,
            }
        )

    reachable_rows = [row for row in signature_rows if row["reachable"]]
    unreachable_rows = [row for row in signature_rows if not row["reachable"]]
    reachable_keys = [row["signature"] for row in reachable_rows]
    unreachable_keys = [row["signature"] for row in unreachable_rows]
    counts_by_signature = {row["signature"]: row["assignment_count"] for row in signature_rows}
    min_weights_by_signature = {row["signature"]: row["minimal_true_atom_count"] for row in reachable_rows}

    implication_rows = [
        {
            "implication": "role_transfer_theorem_level_closure => selector_qw2191_closure",
            "holds_on_all_assignments": all((not row["role_transfer_theorem_level_closure"]) or row["selector_qw2191_closure"] for row in truth_rows),
            "reason": "role-transfer target requires chi11_selector_source, and selector target is exactly chi11_selector_source.",
        },
        {
            "implication": "toe_closure => bridge_theorem_level_closure",
            "holds_on_all_assignments": all((not row["toe_closure"]) or row["bridge_theorem_level_closure"] for row in truth_rows),
            "reason": "ToE target is bridge AND role_transfer AND selector.",
        },
        {
            "implication": "toe_closure => role_transfer_theorem_level_closure",
            "holds_on_all_assignments": all((not row["toe_closure"]) or row["role_transfer_theorem_level_closure"] for row in truth_rows),
            "reason": "ToE target is bridge AND role_transfer AND selector.",
        },
        {
            "implication": "toe_closure => selector_qw2191_closure",
            "holds_on_all_assignments": all((not row["toe_closure"]) or row["selector_qw2191_closure"] for row in truth_rows),
            "reason": "ToE target is bridge AND role_transfer AND selector.",
        },
    ]

    summary = {
        "target_count": len(TARGETS),
        "all_target_signature_count": len(signature_rows),
        "reachable_target_signature_count": len(reachable_rows),
        "unreachable_target_signature_count": len(unreachable_rows),
        "reachable_signatures": reachable_keys,
        "unreachable_signatures": unreachable_keys,
        "counts_by_reachable_signature": {key: counts_by_signature[key] for key in reachable_keys},
        "minimal_weights_by_reachable_signature": min_weights_by_signature,
        "role_implies_selector": implication_rows[0]["holds_on_all_assignments"],
        "toe_implies_all_targets": all(row["holds_on_all_assignments"] for row in implication_rows[1:]),
        "only_full_signature_has_toe_closure": counts_by_signature["1111"] == truth_table["theorem_frontier_truth_table_summary"]["toe_satisfying_assignment_count"] == 1,
        "current_signature_all_false": signature_key(signature_from_row(truth_rows[0])) == "0000" and truth_rows[0]["true_atom_count"] == 0,
        "atom_influence_top_atom_inherited": atom_influence["theorem_frontier_atom_influence_summary"]["top_influence_atoms"] == ["chi11_selector_source"],
        "frontier_cut_open_atoms_match": truth_table["open_atoms"] == list(frontier_cut["open_atoms"].keys()),
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_TARGET_SIGNATURE_LATTICE_CERTIFICATE__NO_CLOSURE_THEOREM",
        "status": "target-signature-lattice-enumerated-six-of-sixteen-signatures-reachable-no-closure",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "target lattice",
                "closure lattice",
                "target signature",
                "reachable target",
                "implication lattice",
            ],
            "finding": "No prior strict-completion report projected the 128 theorem-frontier assignments onto the four closure-target signatures; this report adds that target-lattice audit.",
        },
        "targets": TARGETS,
        "signature_rows": signature_rows,
        "implication_rows": implication_rows,
        "theorem_frontier_target_signature_lattice_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "signature_partition_pass": sum(row["assignment_count"] for row in signature_rows) == truth_table["theorem_frontier_truth_table_summary"]["truth_assignment_count"] == 128,
            "reachable_unreachable_counts_pass": summary["reachable_target_signature_count"] == 6 and summary["unreachable_target_signature_count"] == 10,
            "reachable_signature_counts_pass": summary["counts_by_reachable_signature"] == {"0000": 56, "0010": 49, "0110": 7, "1000": 8, "1010": 7, "1111": 1},
            "minimal_weights_pass": summary["minimal_weights_by_reachable_signature"] == {"0000": 0, "0010": 1, "0110": 4, "1000": 3, "1010": 4, "1111": 7},
            "implications_pass": summary["role_implies_selector"] and summary["toe_implies_all_targets"] and summary["only_full_signature_has_toe_closure"],
            "limits_preserved": summary["current_signature_all_false"] and not summary["bridge_theorem_exported"] and not summary["role_transfer_theorem_exported"] and not summary["selector_closure_exported"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to avoid duplicating a target-signature or implication-lattice audit; none existed for the strict-completion theorem frontier.",
            "projection_step": "Each of the 128 open-atom assignments is projected to a 4-bit target signature ordered as bridge, role-transfer, selector, ToE.",
            "reachability_step": "Exactly 6 of 16 target signatures are reachable: 0000, 0010, 0110, 1000, 1010, and 1111.",
            "count_step": "The reachable signature counts are 56, 49, 7, 8, 7, and 1 respectively, summing to 128.",
            "implication_step": "The lattice proves role-transfer closure implies selector closure, and ToE closure implies bridge, role-transfer, and selector closure.",
            "minimal_step": "Minimal true-atom weights by reachable signature are 0, 1, 4, 3, 4, and 7 respectively; ToE is reachable only at the all-open-atom frontier cut.",
            "scope_step": "This is a finite target-lattice audit only; it exports no theorem atom and proves no bridge, role-transfer, selector, or ToE closure.",
        },
        "hard_limits": [
            "No reachable target signature is promoted to the current theory state.",
            "No missing theorem atom is exported.",
            "No bridge theorem, role-transfer theorem, selector closure, or ToE closure is claimed.",
            "No QW-2191 selector discharge is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict-completion theorem frontier target-signature lattice certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The 128 open-atom assignments project to exactly 6 reachable target",
        "signatures out of 16.  This records closure implications between the",
        "targets without exporting any missing theorem atom.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["theorem_frontier_target_signature_lattice_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Signature rows", ""])
    for row in payload["signature_rows"]:
        lines.append(
            f"- `{row['signature']}`: reachable=`{row['reachable']}`, "
            f"count=`{row['assignment_count']}`, min_weight=`{row['minimal_true_atom_count']}`"
        )
    lines.extend(["", "## Implications", ""])
    for row in payload["implication_rows"]:
        lines.append(f"- `{row['implication']}`: `{row['holds_on_all_assignments']}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Proof certificate", ""])
    for key, value in payload["proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["theorem_frontier_target_signature_lattice_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
