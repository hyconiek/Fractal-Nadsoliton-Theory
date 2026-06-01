#!/usr/bin/env python3
"""Scratch probe: theorem-frontier low-weight extension certificate.

The target-signature lattice says which target signatures are reachable at all.
This probe audits the first finite neighborhoods around the current all-false
frontier state: all singleton and pair extensions of the seven open theorem
atoms.  It tells us exactly what can and cannot be unlocked by one or two new
theorem atoms.

This is a planning/readiness certificate only.  It does not export any theorem
atom, and it does not close the bridge, role-transfer, selector/QW-2191, or ToE
targets in the current theory state.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.md"

SOURCE_REPORTS = {
    "theorem_frontier_truth_table": HERE / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "theorem_frontier_target_signature_lattice": HERE / "bridge_strict_completion_theorem_frontier_target_signature_lattice_certificate_report.json",
    "theorem_frontier_atom_influence": HERE / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.json",
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


def signature_key(row: dict[str, Any]) -> str:
    return "".join("1" if row[target] else "0" for target in TARGETS)


def target_values(assignment: dict[str, bool], target_defs: dict[str, Any]) -> dict[str, bool]:
    bridge = all(assignment[atom] for atom in target_defs["bridge_theorem_level_closure"])
    role = all(assignment[atom] for atom in target_defs["role_transfer_theorem_level_closure"])
    selector = all(assignment[atom] for atom in target_defs["selector_qw2191_closure"])
    toe = bridge and role and selector
    return {
        "bridge_theorem_level_closure": bridge,
        "role_transfer_theorem_level_closure": role,
        "selector_qw2191_closure": selector,
        "toe_closure": toe,
    }


def extension_row(atoms: list[str], target_defs: dict[str, Any], true_atoms: tuple[str, ...]) -> dict[str, Any]:
    assignment = {atom: atom in true_atoms for atom in atoms}
    values = target_values(assignment, target_defs)
    signature = "".join("1" if values[target] else "0" for target in TARGETS)
    return {
        "true_atoms": list(true_atoms),
        "weight": len(true_atoms),
        "target_signature": signature,
        "target_values": values,
        "newly_closed_targets": [target for target in TARGETS if values[target]],
        "exports_now": False,
    }


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    truth_table = loaded["theorem_frontier_truth_table"]
    target_lattice = loaded["theorem_frontier_target_signature_lattice"]
    atom_influence = loaded["theorem_frontier_atom_influence"]

    atoms = truth_table["open_atoms"]
    target_defs = truth_table["target_definitions"]
    singleton_rows = [extension_row(atoms, target_defs, combo) for combo in itertools.combinations(atoms, 1)]
    pair_rows = [extension_row(atoms, target_defs, combo) for combo in itertools.combinations(atoms, 2)]
    low_weight_rows = singleton_rows + pair_rows

    singleton_unlock_rows = [row for row in singleton_rows if row["newly_closed_targets"]]
    pair_unlock_rows = [row for row in pair_rows if row["newly_closed_targets"]]
    pair_bridge_rows = [row for row in pair_rows if row["target_values"]["bridge_theorem_level_closure"]]
    pair_role_rows = [row for row in pair_rows if row["target_values"]["role_transfer_theorem_level_closure"]]
    pair_toe_rows = [row for row in pair_rows if row["target_values"]["toe_closure"]]

    summary = {
        "open_atom_count": len(atoms),
        "singleton_extension_count": len(singleton_rows),
        "pair_extension_count": len(pair_rows),
        "low_weight_extension_count": len(low_weight_rows),
        "singleton_unlock_count": len(singleton_unlock_rows),
        "singleton_unlock_atoms": [row["true_atoms"][0] for row in singleton_unlock_rows],
        "singleton_unlock_signatures": sorted({row["target_signature"] for row in singleton_unlock_rows}),
        "chi11_is_only_singleton_unlock": [row["true_atoms"] for row in singleton_unlock_rows] == [["chi11_selector_source"]],
        "pair_unlock_count": len(pair_unlock_rows),
        "pair_unlocks_are_selector_only": all(row["target_signature"] == "0010" for row in pair_unlock_rows),
        "pair_unlocks_all_contain_chi11": all("chi11_selector_source" in row["true_atoms"] for row in pair_unlock_rows),
        "no_singleton_closes_bridge_role_or_toe": all(not row["target_values"]["bridge_theorem_level_closure"] and not row["target_values"]["role_transfer_theorem_level_closure"] and not row["target_values"]["toe_closure"] for row in singleton_rows),
        "no_pair_closes_bridge": len(pair_bridge_rows) == 0,
        "no_pair_closes_role_transfer": len(pair_role_rows) == 0,
        "no_pair_closes_toe": len(pair_toe_rows) == 0,
        "target_lattice_min_weights_inherited": target_lattice["theorem_frontier_target_signature_lattice_summary"]["minimal_weights_by_reachable_signature"] == {"0000": 0, "0010": 1, "0110": 4, "1000": 3, "1010": 4, "1111": 7},
        "atom_influence_top_atom_inherited": atom_influence["theorem_frontier_atom_influence_summary"]["top_influence_atoms"] == ["chi11_selector_source"],
        "current_signature_all_false": signature_key(truth_table["truth_table_rows"][0]) == "0000" and truth_table["truth_table_rows"][0]["true_atom_count"] == 0,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_THEOREM_FRONTIER_LOW_WEIGHT_EXTENSION_CERTIFICATE__NO_THEOREM_EXPORT",
        "status": "singleton-and-pair-frontier-extensions-enumerated-selector-only-low-weight-unlock-no-closure",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "singleton extension",
                "one-atom",
                "next theorem atom",
                "frontier extension",
                "marginal unlock",
                "unlock matrix",
            ],
            "finding": "No prior strict-completion report enumerated the singleton/pair extension neighborhood of the current all-false theorem frontier; this report adds that low-weight extension audit.",
        },
        "singleton_extension_rows": singleton_rows,
        "pair_extension_rows": pair_rows,
        "theorem_frontier_low_weight_extension_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "enumeration_counts_pass": summary["singleton_extension_count"] == 7 and summary["pair_extension_count"] == 21 and summary["low_weight_extension_count"] == 28,
            "singleton_unlock_pass": summary["singleton_unlock_count"] == 1 and summary["chi11_is_only_singleton_unlock"] and summary["singleton_unlock_signatures"] == ["0010"],
            "pair_unlock_pass": summary["pair_unlock_count"] == 6 and summary["pair_unlocks_are_selector_only"] and summary["pair_unlocks_all_contain_chi11"],
            "low_weight_blockers_pass": summary["no_singleton_closes_bridge_role_or_toe"] and summary["no_pair_closes_bridge"] and summary["no_pair_closes_role_transfer"] and summary["no_pair_closes_toe"],
            "inherited_reports_pass": summary["target_lattice_min_weights_inherited"] and summary["atom_influence_top_atom_inherited"],
            "limits_preserved": summary["current_signature_all_false"] and not summary["bridge_theorem_exported"] and not summary["role_transfer_theorem_exported"] and not summary["selector_closure_exported"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to avoid duplicating a singleton/pair frontier-extension audit; none existed for the strict-completion theorem frontier.",
            "enumeration_step": "All 7 singleton and all 21 pair extensions of the seven open atoms are enumerated from the current all-false state.",
            "singleton_step": "The only singleton extension that closes any target is chi11_selector_source, and it closes only selector/QW-2191 target signature 0010.",
            "pair_step": "Exactly six pair extensions close any target; each contains chi11_selector_source and still closes only selector signature 0010.",
            "blocker_step": "No singleton or pair extension closes bridge theorem, role-transfer theorem, or ToE; bridge still needs three strict-source atoms, role transfer needs four role/selector atoms, and ToE needs all seven.",
            "scope_step": "This is a low-weight planning certificate only; it exports no theorem atom and proves no bridge, role-transfer, selector, or ToE closure.",
        },
        "hard_limits": [
            "No singleton or pair extension is promoted to the current theory state.",
            "No missing theorem atom is exported.",
            "No bridge theorem, role-transfer theorem, selector closure, or ToE closure is claimed.",
            "No QW-2191 selector discharge is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict-completion theorem frontier low-weight extension certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "All singleton and pair extensions of the current all-false theorem frontier",
        "are enumerated.  Low-weight unlocks are selector-only and do not close the",
        "bridge, role-transfer, or ToE targets.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["theorem_frontier_low_weight_extension_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Singleton rows", ""])
    for row in payload["singleton_extension_rows"]:
        lines.append(f"- `{row['true_atoms']}` -> `{row['target_signature']}` closes `{row['newly_closed_targets']}`")
    lines.extend(["", "## Pair unlock rows", ""])
    for row in payload["pair_extension_rows"]:
        if row["newly_closed_targets"]:
            lines.append(f"- `{row['true_atoms']}` -> `{row['target_signature']}` closes `{row['newly_closed_targets']}`")
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
    print(json.dumps(payload["theorem_frontier_low_weight_extension_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
