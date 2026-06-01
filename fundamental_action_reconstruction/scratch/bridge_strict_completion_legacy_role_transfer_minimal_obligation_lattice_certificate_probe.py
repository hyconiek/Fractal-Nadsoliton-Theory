#!/usr/bin/env python3
"""Scratch probe: minimal obligation lattice for legacy role transfer.

The pre-audit blocks every legacy physical role.  This probe turns that blocked
status into a finite obligation lattice: which missing theorem atoms would be
needed to license each role, and which minimal atom set would be needed to even
attempt all listed legacy-role transfers.

It is not a role-transfer theorem.  It is a small exhaustive set-cover
calculation over the currently audited legacy role rows.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.md"

SOURCE_REPORTS = {
    "role_transfer_pre_audit": HERE / "bridge_strict_completion_legacy_role_transfer_pre_audit_certificate_report.json",
    "symbolic_cancellation": HERE / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
}

THEOREM_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
    "chi11_selector_source_theorem",
]

ROLE_REQUIREMENTS = {
    "legacy_weak_mixing_angle": ["alpha_geo_electroweak_role_theorem"],
    "legacy_inverse_alpha_em": ["alpha_geo_electroweak_role_theorem", "beta_tors_strict_role_theorem"],
    "legacy_beta_power_gravity_hierarchy": ["beta_tors_strict_role_theorem", "beta_power_hierarchy_successor_theorem"],
    "legacy_torsion_to_chi11_orientation": ["beta_tors_strict_role_theorem", "chi11_selector_source_theorem"],
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def powerset(items: list[str]) -> list[set[str]]:
    subsets: list[set[str]] = []
    for size in range(len(items) + 1):
        for combo in itertools.combinations(items, size):
            subsets.append(set(combo))
    return subsets


def minimal_sets_for(required: set[str], universe: list[str]) -> list[list[str]]:
    candidates = [subset for subset in powerset(universe) if required <= subset]
    min_size = min(len(subset) for subset in candidates)
    return [[atom for atom in universe if atom in subset] for subset in candidates if len(subset) == min_size]


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    pre_audit = loaded["role_transfer_pre_audit"]
    symbolic = loaded["symbolic_cancellation"]
    component_gap = loaded["component_gap_matrix"]
    guardrail = loaded["legacy_bridge_guardrail"]

    audited_role_ids = [row["role_id"] for row in pre_audit["role_transfer_rows"]]
    all_subsets = powerset(THEOREM_ATOMS)
    role_lattice_rows = []
    for role_id in audited_role_ids:
        required = set(ROLE_REQUIREMENTS[role_id])
        satisfying_subsets = [sorted(subset) for subset in all_subsets if required <= subset]
        minimal = minimal_sets_for(required, THEOREM_ATOMS)
        role_lattice_rows.append(
            {
                "role_id": role_id,
                "required_theorem_atoms": sorted(required),
                "requirement_count": len(required),
                "minimal_satisfying_sets": minimal,
                "satisfying_subset_count": len(satisfying_subsets),
                "blocked_now_because_atoms_exported": [],
                "transfer_allowed_now": False,
            }
        )

    global_required = sorted(set().union(*(set(values) for values in ROLE_REQUIREMENTS.values())))
    global_minimal_sets = minimal_sets_for(set(global_required), THEOREM_ATOMS)
    atom_coverage_rows = []
    for atom in THEOREM_ATOMS:
        covered_roles = sorted(role_id for role_id, reqs in ROLE_REQUIREMENTS.items() if atom in reqs)
        atom_coverage_rows.append(
            {
                "theorem_atom": atom,
                "covered_roles": covered_roles,
                "covered_role_count": len(covered_roles),
                "exported_now": False,
            }
        )

    summary = {
        "role_count": len(audited_role_ids),
        "theorem_atom_count": len(THEOREM_ATOMS),
        "total_subset_count_checked": len(all_subsets),
        "all_pre_audit_roles_loaded": set(audited_role_ids) == set(ROLE_REQUIREMENTS),
        "all_roles_blocked_in_pre_audit": pre_audit["role_transfer_pre_audit_summary"]["all_roles_blocked_now"],
        "all_atoms_missing_now": all(not row["exported_now"] for row in atom_coverage_rows),
        "global_minimal_obligation_count": len(global_minimal_sets),
        "global_minimal_obligation_size": len(global_minimal_sets[0]),
        "global_minimal_obligation_sets": [THEOREM_ATOMS] if set(global_required) == set(THEOREM_ATOMS) else global_minimal_sets,
        "beta_tors_atom_is_shared_by_three_roles": next(row for row in atom_coverage_rows if row["theorem_atom"] == "beta_tors_strict_role_theorem")["covered_role_count"] == 3,
        "symbolic_bridge_still_no_role_transfer": not symbolic["symbolic_cancellation_summary"]["legacy_role_transfer_exported"],
        "component_gap_still_blocks_role_transfer": component_gap["completion_gap_summary"]["role_transfer_blocked_until_full_bridge"],
        "guardrail_requires_audit": guardrail["legacy_kernel_intermediate_bridge_summary"]["role_transfer_audit_required_after_full_bridge"],
        "role_transfer_theorem_exported": False,
        "q2191_discharged": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_ROLE_TRANSFER_MINIMAL_OBLIGATION_LATTICE_CERTIFICATE__NO_TRANSFER_THEOREM",
        "status": "minimal-role-transfer-obligation-lattice-enumerated-all-theorem-atoms-missing",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "role obligation",
                "minimal theorem",
                "role-transfer theorem obligation",
                "alpha_geo role",
                "beta_tors role",
                "chi11 role",
            ],
            "finding": "The repo already had a pre-audit blocking all roles; this report enumerates the minimal missing theorem atoms needed to unblock them without claiming any atom is present.",
        },
        "theorem_atoms": THEOREM_ATOMS,
        "role_requirement_lattice_rows": role_lattice_rows,
        "atom_coverage_rows": atom_coverage_rows,
        "role_transfer_minimal_obligation_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "pre_audit_roles_and_blocks_inherited": summary["all_pre_audit_roles_loaded"] and summary["all_roles_blocked_in_pre_audit"],
            "global_minimal_set_is_all_atoms": summary["global_minimal_obligation_size"] == len(THEOREM_ATOMS) and summary["global_minimal_obligation_sets"] == [THEOREM_ATOMS],
            "beta_tors_shared_obligation_detected": summary["beta_tors_atom_is_shared_by_three_roles"],
            "all_atoms_missing_and_limits_preserved": summary["all_atoms_missing_now"] and not summary["role_transfer_theorem_exported"] and not summary["q2191_discharged"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to distinguish this minimal-obligation lattice from the earlier pre-audit and obstruction notes.",
            "enumeration_step": f"All {len(all_subsets)} subsets of the four theorem atoms were enumerated.",
            "role_step": "Each role is satisfied only by subsets containing its required theorem atoms; with zero atoms exported now, every transfer remains blocked.",
            "global_step": "To cover all four audited roles, the unique minimal obligation set is the full four-atom set: alpha_geo role, beta_tors role, beta-power successor, and chi11 selector/source.",
            "shared_beta_step": "The beta_tors strict-role theorem is shared by alpha_EM, beta^N hierarchy, and beta_tors->chi11 obligations, so it is the widest missing role-transfer bottleneck.",
            "scope_step": "This lattice is a planning certificate only; it exports no role-transfer theorem, no QW-2191 discharge, and no ToE closure.",
        },
        "hard_limits": [
            "No theorem atom in the role-transfer obligation lattice is exported here.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No alpha_geo, beta_tors, beta^N, or chi11 role theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy role-transfer minimal obligation lattice certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "All subsets of the missing role-transfer theorem atoms are enumerated.",
        "The unique minimal set covering all audited legacy roles contains all four",
        "atoms, and none is exported by this certificate.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["role_transfer_minimal_obligation_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Role lattice rows", ""])
    for row in payload["role_requirement_lattice_rows"]:
        lines.append(f"- `{row['role_id']}` requires `{row['required_theorem_atoms']}`")
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
    print(json.dumps(payload["role_transfer_minimal_obligation_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
