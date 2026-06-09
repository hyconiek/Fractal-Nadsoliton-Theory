#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2620_s1570_p2618_p2619_bridge_two_obstruction_cut.json"
MD = GEN / "p2620_s1570_p2618_p2619_bridge_two_obstruction_cut.md"

SOURCE_FILES = {
    "P2604_BRIDGE_READINESS": GEN / "p2604_s1554_strict_damping_post_source_bridge_readiness_matrix.json",
    "P2616_ROLE_OBSTRUCTION": GEN / "p2616_s1566_p2608_role_acceptance_obstruction_after_source_revalidation.json",
    "P2618_ANALYTIC_OBSTRUCTION": GEN / "p2618_s1568_analytic_legacy_to_strict_completion_obstruction.json",
    "P2619_SELECTOR_LATTICE": GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.json",
}

BRIDGE_REPAIR_ATOMS = [
    "nonlinear_damping_completion_source",
    "orientation_odd_selector_source",
]

NEGATIVE_EXPORT_FLAGS = [
    "full_bridge_revalidated",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "beta_tors_chi11_route_reopened",
    "gf2_bridge_revalidated",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2620|S1570|two-obstruction cut|bridge two obstruction|P2618 P2619 bridge cut",
        "nondup_bridge_matrices": "P2604|bridge readiness matrix|completion gate|strict-side residual additions|role-transfer theorem",
        "p2618_p2619_precursors": "P2618|P2619|nonlinear_damping_completion_source|orientation_odd_selector_source|beta_tors.*chi11",
        "guardrails": "K_legacy_ont|K_strict_gate|role-bearing L_total|QW-2191|ToE closure|legacy physical-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def bridge_cut_rows() -> list[dict[str, Any]]:
    rows = []
    for values in itertools.product([False, True], repeat=len(BRIDGE_REPAIR_ATOMS)):
        assignment = dict(zip(BRIDGE_REPAIR_ATOMS, values))
        bridge_repaired = all(values)
        missing = [atom for atom, value in assignment.items() if not value]
        rows.append({
            "assignment": assignment,
            "bridge_source_cut_repaired": bridge_repaired,
            "missing_repair_atoms": missing,
            "missing_repair_atom_count": len(missing),
        })
    return rows


def shortcut_rejection_rows() -> list[dict[str, Any]]:
    shortcuts = [
        {
            "shortcut": "eta_9_5_exponent_source_only",
            "has_nonlinear_damping_completion_source": False,
            "has_orientation_odd_selector_source": False,
            "reason": "P2618 retains the exponent source but blocks exact scalar damping completion and selector source.",
        },
        {
            "shortcut": "beta_tors_scalar_renormalization",
            "has_nonlinear_damping_completion_source": False,
            "has_orientation_odd_selector_source": False,
            "reason": "P2618 derivative obstruction blocks scalar damping completion; P2619 treats beta_tors as orientation-invariant.",
        },
        {
            "shortcut": "axis_only_selector_up_to_Z2",
            "has_nonlinear_damping_completion_source": False,
            "has_orientation_odd_selector_source": False,
            "reason": "Axis-only data can reduce continuous degeneracy but leaves the strict odd sign unresolved.",
        },
        {
            "shortcut": "GF2_bookkeeping_rank_or_cycle_data",
            "has_nonlinear_damping_completion_source": False,
            "has_orientation_odd_selector_source": False,
            "reason": "Combinatorial bookkeeping may classify constraints but is not a physical orientation/source premise under P2612/P2619.",
        },
        {
            "shortcut": "damping_source_plus_no_selector",
            "has_nonlinear_damping_completion_source": True,
            "has_orientation_odd_selector_source": False,
            "reason": "Even a future damping completion source alone leaves the P2619 C2 selector obstruction.",
        },
        {
            "shortcut": "selector_source_plus_no_damping_completion",
            "has_nonlinear_damping_completion_source": False,
            "has_orientation_odd_selector_source": True,
            "reason": "Even a future selector source alone leaves the P2618 nonlinear damping completion obstruction.",
        },
        {
            "shortcut": "both_independent_sources_supplied",
            "has_nonlinear_damping_completion_source": True,
            "has_orientation_odd_selector_source": True,
            "reason": "This is the minimal bridge-source cut repair, still not a role-transfer theorem.",
        },
    ]
    for row in shortcuts:
        row["passes_bridge_source_cut"] = row["has_nonlinear_damping_completion_source"] and row["has_orientation_odd_selector_source"]
    return shortcuts


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    payloads = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2604 = theorem(payloads["P2604_BRIDGE_READINESS"], "strict_damping_post_source_bridge_readiness_matrix")
    p2616 = theorem(payloads["P2616_ROLE_OBSTRUCTION"], "p2608_role_acceptance_obstruction_after_source_revalidation")
    p2618 = theorem(payloads["P2618_ANALYTIC_OBSTRUCTION"], "analytic_legacy_to_strict_completion_obstruction")
    p2619 = theorem(payloads["P2619_SELECTOR_LATTICE"], "p2618_selector_source_obligation_lattice")
    rows = bridge_cut_rows()
    shortcuts = shortcut_rejection_rows()
    accepting = [row for row in rows if row["bridge_source_cut_repaired"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2620_T1_p2618_p2619_bridge_two_obstruction_cut",
        "result_type": "exact_two_atom_bridge_source_cut_after_P2618_P2619",
        "inherits_p2604_bridge_role_matrix": p2604.get("current_role_bearing_ltotal_accepts") is False,
        "inherits_p2616_role_obstruction": p2616.get("current_assignment_accepts_role_bearing_ltotal") is False or p2616.get("p2608_role_bearing_ltotal_reenabled") is False,
        "inherits_p2618_full_bridge_block": p2618.get("full_completion_map_exported") is False,
        "inherits_p2619_selector_source_block": p2619.get("strict_selector_source_exported") is False,
        "formal_cut_theorem": {
            "statement": "After P2618/P2619, any non-role-bearing legacy-to-strict bridge-source repair must supply both an independent nonlinear damping completion source and an independent orientation-odd selector source. Either one alone leaves a named obstruction alive.",
            "proof_steps": [
                "P2618 proves that the strict exponent source eta=9/5 does not provide an exact scalar beta_tors -> beta damping completion map.",
                "P2619 proves that orientation-invariant legacy scalar or axis-only data cannot provide a C2-equivariant strict odd phase-sign selector.",
                "The damping obstruction and selector obstruction live in different output coordinates of the completion map: denominator compression versus phase/topological sign.",
                "Therefore repairing only the damping coordinate leaves the selector obstruction, and repairing only the selector coordinate leaves the damping obstruction.",
                "The bridge-source cut is the conjunction of the two independent repair atoms; the finite truth table has exactly one accepting row.",
            ],
        },
        "bridge_repair_atoms": BRIDGE_REPAIR_ATOMS,
        "bridge_cut_truth_table_rows": rows,
        "bridge_cut_truth_table_row_count": len(rows),
        "bridge_cut_accepting_row_count": len(accepting),
        "bridge_cut_accepting_row": accepting[0],
        "current_assignment": {atom: False for atom in BRIDGE_REPAIR_ATOMS},
        "current_bridge_source_cut_repaired": False,
        "shortcut_rejection_rows": shortcuts,
        "accepted_shortcut_count": sum(1 for row in shortcuts if row["passes_bridge_source_cut"]),
        "recommended_next_honest_step": "Pick one real source target, not a shortcut: either derive a nonlinear damping completion source beyond scalar beta_tors renormalization, or derive an orientation-odd selector source; full bridge-source repair needs both before role transfer can be rerun.",
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2604_matrix_inherited": theorem_export["inherits_p2604_bridge_role_matrix"],
        "p2616_role_block_inherited": theorem_export["inherits_p2616_role_obstruction"],
        "p2618_bridge_block_inherited": theorem_export["inherits_p2618_full_bridge_block"],
        "p2619_selector_block_inherited": theorem_export["inherits_p2619_selector_source_block"],
        "truth_table_has_four_rows": theorem_export["bridge_cut_truth_table_row_count"] == 4,
        "truth_table_has_one_accepting_row": theorem_export["bridge_cut_accepting_row_count"] == 1,
        "current_assignment_rejected": theorem_export["current_bridge_source_cut_repaired"] is False,
        "one_atom_shortcuts_rejected": all(not row["passes_bridge_source_cut"] for row in shortcuts if row["shortcut"] in {"damping_source_plus_no_selector", "selector_source_plus_no_damping_completion"}),
        "both_sources_accepts_cut_only": [row for row in shortcuts if row["passes_bridge_source_cut"]] == [shortcuts[-1]],
        "no_full_bridge_revalidated": theorem_export["full_bridge_revalidated"] is False,
        "no_role_transfer_revalidated": theorem_export["role_transfer_revalidated"] is False,
        "no_ltotal_reenabled": theorem_export["role_bearing_ltotal_reenabled"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2620",
        "stage_id": "S1570",
        "status": "P2620_TWO_OBSTRUCTION_BRIDGE_SOURCE_CUT_EXACT_TRUTH_TABLE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "p2618_p2619_bridge_two_obstruction_cut": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in payloads.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2618_p2619_bridge_two_obstruction_cut"]["theorem_export"]
    lines = [
        "# P2620/S1570 P2618/P2619 bridge two-obstruction cut", "",
        f"Status: `{payload['status']}`", "", "## Theorem", "",
        t["formal_cut_theorem"]["statement"], "", "## Proof", "",
    ]
    for step in t["formal_cut_theorem"]["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend(["", "## Computed truth table", ""])
    for row in t["bridge_cut_truth_table_rows"]:
        lines.append(f"- `{row['assignment']}` -> bridge-source cut repaired `{row['bridge_source_cut_repaired']}`, missing `{row['missing_repair_atoms']}`.")
    lines.extend(["", "## Shortcut rejection matrix", ""])
    for row in t["shortcut_rejection_rows"]:
        lines.append(f"- `{row['shortcut']}`: passes `{row['passes_bridge_source_cut']}` — {row['reason']}")
    lines.extend([
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"],
        "", "## Scope guards", "",
        "No full bridge revalidation, no role-transfer revalidation, no role-bearing `L_total`, no `beta_tors -> chi11` reopening, no GF(2) bridge revalidation, no `QW-2191` discharge, and no ToE closure are exported.",
        "", "## Fingerprint", "",
        f"`{payload['p2618_p2619_bridge_two_obstruction_cut']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2620/S1570 P2618/P2619 bridge two-obstruction cut

`P2620/S1570` converts P2618 and P2619 into an exact two-atom bridge-source cut.  A non-role-bearing bridge repair must supply both an independent nonlinear damping completion source and an independent orientation-odd selector source; either one alone leaves the other named obstruction alive.  The computed truth table has four rows and exactly one accepting bridge-source-cut row, and the accepted row is still not a role-transfer theorem.
""".strip()
    lag_section = """
## P2620/S1570 two-obstruction Ltotal guard

`P2620/S1570` keeps role-bearing `L_total` closed after the P2618/P2619 refinements.  Neither the `eta=9/5` exponent source, `beta_tors` scalar renormalization, axis-only/Z2 selector data, nor GF(2) bookkeeping repairs the bridge-source cut; both nonlinear damping completion and orientation-odd selector sources are required before bridge validity and role transfer may be rerun.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2620/S1570 P2618/P2619 bridge two-obstruction cut", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2620/S1570 two-obstruction Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
