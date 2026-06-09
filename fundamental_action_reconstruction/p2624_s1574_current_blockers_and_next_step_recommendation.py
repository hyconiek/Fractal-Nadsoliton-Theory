#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2624_s1574_current_blockers_and_next_step_recommendation.json"
MD = GEN / "p2624_s1574_current_blockers_and_next_step_recommendation.md"

SOURCE_FILES = {
    "K1_KERNEL_SPLIT_NOTE": ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "P2618_COMPLETION_OBSTRUCTION": GEN / "p2618_s1568_analytic_legacy_to_strict_completion_obstruction.json",
    "P2619_SELECTOR_LATTICE": GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.json",
    "P2620_TWO_OBSTRUCTION_CUT": GEN / "p2620_s1570_p2618_p2619_bridge_two_obstruction_cut.json",
    "P2621_CHIRAL_SCHEMA": GEN / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.json",
    "P2622_QW_NONPROMOTION": GEN / "p2622_s1572_qw636_qw1026_physical_rigor_nonpromotion_audit.json",
    "P2623_WILSON_BOUNDARY": GEN / "p2623_s1573_wilson_loop_flux_orientation_source_boundary.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "strict_kernel_promoted_to_full_kernel",
    "toe_closure_claimed",
    "qw2191_discharged",
    "orientation_odd_selector_source_exported_unconditionally",
    "nonlinear_damping_completion_source_exported",
    "full_legacy_to_strict_bridge_revalidated",
    "legacy_role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:220]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first patterns: these intentionally search research concepts, not just P/QW labels.
    patterns = {
        "selector_orientation_chirality_content": (
            "orientation-odd|orientation source|selector source|phase sign|sign selector|"
            "chiral anomaly|gamma5|Hopfion|parity-odd|Wilson loop|holonomy|oriented cycle|cycle orientation"
        ),
        "nonlinear_damping_completion_content": (
            "nonlinear damping|nonlinear compression|d\\^eta|fractal projection|compressible|hydrodynamic|"
            "beta_tors|strict damping|attenuation|denominator|completion map"
        ),
        "bridge_role_transfer_content": (
            "legacy-to-strict|legacy -> strict|completion bridge|completion-map|role-transfer|"
            "role-bearing L_total|L_total|strict bridge|physical-role"
        ),
        "toe_full_kernel_content": (
            "ToE closure|full kernel|complete kernel|kernel finality|kernel identity|strict-core closure|"
            "QW-2191|selector closure|global closure"
        ),
        "methodology_guard_content": (
            "no numerical fit|no fit|a priori|gauge invariant|basis-independent|symmetry-breaking|"
            "nonpromotion|obstruction|quarantine|guard"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic anti-duplication audit; ticket numbers are not the search key", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def truth_table(atoms: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for values in itertools.product([False, True], repeat=len(atoms)):
        assignment = dict(zip(atoms, values))
        bridge_ready = assignment["nonlinear_damping_completion_source"] and assignment["orientation_odd_selector_source"]
        ltotal_ready = bridge_ready and assignment["role_transfer_theorem"]
        toe_ready = ltotal_ready and assignment["qw2191_selector_discharge"] and assignment["global_kernel_finality_theorem"]
        rows.append({
            "assignment": assignment,
            "bridge_source_cut_repaired": bridge_ready,
            "role_bearing_ltotal_reenabled": ltotal_ready,
            "strict_kernel_is_full_toe_kernel": toe_ready,
            "missing_for_toe": [name for name, value in assignment.items() if not value],
        })
    return rows


def blocker_vector() -> dict[str, Any]:
    atoms = [
        "nonlinear_damping_completion_source",
        "orientation_odd_selector_source",
        "role_transfer_theorem",
        "qw2191_selector_discharge",
        "global_kernel_finality_theorem",
    ]
    current = {
        "nonlinear_damping_completion_source": False,
        "orientation_odd_selector_source": False,
        "role_transfer_theorem": False,
        "qw2191_selector_discharge": False,
        "global_kernel_finality_theorem": False,
    }
    rows = truth_table(atoms)
    current_row = next(row for row in rows if row["assignment"] == current)
    return {
        "atoms": atoms,
        "current_assignment": current,
        "truth_table_row_count": len(rows),
        "bridge_repair_accepting_rows": sum(1 for row in rows if row["bridge_source_cut_repaired"]),
        "ltotal_accepting_rows": sum(1 for row in rows if row["role_bearing_ltotal_reenabled"]),
        "toe_accepting_rows": sum(1 for row in rows if row["strict_kernel_is_full_toe_kernel"]),
        "current_row": current_row,
    }


def symptom_matrix() -> list[dict[str, Any]]:
    return [
        {
            "symptom": "strict RG normalization and prime-log character are repaired",
            "status": "positive_but_local",
            "why_not_full_kernel": "They constrain scale transport/attenuation normalization, not phase-source uniqueness, legacy bridge completion, or global finality.",
        },
        {
            "symptom": "strict exponent eta=9/5 has a fractal-codimension reading",
            "status": "positive_but_partial",
            "why_not_full_kernel": "It supplies exponent semantics, but not an exact beta_tors -> beta completion map and not a nonlinear damping source theorem strong enough to repair P2620.",
        },
        {
            "symptom": "chiral/Hopfion/Wilson mechanisms show ToE-like parity and topology diagnostics",
            "status": "diagnostic_only",
            "why_not_full_kernel": "Their signs remain conditional on typed orientation, gauge-safe cycle, gamma5, or nonzero-flux premises; prior art alone does not export the strict selector source.",
        },
        {
            "symptom": "K_strict_gate is a powerful completed/enriched working kernel candidate",
            "status": "not_a_full_kernel",
            "why_not_full_kernel": "The current guardrails require explicit bridge completion, role-transfer audit, QW-2191 discharge, and global finality before calling it a full ToE kernel.",
        },
    ]


def recommendation() -> dict[str, Any]:
    return {
        "recommended_next_packet": "P2625 nonlinear damping completion-source classification",
        "recommended_focus": "stop cycling through selector/chirality prior art for one step; attack the independent nonlinear_damping_completion_source atom in P2620",
        "proof_goal": (
            "Classify whether a local/scale-covariant fractal hydrodynamic attenuation law can derive the strict denominator "
            "1+beta*d^(9/5) from legacy linear torsion data plus an explicitly typed compressibility/source field, or prove that an independent beta source is unavoidable."
        ),
        "computational_goal": (
            "Enumerate exact finite countermodels/positive models for candidate completion laws: scalar beta_tors renormalization, scale-dependent beta_eff(d), "
            "fractional-measure pushforward, and compressibility-field source.  Accept only models that are a priori, source-typed, non-fit, and not reducible to scalar rescaling."
        ),
        "why_this_is_more_honest_than_more_selector_search": (
            "P2621-P2623 already audited the selector/chirality/Wilson lane and found only conditional diagnostics; P2620 still has a second independent open atom, so nonlinear damping is now the least-duplicative high-value target."
        ),
        "forbidden_shortcuts": [
            "do not call eta=9/5 alone a full completion map",
            "do not fit beta numerically from beta_tors",
            "do not use QW-636/QW-1026/Wilson diagnostics as an unconditional selector source",
            "do not reopen role-bearing L_total until bridge and role-transfer predicates are rerun",
            "do not call K_strict_gate a full ToE kernel from symptoms alone",
        ],
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    blockers = blocker_vector()
    payload = {
        "packet_id": "P2624",
        "slice_id": "S1574",
        "status": "P2624_CURRENT_BLOCKERS_AND_NEXT_STEP_RECOMMENDATION_NO_FULL_KERNEL_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "inherited_artifact_status": {name: art.get("status", art.get("packet_id", "UNKNOWN")) for name, art in artifacts.items()},
        "current_blocker_vector": blockers,
        "strict_kernel_toe_symptom_matrix": symptom_matrix(),
        "answer_to_full_kernel_question": {
            "short_answer": "No: ToE-like symptoms of K_strict_gate do not make it a full kernel under current repo guardrails.",
            "necessary_condition_failed": blockers["current_row"]["missing_for_toe"],
            "precise_status": "K_strict_gate is a strict working/completed-enriched candidate kernel, not an admitted full ToE kernel.",
        },
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["certificate_hash"] = sha256_json({k: v for k, v in payload.items() if k != "certificate_hash"})
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    blockers = payload["current_blocker_vector"]
    rec = payload["recommended_next_honest_step"]
    answer = payload["answer_to_full_kernel_question"]
    lines = [
        "# P2624/S1574 Current blockers and next-step recommendation",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication grep audit",
        "",
        f"Mode: `{payload['semantic_rg_antiduplication_audit']['mode']}`.",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits; samples retained in JSON certificate.")
    lines.extend([
        "",
        "## Current closure blockers",
        "",
        f"Current assignment: `{blockers['current_assignment']}`.",
        f"Truth table rows: `{blockers['truth_table_row_count']}`.",
        f"Bridge-repair accepting rows: `{blockers['bridge_repair_accepting_rows']}`.",
        f"Role-bearing L_total accepting rows: `{blockers['ltotal_accepting_rows']}`.",
        f"Full-ToE-kernel accepting rows: `{blockers['toe_accepting_rows']}`.",
        f"Current missing atoms for ToE/kernel finality: `{blockers['current_row']['missing_for_toe']}`.",
        "",
        "## Is the strict kernel a full kernel if it shows ToE-like symptoms?",
        "",
        f"**Short answer:** {answer['short_answer']}",
        "",
        f"Precise status: `{answer['precise_status']}`.",
        "",
        "Symptom matrix:",
    ])
    for row in payload["strict_kernel_toe_symptom_matrix"]:
        lines.append(f"- `{row['status']}`: {row['symptom']} — {row['why_not_full_kernel']}")
    lines.extend([
        "",
        "## Recommended next honest step",
        "",
        f"Recommended packet: `{rec['recommended_next_packet']}`.",
        f"Focus: {rec['recommended_focus']}.",
        "",
        f"Proof goal: {rec['proof_goal']}",
        "",
        f"Computational goal: {rec['computational_goal']}",
        "",
        f"Why this is least-duplicative: {rec['why_this_is_more_honest_than_more_selector_search']}",
        "",
        "Forbidden shortcuts:",
    ])
    for item in rec["forbidden_shortcuts"]:
        lines.append(f"- {item}")
    lines.extend([
        "",
        "## Negative export flags",
        "",
    ])
    for flag, value in payload["negative_export_flags"].items():
        lines.append(f"- `{flag}`: `{value}`")
    lines.append("")
    MD.write_text("\n".join(lines), encoding="utf-8")


def update_docs() -> None:
    equation_note = """
## P2624/S1574 current blocker and next-step recommendation guard

`P2624/S1574` answers the current-status question without promoting symptoms into closure.  The strict gate kernel may exhibit ToE-like local symptoms (RG normalization, `eta=9/5` exponent semantics, and conditional topological/parity diagnostics), but it is not a full ToE kernel while `nonlinear_damping_completion_source`, unconditional `orientation_odd_selector_source`, role-transfer, `QW-2191` discharge, and global finality remain false.  The least-duplicative next proof target is therefore a content-first nonlinear damping completion-source classification, not another selector/chirality loop.
"""
    ltotal_note = """
## P2624/S1574 current blocker Ltotal guard

`P2624/S1574` keeps role-bearing `L_total` closed.  It records that ToE-like symptoms of `K_strict_gate` are not enough to promote it to a full kernel: bridge-source repair, a separate role-transfer theorem, `QW-2191` discharge, and global kernel-finality theorem are still missing.  The recommended next admissible work item is a non-fit analytic/computational classification of the independent nonlinear damping completion-source atom.
"""
    append_once(ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md", "P2624/S1574 current blocker and next-step recommendation guard", equation_note)
    append_once(ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md", "P2624/S1574 current blocker Ltotal guard", ltotal_note)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    update_docs()
    print(json.dumps({
        "packet_id": payload["packet_id"],
        "status": payload["status"],
        "current_missing_for_toe": payload["current_blocker_vector"]["current_row"]["missing_for_toe"],
        "recommended_next_packet": payload["recommended_next_honest_step"]["recommended_next_packet"],
        "strict_kernel_full_kernel": payload["answer_to_full_kernel_question"]["short_answer"],
        "certificate_hash": payload["certificate_hash"],
    }, indent=2, sort_keys=True, ensure_ascii=False))


if __name__ == "__main__":
    main()
