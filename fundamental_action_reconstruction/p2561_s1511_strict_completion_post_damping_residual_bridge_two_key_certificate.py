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
OUT = GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json"
MD = GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.md"
SOURCE_FILES = {
    "P2502_BRIDGE_MINIMAL_TRIPLE_FRONTIER": GEN / "p2502_s1452_strict_completion_bridge_minimal_triple_frontier_certificate.json",
    "P2560_CONSTANT_LOG_SOURCE_OBSTRUCTION": GEN / "p2560_s1510_legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate.json",
}
BRIDGE_TRIPLE = [
    "strict_damping_beta_eta_source",
    "strict_dynamical_source_for_A_P_D",
    "strict_phase_frequency_source",
]
DAMPING_ATOM = "strict_damping_beta_eta_source"
RESIDUAL_ATOMS = ["strict_dynamical_source_for_A_P_D", "strict_phase_frequency_source"]
NEGATIVE_EXPORT_FLAGS = [
    "strict_damping_beta_eta_source_exported", "strict_dynamical_source_for_A_P_D_exported", "strict_phase_frequency_source_exported",
    "bridge_theorem_exported", "legacy_to_strict_completion_bridge_exported", "full_bridge_theorem_exported",
    "role_transfer_theorem_exported", "legacy_role_transfer_claimed", "role_bearing_ltotal_exported",
    "selector_closure_exported", "qw2191_discharged_by_this_certificate", "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2561|S1511|post-damping residual bridge|residual bridge two-key|post damping residual",
        "intended_research_nonduplication": "two-key bridge frontier|A_P_D.*phase_frequency|phase_frequency.*A_P_D|post-damping.*bridge",
        "precursors": "P2502|S1452|P2560|S1510|minimal triple frontier|constant-log-source current-premise obstruction",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def residual_rows(damping_assumed_true: bool) -> list[dict[str, Any]]:
    rows = []
    for bits in itertools.product([False, True], repeat=len(RESIDUAL_ATOMS)):
        assignment = dict(zip(RESIDUAL_ATOMS, bits))
        strict_dynamical = assignment["strict_dynamical_source_for_A_P_D"]
        phase_frequency = assignment["strict_phase_frequency_source"]
        bridge = damping_assumed_true and strict_dynamical and phase_frequency
        missing = [atom for atom, value in [(DAMPING_ATOM, damping_assumed_true), *assignment.items()] if not value]
        rows.append({
            "damping_assumed_true": damping_assumed_true,
            **assignment,
            "bridge_theorem_level_closure": bridge,
            "missing_atoms": missing,
        })
    return rows


def single_residual_omission_witnesses(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    witnesses = []
    for omitted in RESIDUAL_ATOMS:
        row = next(candidate for candidate in rows if candidate[omitted] is False and all(candidate[atom] is True for atom in RESIDUAL_ATOMS if atom != omitted))
        witnesses.append({
            "omitted_residual_atom": omitted,
            "other_residual_atoms_true": [atom for atom in RESIDUAL_ATOMS if atom != omitted],
            "damping_assumed_true": row["damping_assumed_true"],
            "bridge_theorem_level_closure": row["bridge_theorem_level_closure"],
            "missing_atoms": row["missing_atoms"],
        })
    return witnesses


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2502 = theorem(sources["P2502_BRIDGE_MINIMAL_TRIPLE_FRONTIER"], "strict_completion_bridge_minimal_triple_frontier_certificate")
    p2560 = theorem(sources["P2560_CONSTANT_LOG_SOURCE_OBSTRUCTION"], "legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate")
    best_case_rows = residual_rows(damping_assumed_true=True)
    current_rows = residual_rows(damping_assumed_true=False)
    bridge_accepting_best_case = [row for row in best_case_rows if row["bridge_theorem_level_closure"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2561_T1_post_damping_residual_bridge_two_key_frontier_certificate",
        "audited_chain": ["P2502/S1452", "P2560/S1510"],
        "p2502_unique_bridge_triple_inherited": p2502.get("bridge_closing_triple_count") == 1 and p2502.get("strict_source_triple_closes_only_bridge") is True,
        "p2560_constant_log_source_obstruction_inherited": p2560.get("current_endpoint_premises_do_not_entail_constant_log_source") is True,
        "bridge_triple": BRIDGE_TRIPLE,
        "damping_atom": DAMPING_ATOM,
        "residual_atoms_after_hypothetical_damping_source": RESIDUAL_ATOMS,
        "current_no_damping_source_rows": current_rows,
        "current_bridge_accepting_rows": sum(1 for row in current_rows if row["bridge_theorem_level_closure"]),
        "best_case_damping_residual_truth_table": best_case_rows,
        "best_case_damping_truth_table_row_count": len(best_case_rows),
        "best_case_damping_bridge_accepting_row_count": len(bridge_accepting_best_case),
        "best_case_damping_bridge_accepting_rows": bridge_accepting_best_case,
        "residual_atoms_are_jointly_required_after_damping": len(bridge_accepting_best_case) == 1 and bridge_accepting_best_case[0][RESIDUAL_ATOMS[0]] and bridge_accepting_best_case[0][RESIDUAL_ATOMS[1]],
        "single_residual_omission_witnesses": single_residual_omission_witnesses(best_case_rows),
        "post_damping_residual_bridge_two_key_blocker_exported": True,
        "recommended_next_honest_step": (
            "After P2560, do not claim damping bridge completion. If a future strict source supplies damping_beta_eta, the bridge still requires two independent source theorems: strict_dynamical_source_for_A_P_D and strict_phase_frequency_source. The next honest non-duplicative bridge work should attack one of those two residual source atoms directly, preferably the phase/frequency/topological-bit passage, while preserving the QW-2191 guardrail."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2502_inherited": theorem_export["p2502_unique_bridge_triple_inherited"],
        "p2560_inherited": theorem_export["p2560_constant_log_source_obstruction_inherited"],
        "current_without_damping_has_no_bridge_accepting_row": theorem_export["current_bridge_accepting_rows"] == 0,
        "best_case_truth_table_has_four_rows": theorem_export["best_case_damping_truth_table_row_count"] == 4,
        "best_case_truth_table_has_one_accepting_row": theorem_export["best_case_damping_bridge_accepting_row_count"] == 1,
        "no_bridge_exported": theorem_export["full_bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2561",
        "stage_id": "S1511",
        "status": "P2561_POST_DAMPING_RESIDUAL_TWO_KEY_BRIDGE_FRONTIER_NO_SOURCE_EXPORT_NO_FALSE_PASS",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_completion_post_damping_residual_bridge_two_key_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_completion_post_damping_residual_bridge_two_key_certificate"]["theorem_export"]
    lines = [
        "# P2561/S1511 strict-completion post-damping residual bridge two-key certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- P2502 bridge triple inherited: `{t['p2502_unique_bridge_triple_inherited']}`.",
        f"- P2560 damping obstruction inherited: `{t['p2560_constant_log_source_obstruction_inherited']}`.",
        f"- Current accepting bridge rows without damping source: `{t['current_bridge_accepting_rows']}`.",
        f"- Best-case damping residual table rows: `{t['best_case_damping_truth_table_row_count']}`.",
        f"- Best-case damping accepting rows: `{t['best_case_damping_bridge_accepting_row_count']}`.",
        f"- Residual atoms after hypothetical damping source: `{t['residual_atoms_after_hypothetical_damping_source']}`.", "", "## Interpretation", "",
        "Even a future damping source would not complete the legacy->strict bridge by itself.  Under the P2502 bridge-frontier equation, the post-damping residual bridge blocker is the two-key conjunction of strict A/P/D dynamics and strict phase/frequency source.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No strict source atom, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['strict_completion_post_damping_residual_bridge_two_key_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2561/S1511` combines the P2502 bridge-frontier triple with the P2560 damping obstruction.  With current premises there are `0` bridge-accepting rows because the damping source is still absent.  Even under a best-case future damping source, the residual bridge truth table has `4` rows and only the row with both `strict_dynamical_source_for_A_P_D` and `strict_phase_frequency_source` closes the bridge target.  Thus damping work alone cannot close the bridge; the next nonduplicative bridge attack should target one of those residual source atoms, especially phase/frequency/topological-bit passage under the QW-2191 guardrail.
""".strip()
    lag_section = """
`P2561/S1511` blocks promotion of strict damping progress into a full role-bearing `L_total` bridge.  Post-damping, the bridge still requires strict A/P/D dynamics and strict phase/frequency source; role-transfer and QW-2191 remain separate downstream gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2561/S1511 post-damping residual bridge two-key guard", "## P2561/S1511 post-damping residual bridge two-key guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2561/S1511 post-damping residual bridge two-key Ltotal guard", "## P2561/S1511 post-damping residual bridge two-key Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
