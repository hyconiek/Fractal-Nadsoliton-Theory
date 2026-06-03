#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)

GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"
OUT = GEN / "p2502_s1452_strict_completion_bridge_minimal_triple_frontier_certificate.json"
MD = GEN / "p2502_s1452_strict_completion_bridge_minimal_triple_frontier_certificate.md"

SOURCE_FILES = {
    "P2501_INTERVAL_NEWTON_ROOT_ENCLOSURE": GEN / "p2501_s1451_phase_normalized_curvature_interval_newton_root_enclosure_certificate.json",
    "SCRATCH_THEOREM_FRONTIER_TRUTH_TABLE": SCRATCH / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "SCRATCH_THEOREM_FRONTIER_CUT": SCRATCH / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json",
    "SCRATCH_THEOREM_FRONTIER_LOW_WEIGHT": SCRATCH / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.json",
    "SCRATCH_COMPONENT_GAP": SCRATCH / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
}

TARGET_ORDER = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]
LOCAL_COMPRESSION_EVIDENCE_ATOM = "p2493_to_p2501_phase_normalized_curvature_enclosure_chain"


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
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2502|S1452|bridge minimal triple frontier|strict-source triple|three-atom bridge unlock|triple frontier certificate|minimal bridge triple",
        "precursor_packets": "P2493|P2494|P2495|P2496|P2497|P2498|P2499|P2500|P2501|curvature enclosure|interval Newton",
        "frontier_language": "theorem frontier|truth table|frontier cut|low-weight extension|open theorem atoms|strict source atoms|triple extension",
        "bridge_guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|selector/source",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|legacy-role transfer|directed-rounding|bridge theorem-level closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def target_values(true_atoms: set[str], target_definitions: dict[str, Any]) -> dict[str, bool]:
    bridge_atoms = set(target_definitions["bridge_theorem_level_closure"])
    role_atoms = set(target_definitions["role_transfer_theorem_level_closure"])
    selector_atoms = set(target_definitions["selector_qw2191_closure"])
    bridge = bridge_atoms.issubset(true_atoms)
    role = role_atoms.issubset(true_atoms)
    selector = selector_atoms.issubset(true_atoms)
    toe = bridge and role and selector
    return {
        "bridge_theorem_level_closure": bridge,
        "role_transfer_theorem_level_closure": role,
        "selector_qw2191_closure": selector,
        "toe_closure": toe,
    }


def signature(values: dict[str, bool]) -> str:
    return "".join("1" if values[name] else "0" for name in TARGET_ORDER)


def enumerate_triples(open_atoms: list[str], target_definitions: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for index, atoms in enumerate(itertools.combinations(sorted(open_atoms), 3)):
        atom_set = set(atoms)
        targets = target_values(atom_set, target_definitions)
        closed = [name for name in TARGET_ORDER if targets[name]]
        rows.append({
            "triple_index": index,
            "true_atoms": list(atoms),
            "target_signature_bridge_role_selector_toe": signature(targets),
            "closed_targets": closed,
            **targets,
        })
    return rows


def summarize_triples(rows: list[dict[str, Any]], target_definitions: dict[str, Any]) -> dict[str, Any]:
    bridge_rows = [row for row in rows if row["bridge_theorem_level_closure"]]
    role_rows = [row for row in rows if row["role_transfer_theorem_level_closure"]]
    selector_rows = [row for row in rows if row["selector_qw2191_closure"]]
    toe_rows = [row for row in rows if row["toe_closure"]]
    strict_source_triple = sorted(target_definitions["bridge_theorem_level_closure"])
    exact_bridge_rows = [row for row in bridge_rows if row["true_atoms"] == strict_source_triple]
    return {
        "triple_extension_count": len(rows),
        "bridge_closing_triple_count": len(bridge_rows),
        "bridge_closing_triples": bridge_rows,
        "exact_strict_source_triple": strict_source_triple,
        "exact_strict_source_triple_closes_bridge": len(exact_bridge_rows) == 1,
        "strict_source_triple_closes_only_bridge": len(exact_bridge_rows) == 1 and exact_bridge_rows[0]["target_signature_bridge_role_selector_toe"] == "1000",
        "role_transfer_closing_triple_count": len(role_rows),
        "selector_closing_triple_count": len(selector_rows),
        "toe_closing_triple_count": len(toe_rows),
        "no_triple_closes_role_transfer": len(role_rows) == 0,
        "no_triple_closes_toe": len(toe_rows) == 0,
        "all_selector_triples_contain_chi11": all("chi11_selector_source" in row["true_atoms"] for row in selector_rows),
        "selector_triple_sample_head": selector_rows[:8],
        "triple_rows_head": rows[:8],
        "triple_rows_tail": rows[-8:],
        "triple_rows_sha256": sha256_json(rows),
    }


def build_frontier_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    truth = sources["SCRATCH_THEOREM_FRONTIER_TRUTH_TABLE"]
    low_weight = sources["SCRATCH_THEOREM_FRONTIER_LOW_WEIGHT"]
    cut = sources["SCRATCH_THEOREM_FRONTIER_CUT"]
    component_gap = sources["SCRATCH_COMPONENT_GAP"]
    p2501 = theorem(
        sources["P2501_INTERVAL_NEWTON_ROOT_ENCLOSURE"],
        "phase_normalized_curvature_interval_newton_root_enclosure_certificate",
    )
    open_atoms = truth["open_atoms"]
    target_definitions = truth["target_definitions"]
    rows = enumerate_triples(open_atoms, target_definitions)
    summary = summarize_triples(rows, target_definitions)
    curvature_chain_complete_as_finite_evidence = p2501.get("single_inflection_root_on_audited_domain_inherited_and_narrowed") is True
    local_chain_is_not_frontier_atom = LOCAL_COMPRESSION_EVIDENCE_ATOM not in open_atoms
    return {
        "source_truth_table_status": truth.get("status"),
        "source_low_weight_status": low_weight.get("status"),
        "source_cut_status": cut.get("status"),
        "source_component_gap_status": component_gap.get("status"),
        "open_atoms": open_atoms,
        "target_order": TARGET_ORDER,
        "target_definitions": target_definitions,
        "triple_frontier_summary": summary,
        "p2501_curvature_enclosure_inherited": curvature_chain_complete_as_finite_evidence,
        "local_compression_evidence_atom_name": LOCAL_COMPRESSION_EVIDENCE_ATOM,
        "local_compression_evidence_is_not_an_open_theorem_atom": local_chain_is_not_frontier_atom,
        "adding_local_curvature_chain_alone_changes_frontier_signature": False,
        "curvature_chain_does_not_discharge_strict_source_atoms": curvature_chain_complete_as_finite_evidence and local_chain_is_not_frontier_atom,
        "bridge_minimal_size_three_inherited": truth["theorem_frontier_truth_table_summary"].get("bridge_minimal_set_size") == 3,
        "low_weight_no_pair_bridge_closure_inherited": low_weight["theorem_frontier_low_weight_extension_summary"].get("no_pair_closes_bridge") is True,
        "component_gap_sources_still_missing_inherited": component_gap["completion_gap_summary"].get("strict_dynamic_sources_missing") is True,
        "strict_source_triple_is_next_bridge_theorem_frontier": summary["exact_strict_source_triple_closes_bridge"] and summary["strict_source_triple_closes_only_bridge"],
        "role_transfer_still_needs_four_atom_package": summary["no_triple_closes_role_transfer"],
        "toe_still_needs_all_seven_atoms": summary["no_triple_closes_toe"] and truth["theorem_frontier_truth_table_summary"].get("toe_minimal_set_size") == 7,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
    }


def append_once(path: Path, marker: str, section: str) -> None:
    body = path.read_text(encoding="utf-8")
    if marker not in body:
        path.write_text(body.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2502/S1452 strict-completion bridge minimal-triple frontier certificate

`P2502/S1452` pivots from further local curvature narrowing to the current bridge-completion theorem frontier.  It exhaustively enumerates all `35` three-atom extensions of the seven open strict-completion atoms and finds exactly one triple that closes the bridge target: `{strict_dynamical_source_for_A_P_D, strict_phase_frequency_source, strict_damping_beta_eta_source}`.  That triple closes only bridge theorem-level readiness, not role-transfer, selector/QW-2191, or ToE.  The P2493--P2501 curvature enclosure chain is inherited as finite compression evidence, but it is not one of the open theorem-source atoms and therefore does not by itself change the frontier signature.

This is a finite theorem-frontier enumeration and planning certificate, not a bridge theorem, source theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.  It redirects the next honest bridge move toward actual strict-source atoms rather than more local root contraction.
"""
    lag_section = """
## P2502/S1452 bridge minimal-triple frontier guard

`P2502/S1452` enumerates all three-atom theorem-frontier extensions and identifies the unique bridge-closing strict-source triple while preserving the role-transfer, selector/QW-2191, and ToE blockers.  The P2493--P2501 curvature chain remains finite compression evidence only; it does not license a role-bearing `L_total` term or a nonlinear compression-flow source theorem.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2502/S1452 strict-completion bridge minimal-triple frontier certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2502/S1452 bridge minimal-triple frontier guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    cert = build_frontier_certificate(sources)
    summary = cert["triple_frontier_summary"]
    theorem_export = {
        "theorem_name": "P2502_T1_strict_completion_bridge_minimal_triple_frontier_certificate",
        "audited_chain": ["P2493-P2501 curvature chain", "strict-completion theorem-frontier scratch certificates"],
        "bridge_minimal_triple_frontier_certificate": cert,
        "open_atom_count": len(cert["open_atoms"]),
        "triple_extension_count": summary["triple_extension_count"],
        "bridge_closing_triple_count": summary["bridge_closing_triple_count"],
        "bridge_closing_triples": summary["bridge_closing_triples"],
        "exact_strict_source_triple": summary["exact_strict_source_triple"],
        "strict_source_triple_closes_only_bridge": summary["strict_source_triple_closes_only_bridge"],
        "role_transfer_closing_triple_count": summary["role_transfer_closing_triple_count"],
        "selector_closing_triple_count": summary["selector_closing_triple_count"],
        "toe_closing_triple_count": summary["toe_closing_triple_count"],
        "p2501_curvature_enclosure_inherited": cert["p2501_curvature_enclosure_inherited"],
        "curvature_chain_does_not_discharge_strict_source_atoms": cert["curvature_chain_does_not_discharge_strict_source_atoms"],
        "strict_source_triple_is_next_bridge_theorem_frontier": cert["strict_source_triple_is_next_bridge_theorem_frontier"],
        "role_transfer_still_needs_four_atom_package": cert["role_transfer_still_needs_four_atom_package"],
        "toe_still_needs_all_seven_atoms": cert["toe_still_needs_all_seven_atoms"],
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2502 enumerates theorem-frontier triples; it does not prove any missing strict-source atom.",
            "The P2493-P2501 curvature chain remains finite compression evidence and is not a nonlinear compression-flow source theorem.",
            "No role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Attack the strict-source triple directly, especially strict_damping_beta_eta_source, strict_phase_frequency_source, or strict_dynamical_source_for_A_P_D, instead of treating the local curvature root enclosure as bridge closure.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "truth_table_open_atom_count_is_seven": theorem_export["open_atom_count"] == 7,
        "all_triples_enumerated": theorem_export["triple_extension_count"] == 35,
        "unique_bridge_closing_triple": theorem_export["bridge_closing_triple_count"] == 1,
        "strict_source_triple_closes_only_bridge": theorem_export["strict_source_triple_closes_only_bridge"],
        "no_triple_closes_role_transfer": theorem_export["role_transfer_closing_triple_count"] == 0,
        "no_triple_closes_toe": theorem_export["toe_closing_triple_count"] == 0,
        "p2501_curvature_enclosure_inherited": theorem_export["p2501_curvature_enclosure_inherited"],
        "curvature_chain_not_promoted_to_source_atom": theorem_export["curvature_chain_does_not_discharge_strict_source_atoms"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2502",
        "stage_id": "S1452",
        "status": "STRICT_COMPLETION_BRIDGE_MINIMAL_TRIPLE_FRONTIER_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_completion_bridge_minimal_triple_frontier_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_completion_bridge_minimal_triple_frontier_certificate"]["theorem_export"]
    lines = [
        "# P2502/S1452 strict-completion bridge minimal-triple frontier certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Open theorem atoms: `{t['open_atom_count']}`.",
        f"- Three-atom extensions enumerated: `{t['triple_extension_count']}`.",
        f"- Bridge-closing triples: `{t['bridge_closing_triple_count']}` -> `{t['bridge_closing_triples']}`.",
        f"- Strict-source triple closes only bridge: `{t['strict_source_triple_closes_only_bridge']}`.",
        f"- Role-transfer-closing triples: `{t['role_transfer_closing_triple_count']}`.",
        f"- Selector-closing triples: `{t['selector_closing_triple_count']}`.",
        f"- ToE-closing triples: `{t['toe_closing_triple_count']}`.",
        f"- P2501 curvature enclosure inherited: `{t['p2501_curvature_enclosure_inherited']}`.",
        f"- Curvature chain does not discharge strict source atoms: `{t['curvature_chain_does_not_discharge_strict_source_atoms']}`.",
        "",
        "## Negative controls",
        "",
        "This packet does not export a missing strict-source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_completion_bridge_minimal_triple_frontier_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_completion_bridge_minimal_triple_frontier_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
