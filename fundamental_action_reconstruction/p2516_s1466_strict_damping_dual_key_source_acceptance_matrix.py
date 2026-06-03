#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2516_s1466_strict_damping_dual_key_source_acceptance_matrix.json"
MD = GEN / "p2516_s1466_strict_damping_dual_key_source_acceptance_matrix.md"

SOURCE_FILES = {
    "P2414_PARAMETER_IDENTIFIABILITY": GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json",
    "P2515_OPERATOR_SIGNATURE_ACCEPTANCE": GEN / "p2515_s1465_strict_damping_rg_operator_order_signature_acceptance_audit.json",
}

ATOMS = ["beta_eta_numeric_source", "m2_operator_signature_source"]


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
        "new_packet": "P2516|S1466|dual-key source acceptance|strict damping source acceptance matrix|beta eta operator signature|two-key damping source",
        "precursor_packets": "P2414|S1364|strict damping parameter identifiability|P2515|operator-order signature acceptance|P2514|higher-order selector nonidentifiability",
        "bridge_language": "legacy -> strict completion bridge|damping-compression bridge|beta_tors.*beta/eta|operator signature|biharmonic",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def boolean_truth_table() -> dict[str, Any]:
    rows = []
    for mask in range(1 << len(ATOMS)):
        active = [ATOMS[index] for index in range(len(ATOMS)) if mask & (1 << index)]
        beta_eta = "beta_eta_numeric_source" in active
        m2_signature = "m2_operator_signature_source" in active
        strict_source = beta_eta and m2_signature
        rows.append({
            "mask": format(mask, "02b"),
            "active_atoms": active,
            "beta_eta_numeric_source": beta_eta,
            "m2_operator_signature_source": m2_signature,
            "strict_damping_beta_eta_source_acceptance": strict_source,
            "missing_atoms": [atom for atom in ATOMS if atom not in active],
            "interpretation": "accepted only when both numeric beta/eta source and m=2 operator signature source are present" if strict_source else "proper subset obstruction: source acceptance remains blocked",
        })
    accepted = [row for row in rows if row["strict_damping_beta_eta_source_acceptance"]]
    rejected = [row for row in rows if not row["strict_damping_beta_eta_source_acceptance"]]
    return {
        "atoms": ATOMS,
        "boolean_normal_form": "strict_damping_beta_eta_source = beta_eta_numeric_source AND m2_operator_signature_source",
        "anf_over_gf2": "beta_eta_numeric_source * m2_operator_signature_source",
        "rows": rows,
        "accepted_row_count": len(accepted),
        "rejected_proper_subset_count": len(rejected),
        "all_proper_subsets_rejected": all(len(row["active_atoms"]) < len(ATOMS) and not row["strict_damping_beta_eta_source_acceptance"] for row in rejected),
        "unique_minimal_accepting_set": accepted[0]["active_atoms"] if len(accepted) == 1 else [],
    }


def proper_subset_obstruction_lattice() -> dict[str, Any]:
    subsets = []
    for size in range(len(ATOMS)):
        for combo in combinations(ATOMS, size):
            missing = [atom for atom in ATOMS if atom not in combo]
            subsets.append({
                "active_atoms": list(combo),
                "missing_atoms": missing,
                "strict_damping_beta_eta_source_acceptance": False,
                "obstruction": f"missing {', '.join(missing)}",
            })
    return {
        "proper_subset_rows": subsets,
        "proper_subset_count": len(subsets),
        "every_proper_subset_blocked": all(not row["strict_damping_beta_eta_source_acceptance"] for row in subsets),
    }


def build_dual_key_certificate(p2414: dict[str, Any], p2515: dict[str, Any]) -> dict[str, Any]:
    truth = boolean_truth_table()
    lattice = proper_subset_obstruction_lattice()
    p2414_identifies_numeric = p2414.get("strict_beta_eta_identified_from_samples") is True and p2414.get("strict_beta_eta_source_exported") is False
    p2515_identifies_signature = p2515.get("p2506_roughness_m2_signature_identified_not_sourced") is True and p2515.get("strict_damping_beta_eta_source_exported") is False
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2414_numeric_beta_eta_identified_but_unsourced": p2414_identifies_numeric,
        "p2515_m2_operator_signature_identified_but_unsourced": p2515_identifies_signature,
        "dual_key_atoms": ATOMS,
        "boolean_truth_table": truth,
        "proper_subset_obstruction_lattice": lattice,
        "dual_key_acceptance_normal_form_exported": truth["accepted_row_count"] == 1 and truth["all_proper_subsets_rejected"] and lattice["every_proper_subset_blocked"],
        "numeric_key_alone_insufficient": True,
        "operator_signature_key_alone_insufficient": True,
        "both_keys_required_for_future_source_theorem": True,
        "current_repo_has_both_keys_as_identified_targets_not_sources": p2414_identifies_numeric and p2515_identifies_signature,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2516/S1466 strict damping dual-key source acceptance matrix

`P2516/S1466` combines the P2414 numeric damping result with the P2515 operator-signature result into a two-key source acceptance normal form.  P2414 identifies the strict denominator target `beta=1, eta=9/5` from accepted samples and proves it is not legacy linear `beta_tors` absorption, but it does not source the numbers.  P2515 identifies the P2506 roughness selector as the `m=2` biharmonic/fourth-order operator signature, but it does not source that signature.  The acceptance matrix is therefore `strict_damping_beta_eta_source = beta_eta_numeric_source AND m2_operator_signature_source`; every proper subset is rejected.

This is a sharper source-target certificate, not a source theorem.  It blocks two common false closures: numeric `beta/eta` without an operator signature, and an operator signature without strict numeric damping source.  It exports no damping-compression bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2516/S1466 dual-key strict damping source guard

`P2516/S1466` records that a future strict damping action/source must provide two independent keys: the numeric strict damping target `beta=1, eta=9/5` and the `m=2` biharmonic operator signature.  Either key alone is a proper-subset obstruction, so the current conditional selector chain still cannot license a nonlinear compression-flow source theorem or role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2516/S1466 strict damping dual-key source acceptance matrix", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2516/S1466 dual-key strict damping source guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2414 = theorem(sources["P2414_PARAMETER_IDENTIFIABILITY"], "strict_damping_parameter_identifiability_nonabsorption_certificate")
    p2515 = theorem(sources["P2515_OPERATOR_SIGNATURE_ACCEPTANCE"], "strict_damping_rg_operator_order_signature_acceptance_audit")
    cert = build_dual_key_certificate(p2414, p2515)
    truth = cert["boolean_truth_table"]
    lattice = cert["proper_subset_obstruction_lattice"]
    theorem_export = {
        "theorem_name": "P2516_T1_strict_damping_dual_key_source_acceptance_matrix",
        "audited_chain": ["P2414/S1364", "P2515/S1465"],
        "strict_damping_dual_key_source_acceptance_matrix": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2414_numeric_beta_eta_identified_but_unsourced": cert["p2414_numeric_beta_eta_identified_but_unsourced"],
        "p2515_m2_operator_signature_identified_but_unsourced": cert["p2515_m2_operator_signature_identified_but_unsourced"],
        "boolean_normal_form": truth["boolean_normal_form"],
        "accepted_row_count": truth["accepted_row_count"],
        "rejected_proper_subset_count": truth["rejected_proper_subset_count"],
        "all_proper_subsets_rejected": truth["all_proper_subsets_rejected"],
        "unique_minimal_accepting_set": truth["unique_minimal_accepting_set"],
        "proper_subset_count": lattice["proper_subset_count"],
        "dual_key_acceptance_normal_form_exported": cert["dual_key_acceptance_normal_form_exported"],
        "numeric_key_alone_insufficient": cert["numeric_key_alone_insufficient"],
        "operator_signature_key_alone_insufficient": cert["operator_signature_key_alone_insufficient"],
        "both_keys_required_for_future_source_theorem": cert["both_keys_required_for_future_source_theorem"],
        "current_repo_has_both_keys_as_identified_targets_not_sources": cert["current_repo_has_both_keys_as_identified_targets_not_sources"],
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2516 exports a dual-key acceptance normal form, not a strict damping source theorem.",
            "Numeric beta/eta without an m=2 operator signature is rejected.",
            "An m=2 operator signature without strict numeric beta/eta sourcing is rejected.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Derive both beta=1, eta=9/5 and the m=2 biharmonic operator signature from strict nadsoliton dynamics, or explicitly mark either missing key as an axiom/non-strict premise.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2414_numeric_target_inherited": theorem_export["p2414_numeric_beta_eta_identified_but_unsourced"],
        "p2515_signature_target_inherited": theorem_export["p2515_m2_operator_signature_identified_but_unsourced"],
        "dual_key_normal_form_ok": theorem_export["dual_key_acceptance_normal_form_exported"],
        "proper_subsets_rejected": theorem_export["all_proper_subsets_rejected"] and theorem_export["proper_subset_count"] == 3,
        "source_not_exported": not theorem_export["strict_damping_beta_eta_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "strict_damping_beta_eta_source_exported",
            "damping_compression_bridge_component_ready",
            "full_bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "role_bearing_ltotal_exported",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2516",
        "stage_id": "S1466",
        "status": "STRICT_DAMPING_DUAL_KEY_SOURCE_ACCEPTANCE_MATRIX_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_dual_key_source_acceptance_matrix": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_dual_key_source_acceptance_matrix"]["theorem_export"]
    lines = [
        "# P2516/S1466 strict damping dual-key source acceptance matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2414 numeric beta/eta identified but unsourced: `{t['p2414_numeric_beta_eta_identified_but_unsourced']}`.",
        f"- P2515 m=2 operator signature identified but unsourced: `{t['p2515_m2_operator_signature_identified_but_unsourced']}`.",
        f"- Boolean normal form: `{t['boolean_normal_form']}`.",
        f"- Accepted row count: `{t['accepted_row_count']}`.",
        f"- Proper subset count: `{t['proper_subset_count']}`.",
        f"- All proper subsets rejected: `{t['all_proper_subsets_rejected']}`.",
        f"- Unique minimal accepting set: `{t['unique_minimal_accepting_set']}`.",
        f"- Strict source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports a source-acceptance normal form only. It does not derive either key from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_dual_key_source_acceptance_matrix']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_dual_key_source_acceptance_matrix"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
