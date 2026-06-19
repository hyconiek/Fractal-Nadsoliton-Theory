#!/usr/bin/env python3
"""P2913/S1863: Gamma source-provenance obstruction audit.

P2912 supplied the finite variational Jacobian skeleton for the Gamma_9_5 lane.
The next honest proof-grade question is not another finite matrix: it is whether
current artifacts export a strict source/provenance theorem that gives the
unit-bearing value of Gamma_9_5.

This script constructs that missing theorem as an acceptance schema, builds a
small finite family of candidate source objects, and scans generated JSON
artifacts for positive provenance exports.  The result is an obstruction audit:
dimension and finite localization readiness exist, but no current artifact
exports a strict nadsoliton-derived action-unit source for Gamma_9_5.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2912 = GEN / "p2912_s1862_gamma_variational_chain_rule_skeleton_gate.json"
OUT = GEN / "p2913_s1863_gamma_source_provenance_obstruction_audit.json"
MD = GEN / "p2913_s1863_gamma_source_provenance_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

POSITIVE_KEYS = (
    "strict_gamma_9_5_source_exported",
    "strict_gamma_source_exported",
    "gamma_9_5_strict_provenance_exported",
    "gamma_9_5_source_theorem_exported",
    "unit_bearing_gamma_9_5_exported",
)
MENTION_TERMS = ("Gamma_9_5", "gamma_9_5", "Γ_9_5", "unit-bearing coupling-token")
NEGATIVE_TERMS = ("no strict", "not export", "unexported", "readiness", "candidate", "symbolic", "obstruction")


def candidate_source_objects() -> list[dict[str, Any]]:
    """Construct the finite candidate family for the missing Gamma source theorem."""
    return [
        {
            "name": "Gamma_9_5_imported_symbol",
            "construction": "reuse the P2910 symbol with required Action dimension",
            "dimension_action": True,
            "nonzero_value_supplied": False,
            "strict_nadsoliton_source_law": False,
            "compatible_with_p2911_p2912": True,
            "accepted": False,
            "failure": "symbol import gives dimensional readiness but no source law or value",
        },
        {
            "name": "U_9_5_symbolic_density_constant",
            "construction": "identify Gamma_9_5 with the earlier symbolic U_9_5 local derivative",
            "dimension_action": False,
            "nonzero_value_supplied": False,
            "strict_nadsoliton_source_law": False,
            "compatible_with_p2911_p2912": False,
            "accepted": False,
            "failure": "U_9_5 remains symbolic and was never exported as unit-bearing action token",
        },
        {
            "name": "Gamma_9_5_zero_token",
            "construction": "set Gamma_9_5 = 0 to avoid importing a value",
            "dimension_action": True,
            "nonzero_value_supplied": False,
            "strict_nadsoliton_source_law": False,
            "compatible_with_p2911_p2912": True,
            "accepted": False,
            "failure": "zero token kills the local density and is not a nonzero strict source",
        },
        {
            "name": "Gamma_9_5_free_scale_parameter",
            "construction": "allow an arbitrary positive action-scale parameter",
            "dimension_action": True,
            "nonzero_value_supplied": True,
            "strict_nadsoliton_source_law": False,
            "compatible_with_p2911_p2912": True,
            "accepted": False,
            "failure": "free scale is conventional/imported and not nadsoliton-derived",
        },
        {
            "name": "Gamma_9_5_pullback_normalized_sum",
            "construction": "normalize Gamma by total endpoint-average Jacobian weight",
            "dimension_action": False,
            "nonzero_value_supplied": True,
            "strict_nadsoliton_source_law": False,
            "compatible_with_p2911_p2912": True,
            "accepted": False,
            "failure": "finite weights are dimensionless and cannot source an action unit",
        },
        {
            "name": "Gamma_9_5_strict_action_unit_theorem_required",
            "construction": "a missing theorem deriving a nonzero action unit from strict nadsoliton data and coupling it to Lambda/Jacobian slots",
            "dimension_action": True,
            "nonzero_value_supplied": False,
            "strict_nadsoliton_source_law": False,
            "compatible_with_p2911_p2912": True,
            "accepted": False,
            "failure": "this is the exact missing theorem, not a current export",
        },
    ]


def json_files() -> list[Path]:
    return sorted(path for path in GEN.glob("*.json") if path.name < OUT.name)


def scan_generated_json() -> dict[str, Any]:
    files = json_files()
    mention_hits: list[dict[str, Any]] = []
    positive_hits: list[dict[str, Any]] = []
    for path in files:
        raw = path.read_text(encoding="utf-8")
        lowered = raw.lower()
        if any(term.lower() in lowered for term in MENTION_TERMS):
            mention_hits.append({
                "file": str(path.relative_to(ROOT)),
                "has_negative_boundary_language": any(term in lowered for term in NEGATIVE_TERMS),
            })
        try:
            data = json.loads(raw)
        except json.JSONDecodeError:
            continue
        encoded = json.dumps(data, sort_keys=True, ensure_ascii=False)
        for key in POSITIVE_KEYS:
            if f'"{key.lower()}": true' in encoded.lower():
                positive_hits.append({"file": str(path.relative_to(ROOT)), "key": key})
    return {
        "json_files_scanned": len(files),
        "mention_hit_count": len(mention_hits),
        "mention_hits": mention_hits,
        "positive_provenance_hit_count": len(positive_hits),
        "positive_provenance_hits": positive_hits,
    }


def build_payload(p2912: dict[str, Any]) -> dict[str, Any]:
    candidates = candidate_source_objects()
    scan = scan_generated_json()
    accepted = [candidate for candidate in candidates if candidate["accepted"]]
    return {
        "status": "P2913_GAMMA_SOURCE_PROVENANCE_OBSTRUCTION_AUDIT_NO_EXPORT",
        "input_hashes": {"P2912": hashlib.sha256(P2912.read_bytes()).hexdigest() if P2912.exists() else None},
        "constructed_theoretical_objects": {
            "missing_theorem_name": "Strict_Gamma_9_5_Action_Unit_Source_Theorem",
            "acceptance_schema": [
                "derive Gamma_9_5 from strict nadsoliton data rather than import it",
                "export a nonzero signed/action-unit value with Action dimension",
                "couple the value to the P2911 Lambda_edge_to_site localization",
                "couple the value to the P2912 finite Jacobian/variational slot",
                "avoid Xi/J/defect-placement replay and role-transfer/ToE promotion",
            ],
            "candidate_source_objects": candidates,
            "generated_artifact_scan": scan,
        },
        "acceptance_matrix": {
            "p2912_rechecked_variational_readiness": p2912.get("acceptance_matrix", {}).get("finite_variational_chain_rule_skeleton_constructed") is True,
            "candidate_source_object_count": len(candidates),
            "accepted_candidate_count": len(accepted),
            "generated_json_files_scanned": scan["json_files_scanned"],
            "gamma_mention_hit_count": scan["mention_hit_count"],
            "positive_gamma_source_provenance_hit_count": scan["positive_provenance_hit_count"],
            "strict_gamma_9_5_source_theorem_exported": False,
            "nonzero_unit_bearing_gamma_value_exported": False,
            "p2911_p2912_coupling_theorem_exported": False,
            "accepted_as_nonproxy_ltotal_coupling_source": False,
        },
        "decision": {
            "positive_witnesses": {
                "missing_theorem_acceptance_schema_constructed": True,
                "finite_candidate_family_constructed": True,
                "generated_artifact_provenance_scan_executed": True,
            },
            "negative_export_flags": {
                "strict_gamma_9_5_source_exported": False,
                "nonzero_unit_bearing_gamma_value_exported": False,
                "strict_action_unit_source_theorem_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2913 constructs the exact missing Gamma_9_5 source theorem as an acceptance schema and audits six finite candidate source objects plus generated JSON provenance.  Dimension/readiness candidates exist, but accepted candidates remain zero and the scan finds no positive strict Gamma_9_5 source/provenance export.  Therefore P2911/P2912 readiness still cannot be promoted to unit-bearing nonproxy L_total, EOM, Hamiltonian, bridge, role transfer, or ToE closure.",
            "next_honest_step": "Do not add another symbolic Gamma/U inventory.  The next proof-grade step must either supply one explicit strict nadsoliton-derived action-unit theorem computing a nonzero Gamma_9_5 and coupling it to the P2911/P2912 slots, or pivot to the other missing theorem: strict field-variable/continuum-measure provenance for the finite Jacobian.  If neither theorem is supplied, preserve no-new-live-frontier for the Gamma lane.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2913/S1863 Gamma source-provenance obstruction audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Source theorem gate",
        f"- candidate source objects: `{acc['candidate_source_object_count']}`",
        f"- accepted candidates: `{acc['accepted_candidate_count']}`",
        f"- generated JSON files scanned: `{acc['generated_json_files_scanned']}`",
        f"- Gamma mention hits: `{acc['gamma_mention_hit_count']}`",
        f"- positive Gamma source/provenance hits: `{acc['positive_gamma_source_provenance_hit_count']}`",
        f"- strict Gamma_9_5 source theorem exported: `{acc['strict_gamma_9_5_source_theorem_exported']}`",
        f"- accepted as nonproxy L_total coupling source: `{acc['accepted_as_nonproxy_ltotal_coupling_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2912))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2913/S1863 Gamma source-provenance obstruction audit", "## P2913/S1863 Gamma source-provenance obstruction audit\n\n`P2913/S1863` constructs the missing `Strict_Gamma_9_5_Action_Unit_Source_Theorem` as an acceptance schema and audits a finite six-object candidate family plus generated JSON provenance.  The audit finds `0` accepted source candidates and `0` positive strict `Gamma_9_5` source/provenance hits.  Thus P2911/P2912 localization/Jacobian readiness remains nonpromotable: no nonzero unit-bearing `Gamma_9_5`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2913/S1863 Gamma source-provenance `L_total` guard", "## P2913/S1863 Gamma source-provenance `L_total` guard\n\n`P2913/S1863` audits whether current artifacts export the strict action-unit source theorem needed for `Gamma_9_5`.  Six candidate source objects and generated provenance are checked, but no accepted strict nonzero unit-bearing source is found.  Therefore the P2911/P2912 pullback/Jacobian skeleton cannot yet become a nonproxy `L_total` term, EOM, Hamiltonian, role-transfer theorem, bridge closure, or ToE closure.\n")
    append_once(AGENTS, "Current Gamma source-provenance obstruction guardrail (P2913/S1863, 2026-06-19)", "## Current Gamma source-provenance obstruction guardrail (P2913/S1863, 2026-06-19)\n\n- P2913 constructs the missing `Strict_Gamma_9_5_Action_Unit_Source_Theorem` as an acceptance schema and audits six candidate source objects plus generated JSON provenance.\n- The audit finds `0` accepted candidates and `0` positive strict `Gamma_9_5` source/provenance hits; dimension/readiness evidence remains insufficient.\n- Do not promote `Gamma_9_5`, symbolic `U_9_5`, P2911 pullback readiness, or P2912 Jacobian readiness to nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without an explicit strict nonzero action-unit source theorem and coupling theorem.\n- A next admissible move in this lane must provide that strict source theorem or pivot to strict field-variable/continuum-measure provenance for the finite Jacobian; otherwise preserve no-new-live-frontier for the Gamma lane.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
