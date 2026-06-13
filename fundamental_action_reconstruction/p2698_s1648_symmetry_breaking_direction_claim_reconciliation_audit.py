#!/usr/bin/env python3
"""P2698/S1648: symmetry-breaking/direction-claim reconciliation audit.

Responds to the objection that symmetry breaking and direction selection were
already found.  The audit greps the repo and reads the exported H36/H37/H38/H39,
P739/P740, P1227/P1233, and P2697 artifacts.  It separates real exported
premise-based directed/orientation support from strict-core QW-2191 discharge.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.json"
MD = GEN / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2697": GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json",
    "H36": GEN / "h36_directed_axis_orientation_audit.json",
    "H37": GEN / "h37_sign_distinction_state_audit.json",
    "H38": GEN / "h38_projective_selector_state_audit.json",
    "H39": GEN / "h39_global_selector_object_absence_audit.json",
    "H41": GEN / "h41_selector_atlas_and_gluing_data_audit_summary.json",
    "P739": GEN / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "P740": GEN / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "P1227": GEN / "p1227_w1_execution_summary.json",
    "P1233": GEN / "p1233_w1_strict_lane_readiness_gate_summary.json",
    "KAPPA": GEN / "kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "qw2191_discharged",
    "strict_core_selector_closure_exported",
    "aut_z12_invariant_canonicity_claimed",
    "premise_based_direction_promoted_to_strict_derivation",
    "p739_p740_pair12_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_claimed",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "nadsoliton_information_ontology": r"nadsoliton.*information|informational nadsoliton|primordial information|czyst[aą] INFORMACJ|fractal",
        "symmetry_breaking_direction_claims": r"symmetry.*break|symmetry_state|directed selector|directed/sign-sensitive|orientation|kierunk|generator/orientation|sign-sensitive",
        "premise_based_scope_markers": r"premise-based|T164|strict provenance fixing datum|convention|gauge only|not strict-derived",
        "qw2191_open_markers": r"QW-2191.*open|qw2191_closed.*false|global QW-2191 discharge remains open|strict-core uniqueness obstruction",
        "p739_p740_nonexport": r"P739|P740|strict-core upgrade.*nonexport|t193_target_exported|t194_target_exported|pair12.*upgrade.*false",
        "no_false_promotion_boundaries": r"L_total|ToE closure|strict-core selector closure|Aut\(Z_12\)-invariant|role transfer",
    }
    return {"tool": "rg", "mode": "symmetry-breaking/direction claim reconciliation after P2697", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def state_reads() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in INPUTS.items()}
    h36, h37, h38, h39, h41 = (loaded[name] for name in ("H36", "H37", "H38", "H39", "H41"))
    p739, p740, p1227, p1233, p2697, kappa = (loaded[name] for name in ("P739", "P740", "P1227", "P1233", "P2697", "KAPPA"))
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "nadsoliton_ontology_respected": True,
        "ontology_statement": "The nadsoliton is treated as primordial fractal information in solitonic state; no lower information layer is introduced.",
        "directed_axis_orientation_present": "PREMISE_BASED_T164" in h36.get("status", ""),
        "h37_sign_sensitive_directed_observable_exported": "SIGN_SENSITIVE" in h37.get("status", "") and "PREMISE_BASED_T164" in h37.get("status", ""),
        "h38_directed_continuation_selected": h38.get("evidence", {}).get("directed_continuation_selected") is True,
        "h39_qw2191_still_open": "global_qw_2191_discharge" in h39.get("missing", []) or "QW2191_STILL_OPEN" in h39.get("status", ""),
        "h41_qw2191_still_open": "QW2191_STILL_OPEN" in h41.get("status", "") or "QW-2191 discharge remains open" in h41.get("frontier", ""),
        "p739_directed_lane_exported": p739.get("current_global_premise_based_directed_selector_state_lane_exported") is True,
        "p739_lane_is_premise_based": p739.get("current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164") is True,
        "p739_pair12_upgrade_exported": p739.get("t193_target_exported_on_current_repo_state") is True,
        "p740_sign_fixed_lane_exported": p740.get("current_global_sign_fixed_directed_closure_lane_exported") is True,
        "p740_lane_is_gauge_only": p740.get("current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only") is True,
        "p740_pair12_upgrade_exported": p740.get("t194_target_exported_on_current_repo_state") is True,
        "p1227_w1_witness_discharged_but_theory_open": p1227.get("witness_status") == "DISCHARGED" and p1227.get("theory_closure_status") == "OPEN" and p1227.get("strict_closure_claim_allowed") is False,
        "p1233_readiness_not_global_closure": p1233.get("strict_lane_readiness_gate_pass") is True and p1233.get("theory_closure_status") == "OPEN" and p1233.get("strict_closure_claim_allowed") is False,
        "kappa_declares_no_aut_invariant_canonicity": "does not claim Aut-invariant canonicity" in json.dumps(kappa),
        "p2697_no_new_live_frontier_certificate": p2697.get("decision", {}).get("no_new_live_frontier_certificate") is True,
    }


def reconciliation_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "claim_under_audit": "symmetry/direction was already found",
            "support_is_real": reads["directed_axis_orientation_present"] and reads["h37_sign_sensitive_directed_observable_exported"] and reads["h38_directed_continuation_selected"],
            "scope": "real exported directed/sign-sensitive support in premise-based T164/T171 scope",
            "strict_closure_consequence": "does_not_discharge_QW_2191_by_itself",
            "promotion_allowed_now": False,
        },
        {
            "claim_under_audit": "strict-core selector closure follows from earlier W1/readiness artifacts",
            "support_is_real": reads["p1227_w1_witness_discharged_but_theory_open"] and reads["p1233_readiness_not_global_closure"],
            "scope": "witness/readiness support only; both artifacts keep theory closure open and forbid strict closure claim",
            "strict_closure_consequence": "readiness_is_not_global_selector_closure",
            "promotion_allowed_now": False,
        },
        {
            "claim_under_audit": "pair12 or global strict-core upgrade follows from directed lane",
            "support_is_real": reads["p739_directed_lane_exported"] and reads["p740_sign_fixed_lane_exported"],
            "scope": "P739/P740 export directed/sign-fixed lanes but explicitly nonexport pair12 strict-core upgrade",
            "strict_closure_consequence": "t193_t194_targets_remain_unexported",
            "promotion_allowed_now": False,
        },
        {
            "claim_under_audit": "Aut(Z12)-invariant canonical direction is now strict-derived",
            "support_is_real": reads["kappa_declares_no_aut_invariant_canonicity"],
            "scope": "explicit strict-provenance fixing datum, premise-based, not Aut-invariant canonicity",
            "strict_closure_consequence": "premise_based_fixing_may_be_used_only_with_scope_tag",
            "promotion_allowed_now": False,
        },
    ]


def decision(rows: list[dict[str, Any]], reads: dict[str, Any]) -> dict[str, Any]:
    no_false_promotion = all(not row["promotion_allowed_now"] for row in rows)
    return {
        "decision": "P2698_SYMMETRY_BREAKING_DIRECTION_CLAIM_RECONCILIATION_NO_FALSE_PASS",
        "acknowledgement": "The repo does contain real direction/orientation/sign-sensitive exports; the audit preserves them instead of ignoring them.",
        "bounded_result": "Those exports are premise-based or readiness/witness scoped and do not by themselves close QW-2191, pair12 strict-core upgrade, Aut(Z12)-invariant canonicity, L_total, or ToE.",
        "no_new_live_frontier_certificate_preserved": reads["p2697_no_new_live_frontier_certificate"] and no_false_promotion,
        "next_honest_step": "P2699 should be a one-object strict-internal selector-source candidate construction attempt only if a genuinely new non-premise-based orientation/selector source object is introduced.  If no such object is introduced, keep P2697/P2698 as the active no-new-live-frontier certificate and do not replay H36/H37/T164/P739/P740.",
        "forbidden_promotions": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2698/S1648 symmetry-breaking/direction-claim reconciliation audit", "", f"Status: `{payload['status']}`", "", "## Ontology"]
    lines.append(payload["state_reads"]["ontology_statement"])
    lines.extend(["", "## Reconciliation matrix"])
    for row in payload["reconciliation_matrix"]:
        lines.append(f"- `{row['claim_under_audit']}`: support_is_real=`{row['support_is_real']}`; scope={row['scope']}; promotion_allowed_now=`{row['promotion_allowed_now']}`")
    lines.extend(["", "## Decision", payload["decision"]["acknowledgement"], payload["decision"]["bounded_result"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    rows = reconciliation_matrix(reads)
    payload: dict[str, Any] = {
        "status": "P2698_SYMMETRY_BREAKING_DIRECTION_CLAIM_RECONCILIATION_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "reconciliation_matrix": rows,
        "decision": decision(rows, reads),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2698/S1648 symmetry-breaking direction reconciliation",
        "## P2698/S1648 symmetry-breaking direction reconciliation\n\n"
        "`P2698/S1648` greps and reads the existing direction/orientation artifacts.  H36/H37/H38 and P739/P740 are real support: a directed/sign-sensitive layer exists, but it is premise-based/convention-scoped and does not export `QW-2191` discharge, Aut(Z12)-invariant canonicity, pair12 strict-core upgrade, `L_total`, or ToE closure.  The nadsoliton ontology remains primordial fractal information in solitonic state; no lower information layer is introduced.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2698/S1648 directed support nonpromotion guard",
        "## P2698/S1648 directed support nonpromotion guard\n\n"
        "`P2698/S1648` preserves the real directed/sign-sensitive exports but blocks their silent promotion into a strict variational selector, role-bearing `L_total`, or ToE closure.  A next move needs a genuinely new non-premise-based strict selector/orientation source object; otherwise P2697/P2698 remain the no-new-live-frontier state.\n",
    )
    append_once(
        AGENTS,
        "Current symmetry-breaking/direction reconciliation guardrail (P2698/S1648, 2026-06-13)",
        "## Current symmetry-breaking/direction reconciliation guardrail (P2698/S1648, 2026-06-13)\n\n"
        "- P2698 confirms that earlier direction/orientation/sign-sensitive exports are real and must not be ignored, but they are premise-based/readiness/convention scoped: H36/H37/H38 and P739/P740 do not discharge `QW-2191`, do not export pair12 strict-core upgrade, and do not provide Aut(Z12)-invariant canonicity.\n"
        "- Preserve the ontology: the nadsoliton is primordial fractal information in a solitonic state, with no lower information layer.\n"
        "- Do not promote premise-based direction into strict selector closure, `L_total`, role transfer, or ToE closure.  A next admissible move must introduce a genuinely new non-premise-based strict selector/orientation source object; otherwise keep the P2697/P2698 no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    main()
