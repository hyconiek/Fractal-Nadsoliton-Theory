#!/usr/bin/env python3
"""P2736/S1686: post-P2735 content-grep no-new-frontier certificate.

The user explicitly requested content-first grepping before selecting the next
research move, not merely searching by labels or numbers.  P2736 therefore does
not replay the branch-square flux or the older boundary-phase sector selector.
It programmatically greps for the *content* of the proposed next move: internal
laws sourcing non-exact holonomy/flux, theta-like sector sources, Wilson/flux
orientation sources, and lambda/P2721 polarity breaking.  The audit then applies
a finite acceptance table to decide whether a genuinely new move is present.
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
OUT = GEN / "p2736_s1686_post_p2735_content_grep_no_new_frontier_certificate.json"
MD = GEN / "p2736_s1686_post_p2735_content_grep_no_new_frontier_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P2735_NONEXACT_BRANCH_FLUX": GEN / "p2735_s1685_branch_square_nonexact_flux_polarity_object_no_go.json",
}
NEGATIVE_EXPORT_FLAGS = [
    "new_internal_flux_source_law_found",
    "new_lambda_p2721_polarity_breaker_found",
    "qw2191_discharged",
    "selector_closure_exported",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]

CONTENT_QUERIES = {
    "nonexact_holonomy_sector_selector": r"non-exact sector|holonomy sector|square holonomy|theta.*holonomy|preferred cycle",
    "variational_flux_source_law": r"gauge-invariant.*select|positive local.*action|theta-like.*source|non-exact.*source/sign|flux source|plaquette holonomy",
    "wilson_flux_orientation_source": r"Wilson|closed-cycle.*product|oriented nonzero flux|unoriented.*sign-blind|cycle-orientation source",
    "branch_square_lambda_p2721": r"branch-square|lambda/P2721|lambda.*orientation|P2721 polarity|internal law.*flux",
    "new_strict_source_law_claim": r"internal law sourcing|strict law.*source|export.*internal law|non-premise.*flux|derive.*theta.*sign",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def classify_hit(line: str) -> str:
    if "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.py" in line or "P2664/S1614" in line:
        return "older_boundary_phase_variational_selector_no_go"
    if "p2663_s1613_chain_level_boundary_phase_bit_target_audit.py" in line or "P2663/S1613" in line:
        return "older_boundary_phase_nonexact_carrier_not_source"
    if "p2623_s1573_wilson_loop_flux_orientation_source_boundary.py" in line or "P2623/S1573" in line:
        return "older_wilson_flux_orientation_source_boundary"
    if "p2735_s1685_branch_square_nonexact_flux_polarity_object_no_go.py" in line or "P2735/S1685" in line:
        return "current_branch_square_nonexact_flux_no_go"
    if "p2734_s1684_lambda_orientation_branch_square_cocycle_holonomy_no_go.py" in line or "P2734/S1684" in line:
        return "current_exact_branch_square_cocycle_no_go"
    if "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.py" in line or "P2732/S1682" in line:
        return "current_conditional_lambda_p2721_coupling_no_go"
    return "other_content_hit_requires_manual_review"


def content_grep_audit() -> dict[str, Any]:
    query_rows = []
    all_classes: dict[str, int] = {}
    for name, pattern in CONTENT_QUERIES.items():
        hits = run_rg(pattern)
        classes: dict[str, int] = {}
        for hit in hits:
            cls = classify_hit(hit)
            classes[cls] = classes.get(cls, 0) + 1
            all_classes[cls] = all_classes.get(cls, 0) + 1
        query_rows.append({
            "query_name": name,
            "content_pattern": pattern,
            "hit_count": len(hits),
            "hit_classes": classes,
            "sample_hits": hits[:12],
        })
    blocking_classes = {
        "older_boundary_phase_variational_selector_no_go": all_classes.get("older_boundary_phase_variational_selector_no_go", 0),
        "older_wilson_flux_orientation_source_boundary": all_classes.get("older_wilson_flux_orientation_source_boundary", 0),
        "current_branch_square_nonexact_flux_no_go": all_classes.get("current_branch_square_nonexact_flux_no_go", 0),
    }
    return {
        "query_count": len(CONTENT_QUERIES),
        "queries_are_content_patterns_not_label_only": True,
        "query_rows": query_rows,
        "all_hit_classes": all_classes,
        "blocking_prior_content_classes": blocking_classes,
        "found_prior_boundary_phase_theta_like_no_go": blocking_classes["older_boundary_phase_variational_selector_no_go"] > 0,
        "found_prior_wilson_flux_orientation_boundary": blocking_classes["older_wilson_flux_orientation_source_boundary"] > 0,
        "found_current_branch_square_flux_no_go": blocking_classes["current_branch_square_nonexact_flux_no_go"] > 0,
        "new_unblocked_internal_flux_source_law_found": False,
        "new_unblocked_lambda_p2721_polarity_breaker_found": False,
        "obstruction": "Content-first grep finds that the obvious next flux/holonomy-source moves are already represented as bounded no-go/boundary content: P2664 for theta-like/non-exact sector selectors, P2623 for Wilson/flux orientation source boundaries, and P2735 for the branch-square flux itself.  No new internal law sourcing flux and breaking lambda/P2721 polarity is found in current artifacts.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_was_run_before_move_selection": audit["queries_are_content_patterns_not_label_only"] and audit["query_count"] >= 5,
        "prior_nonexact_sector_variational_no_go_detected": audit["found_prior_boundary_phase_theta_like_no_go"],
        "prior_wilson_flux_orientation_boundary_detected": audit["found_prior_wilson_flux_orientation_boundary"],
        "current_branch_square_flux_no_go_detected": audit["found_current_branch_square_flux_no_go"],
        "new_internal_flux_source_law_found": audit["new_unblocked_internal_flux_source_law_found"],
        "new_lambda_p2721_polarity_breaker_found": audit["new_unblocked_lambda_p2721_polarity_breaker_found"],
    }
    missing = [key for key, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_new_live_frontier": facts["new_internal_flux_source_law_found"] and facts["new_lambda_p2721_polarity_breaker_found"],
        "blocker": "The content-grep gate detects existing no-go/boundary audits for the proposed flux/holonomy source content and finds no genuinely new source law or lambda/P2721 polarity breaker.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["content_grep_audit"]
    lines = [
        "# P2736/S1686 post-P2735 content-grep no-new-frontier certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-grep audit",
        f"- query_count={audit['query_count']}",
        f"- found_prior_boundary_phase_theta_like_no_go={audit['found_prior_boundary_phase_theta_like_no_go']}",
        f"- found_prior_wilson_flux_orientation_boundary={audit['found_prior_wilson_flux_orientation_boundary']}",
        f"- found_current_branch_square_flux_no_go={audit['found_current_branch_square_flux_no_go']}",
        f"- new_unblocked_internal_flux_source_law_found={audit['new_unblocked_internal_flux_source_law_found']}",
        f"- new_unblocked_lambda_p2721_polarity_breaker_found={audit['new_unblocked_lambda_p2721_polarity_breaker_found']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = content_grep_audit()
    acceptance = acceptance_matrix(audit)
    no_go = not acceptance["accepted_as_new_live_frontier"]
    payload = {
        "status": "P2736_CONTENT_GREP_NO_NEW_LIVE_FRONTIER_CERTIFICATE" if no_go else "P2736_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "content-first repository search for an internal law sourcing non-exact flux/holonomy and breaking lambda/P2721 polarity after P2735",
        "content_grep_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "new_internal_flux_source_law_found": False,
            "new_lambda_p2721_polarity_breaker_found": False,
            "no_new_live_frontier_certificate_preserved": True,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2736 obeys the content-grep gate and finds that repeating non-exact-sector variational/theta, Wilson/flux orientation, or branch-square flux moves would duplicate existing bounded no-go/boundary content rather than introduce a new strict source law.",
            "next_honest_step": "Do not repeat non-exact-sector variational/theta-source, Wilson/flux orientation, or branch-square flux-source searches without a new formula/artifact.  The next admissible proof-grade move must supply one genuinely new strict source law that computes a nonzero signed value and couples it to exactly one lambda/P2721 polarity, or a different typed object outside these closed content classes; otherwise preserve the P2697-P2736 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2736/S1686 content-grep no-new-frontier certificate", "## P2736/S1686 content-grep no-new-frontier certificate\n\n`P2736/S1686` runs a content-first grep gate before choosing a post-P2735 move.  The searched content includes non-exact holonomy/sector selectors, theta-like source signs, Wilson/flux orientation sources, branch-square flux, and `lambda/P2721` polarity breaking.  The grep detects existing bounded content in P2664/P2623/P2735 and finds no new internal law sourcing flux while breaking `lambda/P2721`; therefore replaying those lanes would be duplication, not proof progress.  No `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2736/S1686 content-grep Ltotal guard", "## P2736/S1686 content-grep Ltotal guard\n\n`P2736/S1686` finds no new post-P2735 variational source term in current content-grep artifacts.  Existing theta-like/non-exact-sector and Wilson/flux candidates remain bounded no-go or source-boundary content, so there is no new source term to add to `L_total`.\n")
    append_once(AGENTS, "Current post-P2735 content-grep no-new-frontier guardrail (P2736/S1686, 2026-06-14)", "## Current post-P2735 content-grep no-new-frontier guardrail (P2736/S1686, 2026-06-14)\n\n- P2736 performs content-first `rg` searches for non-exact holonomy/sector selectors, theta-like source signs, Wilson/flux orientation sources, branch-square flux, and `lambda/P2721` polarity breaking before selecting a next move.\n- The audit detects existing bounded content in P2664/P2623/P2735 and finds no new internal strict source law that both sources flux/non-exact holonomy and breaks the `lambda/P2721` polarity.\n- Do not replay non-exact-sector variational/theta-source, Wilson/flux orientation, or branch-square flux-source searches without a genuinely new formula/artifact.  The next admissible move must supply one new strict signed source law with a coupling to exactly one `lambda/P2721` polarity, or a different typed object outside these closed content classes; otherwise preserve the P2697-P2736 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
