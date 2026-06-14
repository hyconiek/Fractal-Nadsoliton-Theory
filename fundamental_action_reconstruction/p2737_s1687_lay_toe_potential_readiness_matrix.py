#!/usr/bin/env python3
"""P2737/S1687: lay ToE-potential readiness matrix.

This is the next honest step after P2736 because the user asked for a lay
explanation of ToE potential.  The script keeps it proof-grade by converting the
question into a finite readiness matrix: which necessary ToE conditions are
already supported by content-grep evidence, and which remain blocked by current
guardrails.  It does not claim ToE closure.
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
OUT = GEN / "p2737_s1687_lay_toe_potential_readiness_matrix.json"
MD = GEN / "p2737_s1687_lay_toe_potential_readiness_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {
    "P2736_CONTENT_GREP_GATE": GEN / "p2736_s1686_post_p2735_content_grep_no_new_frontier_certificate.json",
}

CONTENT_PATTERNS = {
    "corrected_nadsoliton_ontology": r"nadsoliton.*primordial information|no separate informational layer|nadsoliton -> light -> matter -> emergent observer",
    "finite_computational_artifacts": r"finite calculation|finite audit|enumerates all|2\^12|2\^4|rank\(d0\)|dim H1|acceptance matrix",
    "selector_blocker": r"QW-2191.*open|QW-2191.*obstruction|no.*QW-2191 discharge|does not discharge `QW-2191`",
    "lambda_p2721_blocker": r"lambda/P2721|P2721 polarity|lambda remains unfixed|strict law fixing `lambda`|polarity remains unsourced",
    "bridge_blocker": r"bridge completion|legacy -> strict|role-transfer theorem|not.*bridge theorem|bridge closure",
    "ltotal_blocker": r"role-bearing `L_total`|cannot promote `L_total`|no new source term to add to `L_total`|L_total.*blocked",
    "toe_blocker": r"ToE closure|TOE closure|does not claim ToE|No ToE closure|ToE lanes remain blocked",
}

NEGATIVE_EXPORT_FLAGS = [
    "toe_closure_exported",
    "selector_closure_exported",
    "qw2191_discharged",
    "lambda_p2721_polarity_selected",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
]


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


def evidence_scan() -> dict[str, Any]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({
            "content_lane": name,
            "pattern": pattern,
            "hit_count": len(hits),
            "sample_hits": hits[:8],
        })
    by_name = {row["content_lane"]: row for row in rows}
    return {
        "content_pattern_count": len(CONTENT_PATTERNS),
        "rows": rows,
        "all_patterns_have_hits": all(row["hit_count"] > 0 for row in rows),
        "hit_counts": {name: by_name[name]["hit_count"] for name in CONTENT_PATTERNS},
    }


def readiness_matrix(scan: dict[str, Any]) -> dict[str, Any]:
    rows = [
        {
            "condition": "corrected_ontology",
            "lay_meaning": "The theory has a clear starting picture: one informational nadsoliton, not a hidden lower information layer.",
            "status": "supported",
            "met": scan["hit_counts"].get("corrected_nadsoliton_ontology", 0) > 0,
        },
        {
            "condition": "finite_computational_discipline",
            "lay_meaning": "Many claims are being tested by finite matrices/enumerations rather than only verbal analogy.",
            "status": "supported",
            "met": scan["hit_counts"].get("finite_computational_artifacts", 0) > 0,
        },
        {
            "condition": "strict_selector_or_orientation_source",
            "lay_meaning": "The theory still needs an internal rule choosing the missing sign/direction, not a human convention.",
            "status": "blocked_open_QW2191",
            "met": False,
        },
        {
            "condition": "lambda_P2721_polarity_source",
            "lay_meaning": "The recent chiral/time-arrow lane still needs one non-premise polarity choice.",
            "status": "blocked_unfixed_lambda_polarity",
            "met": False,
        },
        {
            "condition": "legacy_to_strict_bridge_completion",
            "lay_meaning": "The old intuitive kernel and the later strict kernel are not yet connected by a completed theorem.",
            "status": "blocked_bridge_incomplete",
            "met": False,
        },
        {
            "condition": "role_transfer_theorem",
            "lay_meaning": "Physical roles from the legacy story cannot be inherited until the bridge is finished and audited.",
            "status": "blocked_downstream_of_bridge",
            "met": False,
        },
        {
            "condition": "role_bearing_L_total",
            "lay_meaning": "A full action/Lagrangian carrying the physical roles is not yet licensed.",
            "status": "blocked_no_variational_source",
            "met": False,
        },
        {
            "condition": "toe_closure",
            "lay_meaning": "A true Theory of Everything claim would require all the above, so it is not currently closed.",
            "status": "not_claimed",
            "met": False,
        },
    ]
    met_count = sum(1 for row in rows if row["met"])
    total = len(rows)
    return {
        "rows": rows,
        "met_count": met_count,
        "total_conditions": total,
        "readiness_fraction": f"{met_count}/{total}",
        "closure_allowed": met_count == total,
        "lay_verdict": "High research potential, low closure readiness: the framework has a coherent ontology and strong finite-audit discipline, but the sign/selector source, lambda/P2721 polarity, bridge, role transfer, and role-bearing L_total remain open.",
    }


def acceptance_matrix(matrix: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "ontology_supported": any(row["condition"] == "corrected_ontology" and row["met"] for row in matrix["rows"]),
        "finite_audit_discipline_supported": any(row["condition"] == "finite_computational_discipline" and row["met"] for row in matrix["rows"]),
        "strict_selector_source_closed": any(row["condition"] == "strict_selector_or_orientation_source" and row["met"] for row in matrix["rows"]),
        "lambda_p2721_polarity_closed": any(row["condition"] == "lambda_P2721_polarity_source" and row["met"] for row in matrix["rows"]),
        "bridge_completion_closed": any(row["condition"] == "legacy_to_strict_bridge_completion" and row["met"] for row in matrix["rows"]),
        "role_bearing_ltotal_closed": any(row["condition"] == "role_bearing_L_total" and row["met"] for row in matrix["rows"]),
        "toe_closure_allowed": matrix["closure_allowed"],
    }
    missing = [key for key, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_toe_closure": facts["toe_closure_allowed"],
        "blocker": "The ToE-potential matrix supports research promise but blocks ToE closure because strict selector/polarity, bridge, role-transfer, and role-bearing L_total conditions remain unmet.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["toe_readiness_matrix"]
    lines = [
        "# P2737/S1687 lay ToE-potential readiness matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Readiness summary",
        f"- readiness_fraction={matrix['readiness_fraction']}",
        f"- closure_allowed={matrix['closure_allowed']}",
        f"- lay_verdict={matrix['lay_verdict']}",
        "",
        "## Lay explanation",
        payload["lay_explanation"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    matrix = readiness_matrix(scan)
    acceptance = acceptance_matrix(matrix)
    payload = {
        "status": "P2737_LAY_TOE_POTENTIAL_READINESS_MATRIX_NO_CLOSURE" if not acceptance["accepted_as_toe_closure"] else "P2737_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_question": "For a lay reader, what is the current ToE potential without overclaiming closure?",
        "content_evidence_scan": scan,
        "toe_readiness_matrix": matrix,
        "acceptance_matrix": acceptance,
        "lay_explanation": "For a non-specialist: the program has ToE potential because it tries to derive physics from one primordial informational-soliton structure and it keeps turning vague sign/phase questions into finite checks.  But it is not a ToE yet.  The missing piece is like a compass needle with two equally possible directions: until the theory itself chooses the direction/sign and connects the legacy and strict kernels with role-safe dynamics, the framework remains a promising research program rather than a completed final theory.",
        "decision": {
            "toe_closure_exported": False,
            "research_potential_preserved": True,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": matrix["lay_verdict"],
            "next_honest_step": "For ToE progress, do not make a closure claim.  Supply one new strict signed source law that fixes the selector/lambda-P2721 polarity and then rerun the bridge/role-transfer/L_total readiness matrix; without that new object, preserve the P2697-P2737 no-new-live-frontier and no-ToE-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2737/S1687 lay ToE-potential readiness matrix", "## P2737/S1687 lay ToE-potential readiness matrix\n\n`P2737/S1687` converts the lay question about ToE potential into a finite readiness matrix.  The matrix marks corrected ontology and finite-audit discipline as supported, but keeps strict selector/polarity source, `lambda/P2721` closure, legacy-to-strict bridge completion, role-transfer theorem, role-bearing `L_total`, and ToE closure unmet.  The lay verdict is therefore: high research potential, low closure readiness.  No `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2737/S1687 lay ToE-potential Ltotal guard", "## P2737/S1687 lay ToE-potential Ltotal guard\n\n`P2737/S1687` explicitly blocks promotion from research potential to role-bearing `L_total`: the readiness matrix still lacks a strict signed selector/polarity source, bridge completion, and role-transfer theorem.  Therefore no new Lagrangian term or ToE closure is added.\n")
    append_once(AGENTS, "Current lay ToE-potential readiness guardrail (P2737/S1687, 2026-06-14)", "## Current lay ToE-potential readiness guardrail (P2737/S1687, 2026-06-14)\n\n- P2737 translates the lay ToE-potential question into a finite readiness matrix: corrected ontology and finite-audit discipline are supported, but strict selector/polarity source, `lambda/P2721` closure, legacy-to-strict bridge completion, role-transfer theorem, role-bearing `L_total`, and ToE closure are unmet.\n- The honest lay verdict is high research potential but low closure readiness: the program is promising because it is coherent and computationally disciplined, but it is not a completed Theory of Everything.\n- Do not promote ToE, selector closure, bridge closure, role transfer, or `L_total` from lay-potential language.  A next admissible move must supply one new strict signed source law fixing the selector/`lambda-P2721` polarity, then rerun bridge/role-transfer/`L_total` readiness; otherwise preserve the P2697-P2737 no-new-live-frontier/no-ToE-closure certificate.\n")
    return payload


if __name__ == "__main__":
    main()
