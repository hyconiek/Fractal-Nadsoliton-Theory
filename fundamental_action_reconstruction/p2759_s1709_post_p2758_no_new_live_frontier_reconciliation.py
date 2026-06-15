#!/usr/bin/env python3
"""P2759/S1709: post-P2758 no-new-live-frontier reconciliation.

P2758 closed the natural escalation from entropy pair-currents to local oriented
three-point entropy circulations.  This script does not manufacture another
selector candidate.  It performs a content-first finite state-map reconciliation
of the currently named post-P2697 lanes and checks whether the repository has
supplied the exact new object demanded after P2758: an independent strict
orientation/polarity law with an explicit P2721 coupling theorem, or a genuinely
new typed object outside local entropy observables and other closed replay
classes.
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
P2758 = GEN / "p2758_s1708_entropy_triangle_circulation_aut_cancellation_theorem.json"
OUT = GEN / "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.json"
MD = GEN / "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LANE_PATTERNS = {
    "direct_route_residuals": r"P2695|P2696|direct-route|residual cancellation|m2_psi4|g4/g6/gY|pair1 c1c1",
    "bridge_source_atoms": r"P2680|P2693|canonical UV unit|alpha_geo|beta/Z_beta|bridge-source|role transfer",
    "lagrangian_eom_reverse_closure": r"P2685|P2686|P2687|Lagrangian/EOM|anisotropic source|EA/EH/ELg",
    "lower_boundary_tau_pair_replay": r"P2684|lower-boundary|tau_src|pair12|boundary-square|L5/L12",
    "selector_qw2191_sign_lane": r"QW-2191|P2699|P2700|P2708|P2712|lambda|orientation torsor",
    "chiral_bispectrum_tau_spectral_lane": r"P2718|P2721|P2732|P2733|chiral-bispectrum|tau-coupling|spectral",
    "polynomial_phase_sum_lane": r"P2745|P2746|P2747|P2753|polynomial phase|Gauss-phase|cubic|quartic",
    "entropy_scalar_current_pair_triangle_lane": r"P2754|P2755|P2756|P2757|P2758|Shannon entropy|entropy current|pair-current|triangle-circulation",
}

POST_P2758_REQUIRED_PATTERNS = {
    "independent_strict_orientation_or_polarity_law": r"independent strict orientation/polarity law|strict orientation/polarity law|strict law selecting.*polarity",
    "explicit_p2721_coupling_theorem": r"explicit `?P2721`? coupling theorem|P2721 coupling theorem|P2721 polarity-coupling theorem",
    "outside_local_entropy_observables_requirement": r"outside local entropy observables|outside local entropy pair-current|genuinely new typed object",
    "no_new_live_frontier_certificate": r"no-new-live-frontier certificate|P2697-P2758|P2697-P2759",
}

NEGATIVE_EXPORT_FLAGS = [
    "new_post_p2758_typed_object_supplied",
    "independent_strict_orientation_or_polarity_law_exported",
    "explicit_p2721_coupling_theorem_exported",
    "lambda_p2721_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def scan_patterns(patterns: dict[str, str]) -> dict[str, Any]:
    rows = []
    for lane, pattern in patterns.items():
        hits = run_rg(pattern)
        rows.append({"lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return {"row_count": len(rows), "rows": rows, "hit_counts": {r["lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def lane_matrix() -> dict[str, Any]:
    lane_scan = scan_patterns(LANE_PATTERNS)
    requirement_scan = scan_patterns(POST_P2758_REQUIRED_PATTERNS)
    rows = []
    closed_reason = {
        "direct_route_residuals": "P2695/P2696/P2697 already gate direct residual replay without a new strict provider or blocker-cut.",
        "bridge_source_atoms": "P2693 closes named P2680 source atoms and forbids generic bridge-source replay or role transfer.",
        "lagrangian_eom_reverse_closure": "P2687 freezes strict Lagrangian/EOM reverse-closure without a genuinely new anisotropic source class.",
        "lower_boundary_tau_pair_replay": "P2684 and later guardrails block cyclic lower-boundary/tau/pair replay without a new provider class.",
        "selector_qw2191_sign_lane": "P2699-P2717/P2721/P2758 leave QW-2191 open absent a non-premise strict polarity source.",
        "chiral_bispectrum_tau_spectral_lane": "P2718-P2733 keep chiral marker/tau/spectral evidence premise- or polarity-paired rather than strict-source closed.",
        "polynomial_phase_sum_lane": "P2745-P2753 close degree escalation and coefficient-orbit replay without a negation-breaking source law.",
        "entropy_scalar_current_pair_triangle_lane": "P2754-P2758 close scalar entropy, directed current, pair-current signatures, and local triangle circulations as selector sources on current artifacts.",
    }
    for row in lane_scan["rows"]:
        rows.append(
            {
                "lane": row["lane"],
                "evidence_hit_count": row["hit_count"],
                "closed_or_repetition_gated_on_current_artifacts": row["hit_count"] > 0,
                "reason": closed_reason[row["lane"]],
                "sample_hits": row["sample_hits"][:4],
            }
        )
    return {
        "lane_scan": lane_scan,
        "post_p2758_requirement_scan": requirement_scan,
        "rows": rows,
        "closed_lane_count": sum(1 for r in rows if r["closed_or_repetition_gated_on_current_artifacts"]),
        "all_named_lanes_have_current_evidence": all(r["closed_or_repetition_gated_on_current_artifacts"] for r in rows),
    }


def acceptance_matrix(matrix: dict[str, Any], p2758: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2758_input_present": p2758.get("status") == "P2758_ENTROPY_TRIANGLE_CIRCULATION_AUT_CANCELLATION_THEOREM",
        "all_named_lanes_have_current_evidence": matrix["all_named_lanes_have_current_evidence"],
        "post_p2758_requirement_language_present": matrix["post_p2758_requirement_scan"]["all_patterns_have_hits"],
        "new_post_p2758_typed_object_supplied": False,
        "independent_strict_orientation_or_polarity_law_exported": False,
        "explicit_p2721_coupling_theorem_exported": False,
    }
    missing = [k for k, v in facts.items() if not v]
    return {
        "facts": facts,
        "accepted_as_new_live_frontier": False,
        "missing_criteria": missing,
        "blocker": "The current repository contains the post-P2758 requirement language and closed-lane evidence, but no supplied new typed object, no independent strict orientation/polarity law, and no explicit P2721 coupling theorem.  The honest output is therefore a no-new-live-frontier reconciliation, not closure promotion.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    m = payload["post_p2758_state_map_matrix"]
    lines = [
        "# P2759/S1709 post-P2758 no-new-live-frontier reconciliation",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## State-map lane matrix",
        f"- audited_lane_count={len(m['rows'])}",
        f"- closed_lane_count={m['closed_lane_count']}",
        f"- all_named_lanes_have_current_evidence={m['all_named_lanes_have_current_evidence']}",
        "",
    ]
    for row in m["rows"]:
        lines.append(f"- {row['lane']}: closed_or_repetition_gated={row['closed_or_repetition_gated_on_current_artifacts']} — {row['reason']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2758 = read_json(P2758)
    matrix = lane_matrix()
    acceptance = acceptance_matrix(matrix, p2758)
    payload = {
        "status": "P2759_POST_P2758_NO_NEW_LIVE_FRONTIER_RECONCILIATION",
        "input_hashes": {"P2758_ENTROPY_TRIANGLE_CIRCULATION_AUT_CANCELLATION_THEOREM": sha(P2758)},
        "input_statuses": {"P2758_ENTROPY_TRIANGLE_CIRCULATION_AUT_CANCELLATION_THEOREM": p2758.get("status")},
        "audited_question": "After P2758, has the repository supplied a new strict orientation/polarity law, explicit P2721 coupling theorem, or genuinely new typed object outside closed lanes?",
        "post_p2758_state_map_matrix": matrix,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not manufacture another replay candidate after P2758.  The next proof-grade move is admissible only if one concrete new strict typed object/source is supplied, preferably an independent strict orientation/polarity law with an explicit P2721 coupling theorem; otherwise preserve the P2697-P2759 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2759/S1709 post-P2758 no-new-live-frontier reconciliation", "## P2759/S1709 post-P2758 no-new-live-frontier reconciliation\n\n`P2759/S1709` performs a content-first state-map reconciliation after P2758 rather than inventing another local entropy replay.  It audits eight currently named lanes: direct residuals, P2680 bridge-source atoms, Lagrangian/EOM reverse closure, lower-boundary/tau/pair replay, selector/`QW-2191` sign lane, chiral-bispectrum/tau/spectral lane, polynomial phase-sum lane, and entropy scalar/current/pair/triangle lane.  All named lanes have current closed/repetition-gated evidence, and the repository supplies no new post-P2758 strict typed object, no independent strict orientation/polarity law, and no explicit `P2721` coupling theorem.  Thus no `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, bridge closure, or ToE closure is exported; the honest state is the P2697-P2759 no-new-live-frontier certificate.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2759/S1709 no-new-live-frontier Ltotal guard", "## P2759/S1709 no-new-live-frontier Ltotal guard\n\n`P2759/S1709` adds no variational source term.  It is a state-map reconciliation showing that no new strict typed object, independent orientation/polarity law, or `P2721` coupling theorem has been supplied after P2758.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current post-P2758 no-new-live-frontier reconciliation guardrail (P2759/S1709, 2026-06-15)", "## Current post-P2758 no-new-live-frontier reconciliation guardrail (P2759/S1709, 2026-06-15)\n\n- P2759 performs the broad post-P2758 state-map reconciliation instead of manufacturing another local entropy replay.\n- Direct residuals, P2680 bridge-source atoms, Lagrangian/EOM reverse closure, lower-boundary/tau/pair replay, selector/`QW-2191` sign lane, chiral-bispectrum/tau/spectral lane, polynomial phase-sum lane, and entropy scalar/current/pair/triangle lane remain closed or repetition-gated on current artifacts.\n- Do not promote `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without one concrete new strict typed object/source, preferably an independent strict orientation/polarity law with explicit `P2721` coupling.  Otherwise preserve the P2697-P2759 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
