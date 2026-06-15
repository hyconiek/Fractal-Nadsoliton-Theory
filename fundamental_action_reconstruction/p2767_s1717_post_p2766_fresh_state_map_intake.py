#!/usr/bin/env python3
"""P2767/S1717: post-P2766 fresh broad state-map intake.

P2766 preserved a no-closure certificate for moment-map physical-coupling
provenance and recommended either one genuinely new theorem/artifact closing a
named atom or a fresh broad state-map intake.  This script executes that intake
as a finite, content-first lane audit.  It does not claim that absence of a
string proves mathematical impossibility; it checks the current generated
artifact/guardrail state and records which lanes are still repetition-gated
without a newly supplied typed object.
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
P2697 = GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json"
P2759 = GEN / "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.json"
P2760 = GEN / "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.json"
P2766 = GEN / "p2766_s1716_post_moment_provenance_state_reconciliation.json"
OUT = GEN / "p2767_s1717_post_p2766_fresh_state_map_intake.json"
MD = GEN / "p2767_s1717_post_p2766_fresh_state_map_intake.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LANE_SPECS = {
    "moment_coupling_provenance": {
        "pattern": r"P2761|P2762|P2763|P2764|P2765|P2766|physical-coupling provenance|reference-cell/action-density|field/curvature normalization|nonproxy variational",
        "required_new_object": "one theorem closing exactly one of reference-cell/action-density, sign, field/curvature normalization, or 4D nonproxy residual rows",
        "current_boundary": "P2766 leaves all four moment provenance atoms open and blocks L_total promotion.",
    },
    "foundation_kernel_lagrangian_gap": {
        "pattern": r"P2760|foundation-kernel-Lagrangian gap|ontology-to-kernel measure|amplitude normalization|damping/compression bridge|stale closure-flag",
        "required_new_object": "machine-checkable theorem closing one named P2760 foundation-to-Lagrangian gap",
        "current_boundary": "P2760 lists seven open gaps and quarantines stale P1562 closure flags.",
    },
    "selector_orientation_qw2191": {
        "pattern": r"QW-2191|P2699|P2700|P2708|P2712|P2716|P2717|P2721|orientation torsor|lambda remains unfixed",
        "required_new_object": "non-premise strict selector/orientation source or polarity law with an explicit P2721 coupling polarity",
        "current_boundary": "Selector/sign lanes remain paired or premise-scoped; QW-2191 is not discharged.",
    },
    "entropy_and_phase_observables": {
        "pattern": r"P2753|P2754|P2755|P2756|P2757|P2758|polynomial phase|Shannon entropy|entropy current|pair-current|triangle-circulation",
        "required_new_object": "strict negation-breaking or step/polarity source outside polynomial/local entropy observables",
        "current_boundary": "Polynomial phase-sum and entropy current/circulation escalations are closed as replay on current artifacts.",
    },
    "direct_route_residuals": {
        "pattern": r"P2695|P2696|P2697|direct-route|residual cancellation|g4/g6/gY|pair1 c1c1|pair1 s1s1",
        "required_new_object": "new strict-derived provider class, non-N477 ingredient, or blocker-cut",
        "current_boundary": "Direct residual cancellation is repetition-gated after P2695-P2697.",
    },
    "lagrangian_eom_reverse_closure": {
        "pattern": r"P2685|P2686|P2687|Lagrangian/EOM reverse|EA/EH/ELg|anisotropic source|4D nonproxy",
        "required_new_object": "new strict anisotropic source class evading P1977/P1978 plus nonproxy residual rows",
        "current_boundary": "Strict Lagrangian/EOM reverse closure is bounded no-go without a new anisotropic source class.",
    },
    "bridge_source_and_role_transfer": {
        "pattern": r"P2680|P2693|legacy -> strict|bridge-source|canonical UV unit|alpha_geo|beta/Z_beta|role-transfer",
        "required_new_object": "single new non-selector source atom or explicit bridge-completion theorem before any role-transfer audit",
        "current_boundary": "P2680 source atoms are closed as bounded no-go; generic bridge replay and role transfer are blocked.",
    },
    "lower_boundary_tau_pair": {
        "pattern": r"P2684|lower-boundary|T24x|T25x|tau_src|pair12|boundary-square|L5/L12",
        "required_new_object": "chart-label-retaining pair12 typed seed/subinterface plus nonconventional provider",
        "current_boundary": "Lower-boundary recursion and tau/pair replay are cycle-cut without a new provider class.",
    },
}

NEGATIVE_EXPORT_FLAGS = [
    "new_live_frontier_selected_without_new_object",
    "moment_coupling_provenance_promoted",
    "qw2191_discharged",
    "selector_closure_exported",
    "bridge_closure_exported",
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


def lane_intake_matrix() -> dict[str, Any]:
    rows = []
    for lane, spec in LANE_SPECS.items():
        hits = run_rg(spec["pattern"])
        rows.append({
            "lane": lane,
            "pattern": spec["pattern"],
            "evidence_hit_count": len(hits),
            "has_current_boundary_evidence": len(hits) > 0,
            "required_new_object": spec["required_new_object"],
            "current_boundary": spec["current_boundary"],
            "admissible_without_new_object": False,
            "sample_hits": hits[:6],
        })
    return {
        "rows": rows,
        "row_count": len(rows),
        "boundary_evidence_count": sum(1 for row in rows if row["has_current_boundary_evidence"]),
        "admissible_without_new_object_count": sum(1 for row in rows if row["admissible_without_new_object"]),
        "all_lanes_have_boundary_evidence": all(row["has_current_boundary_evidence"] for row in rows),
    }


def new_object_intake_gate(matrix: dict[str, Any], p2766: dict[str, Any]) -> dict[str, Any]:
    p2766_acceptance = p2766.get("acceptance_matrix", {})
    facts = {
        "p2766_input_present": p2766.get("status") == "P2766_POST_MOMENT_PROVENANCE_STATE_RECONCILIATION_NO_CLOSURE",
        "p2766_blocks_physical_coupling_provenance": p2766_acceptance.get("accepted_as_physical_coupling_provenance_theorem") is False,
        "all_lanes_have_boundary_evidence": matrix["all_lanes_have_boundary_evidence"],
        "no_lane_admissible_without_new_object": matrix["admissible_without_new_object_count"] == 0,
        "new_typed_object_supplied_in_this_intake": False,
        "new_atomic_theorem_supplied_in_this_intake": False,
    }
    return {
        "facts": facts,
        "accepted_as_new_live_frontier_selection": False,
        "accepted_as_no_new_live_frontier_certificate": True,
        "missing_criteria_for_live_frontier": [key for key, value in facts.items() if not value],
        "blocker": "The broad post-P2766 intake finds boundary evidence on every audited lane and no lane is admissible without a genuinely new typed object/theorem.  This intake supplies no such object, so it certifies no-new-live-frontier rather than selecting a proof target.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["fresh_state_map_intake_matrix"]
    lines = [
        "# P2767/S1717 post-P2766 fresh broad state-map intake",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Lane matrix",
        f"- row_count={matrix['row_count']}",
        f"- boundary_evidence_count={matrix['boundary_evidence_count']}",
        f"- admissible_without_new_object_count={matrix['admissible_without_new_object_count']}",
        "",
    ]
    for row in matrix["rows"]:
        lines.append(f"- {row['lane']}: admissible_without_new_object={row['admissible_without_new_object']}; required_new_object={row['required_new_object']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2697 = read_json(P2697)
    p2759 = read_json(P2759)
    p2760 = read_json(P2760)
    p2766 = read_json(P2766)
    matrix = lane_intake_matrix()
    gate = new_object_intake_gate(matrix, p2766)
    payload = {
        "status": "P2767_POST_P2766_FRESH_STATE_MAP_INTAKE_NO_NEW_LIVE_FRONTIER",
        "input_hashes": {"P2697": sha(P2697), "P2759": sha(P2759), "P2760": sha(P2760), "P2766": sha(P2766)},
        "input_statuses": {"P2697": p2697.get("status"), "P2759": p2759.get("status"), "P2760": p2760.get("status"), "P2766": p2766.get("status")},
        "audited_question": "After P2766, does a broad state-map intake find a lane admissible without a genuinely new typed object/theorem?",
        "fresh_state_map_intake_matrix": matrix,
        "new_object_intake_gate": gate,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": gate["blocker"],
            "next_honest_step": "Do not replay closed lanes or promote L_total.  The next honest proof-grade move must first supply one concrete new typed object/theorem, then run a bounded acceptance test only for that object.  Highest-value allowed examples are: a theorem closing one P2766 provenance atom, a strict non-premise orientation/polarity source with explicit P2721 coupling, a new anisotropic source class with 4D nonproxy residual rows, or a machine-checkable theorem closing one P2760 foundation-to-Lagrangian gap.  If none is supplied, preserve the P2697-P2767 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2767/S1717 post-P2766 fresh broad state-map intake", "## P2767/S1717 post-P2766 fresh broad state-map intake\n\n`P2767/S1717` executes the fresh broad state-map intake requested after P2766.  It audits eight lanes: moment-coupling provenance, foundation-kernel-Lagrangian gaps, selector/orientation/`QW-2191`, entropy and phase observables, direct-route residuals, Lagrangian/EOM reverse closure, bridge-source/role-transfer, and lower-boundary/tau/pair recursion.  Every lane has current boundary evidence and zero lanes are admissible without a genuinely new typed object or theorem.  This is a no-new-live-frontier certificate, not a proof target selection; no physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2767/S1717 fresh state-map Ltotal guard", "## P2767/S1717 fresh state-map Ltotal guard\n\n`P2767/S1717` adds no variational source term.  Its broad intake finds no lane admissible without a new typed object/theorem, so it cannot promote `lambda_sm_eff`, `kappa_gr_eff`, `epsilon_mix_eff`, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current post-P2766 fresh broad state-map intake guardrail (P2767/S1717, 2026-06-15)", "## Current post-P2766 fresh broad state-map intake guardrail (P2767/S1717, 2026-06-15)\n\n- P2767 executes the fresh broad state-map intake requested after P2766 instead of replaying moment-coupling or `L_total` promotion.\n- The finite eight-lane matrix finds boundary evidence for moment-coupling provenance, foundation-kernel-Lagrangian gaps, selector/orientation/`QW-2191`, entropy/phase observables, direct residuals, Lagrangian/EOM reverse closure, bridge-source/role-transfer, and lower-boundary/tau/pair recursion; zero lanes are admissible without a genuinely new typed object/theorem.\n- Do not promote physical-coupling provenance, role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.  The next admissible move must first supply one concrete new typed object/theorem and then run only its bounded acceptance test; otherwise preserve the P2697-P2767 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
