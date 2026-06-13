#!/usr/bin/env python3
"""P2702/S1652: selector-on-the-circle lay mechanism and status packet.

Explains, with a finite Z12-circle toy calculation, what the selector was meant
to do: choose one representative/direction from symmetric possibilities.  It
separates the intuitive mechanism (a premise selector can pick +1 and cut the
circle) from the current proof status (no strict-sourced provider in P2701).
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.json"
MD = GEN / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2701": GEN / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.json",
    "P2700": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2699": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
    "P2698": GEN / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "new_selector_source_exported",
    "qw2191_discharged",
    "strict_selector_closure_exported",
    "convexity_breaking_promoted_to_theorem",
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
        "selector_status_chain": r"P2698|P2699|P2700|P2701|selector|direction|orientation",
        "circle_language": r"circle|Z12|Z_12|Aut\(Z12\)|generator|\+1|\-1|orbit",
        "premise_vs_strict": r"premise-based|strict-sourced|non-premise|QW-2191|no-new-live-frontier",
        "convexity_language": r"convex|wypuk|cut|branch|sector|arc",
    }
    return {"tool": "rg", "mode": "P2702 selector-on-circle lay mechanism and proof status", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def z12_circle_models() -> dict[str, Any]:
    units = [u for u in range(12) if math.gcd(u, 12) == 1]
    generator_orbit = sorted({(u * 1) % 12 for u in units})
    minus_one = 11
    models = [
        {
            "model": "no_selector_circle",
            "selected_direction": None,
            "available_unit_directions": units,
            "residual_symmetry_units": units,
            "single_direction_chosen": False,
            "plain_meaning": "Without a selector, the 12-node circle has several equivalent generator directions; no arrow is privileged.",
        },
        {
            "model": "premise_selector_plus_one",
            "selected_direction": 1,
            "available_unit_directions": units,
            "residual_symmetry_units": [u for u in units if (u * 1) % 12 == 1],
            "single_direction_chosen": True,
            "plain_meaning": "A premise selector can declare +1 as the arrow; this breaks the symmetry by assumption, like cutting a circular ambiguity at one marked arrow.",
        },
        {
            "model": "aut_invariant_internal_candidate",
            "selected_direction": None,
            "available_unit_directions": generator_orbit,
            "residual_symmetry_units": units,
            "single_direction_chosen": False,
            "plain_meaning": "If the rule must be Aut(Z12)-invariant, +1 and -1 stay in the same orbit, so the circle is not strictly cut to one arrow.",
        },
    ]
    return {
        "units": units,
        "generator_orbit_of_plus_one": generator_orbit,
        "minus_one_in_plus_one_orbit": minus_one in generator_orbit,
        "models": models,
        "premise_selector_breaks_symmetry_in_toy_model": models[1]["single_direction_chosen"] and len(models[1]["residual_symmetry_units"]) == 1,
        "aut_invariant_candidate_breaks_symmetry_in_toy_model": models[2]["single_direction_chosen"],
        "convexity_lay_note": "The 'convexity break on a circle' is best read as cutting a cyclic/flat degeneracy into one branch or sector.  It is intuitive branch selection, not a strict convexity theorem exported here.",
    }


def state_reads() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in INPUTS.items()}
    p2701 = loaded["P2701"]
    p2700 = loaded["P2700"]
    p2699 = loaded["P2699"]
    p2698 = loaded["P2698"]
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2698_real_direction_support_preserved": "real direction/orientation" in p2698.get("decision", {}).get("acknowledgement", ""),
        "p2699_aut_invariant_no_go": p2699.get("decision", {}).get("bounded_no_go_now") is True,
        "p2700_exhaustive_no_go": p2700.get("decision", {}).get("bounded_no_go_now") is True,
        "p2701_no_provider": p2701.get("decision", {}).get("bounded_no_go_now") is True and p2701.get("provider_inventory", {}).get("accepted_provider_count") == 0,
    }


def lay_explanation_rows(circle: dict[str, Any], reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "question": "What is the selector?",
            "plain_answer": "A selector is a rule that chooses one representative from many symmetric possibilities.  On a 12-node circle it is like choosing one arrow around the circle instead of treating all equivalent arrows as interchangeable.",
            "proof_status": "mechanism_real_as_concept",
        },
        {
            "question": "How does it break symmetry and pick one direction?",
            "plain_answer": "A premise selector can mark +1 as the arrow.  In the toy model this reduces the residual unit symmetry from [1,5,7,11] to [1], so one direction is chosen.",
            "computed_support": circle["models"][1],
            "proof_status": "works_if_the_selector_is_given_as_a_premise",
        },
        {
            "question": "What about convexity on the circle?",
            "plain_answer": circle["convexity_lay_note"],
            "proof_status": "lay_metaphor_or_branch_cut_not_a_new_strict_convexity_theorem",
        },
        {
            "question": "Why does current repo not claim final closure?",
            "plain_answer": "Because P2699/P2700 show Aut-invariant internal information does not choose +1 over -1, and P2701 finds no new strict-sourced provider.  The earlier direction support is real but premise-based.",
            "computed_support": reads,
            "proof_status": "no_false_promotion",
        },
    ]


def decision(reads: dict[str, Any], circle: dict[str, Any]) -> dict[str, Any]:
    no_false_pass = reads["p2698_real_direction_support_preserved"] and reads["p2699_aut_invariant_no_go"] and reads["p2700_exhaustive_no_go"] and reads["p2701_no_provider"]
    return {
        "decision": "P2702_SELECTOR_CIRCLE_LAY_MECHANISM_AND_STATUS_PACKET_NO_FALSE_PASS",
        "premise_selector_can_pick_one_direction_in_toy_model": circle["premise_selector_breaks_symmetry_in_toy_model"],
        "strict_aut_invariant_candidate_picks_one_direction": circle["aut_invariant_candidate_breaks_symmetry_in_toy_model"],
        "no_new_closure_exported": no_false_pass,
        "next_honest_step": "The next proof-grade move cannot be another explanation or Aut-invariant replay.  It must construct an actual new strict-sourced symmetry-breaking provider with non-premise provenance, or introduce a different new typed object outside closed lanes.  Without that, keep the P2697-P2702 no-new-live-frontier certificate.",
        "forbidden_promotions": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2702/S1652 selector-on-circle lay mechanism and status", "", f"Status: `{payload['status']}`", "", "## Lay explanation"]
    for row in payload["lay_explanation_rows"]:
        lines.append(f"- **{row['question']}** {row['plain_answer']}  Status: `{row['proof_status']}`")
    lines.extend(["", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    circle = z12_circle_models()
    rows = lay_explanation_rows(circle, reads)
    payload: dict[str, Any] = {
        "status": "P2702_SELECTOR_CIRCLE_LAY_MECHANISM_AND_STATUS_PACKET_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "z12_circle_models": circle,
        "lay_explanation_rows": rows,
        "decision": decision(reads, circle),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2702/S1652 selector circle lay mechanism status",
        "## P2702/S1652 selector circle lay mechanism status\n\n"
        "`P2702/S1652` explains the selector mechanism with a finite Z12-circle model.  A premise selector can mark `+1` and thereby choose one direction in the toy circle, but the Aut-invariant internal candidate cannot distinguish `+1` from `-1`; P2701 also found no new strict-sourced provider.  Thus the selector story is an intuitive branch/direction mechanism, not a new `QW-2191`, strict selector, `L_total`, or ToE closure export.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2702/S1652 selector explanation Ltotal guard",
        "## P2702/S1652 selector explanation Ltotal guard\n\n"
        "`P2702/S1652` is explanatory and finite-model diagnostic.  It does not promote the premise selector into a variational source, `L_total`, role transfer, bridge closure, or ToE closure.\n",
    )
    append_once(
        AGENTS,
        "Current selector circle explanation guardrail (P2702/S1652, 2026-06-13)",
        "## Current selector circle explanation guardrail (P2702/S1652, 2026-06-13)\n\n"
        "- P2702 explains for humans that a selector is a rule choosing one representative/direction on the Z12 circle.  A premise selector can choose `+1`, but current strict artifacts do not source that choice non-premise; Aut-invariant candidates still fail and P2701 finds no provider.\n"
        "- Treat the circle/convexity language as branch-cut intuition unless a strict convexity/selector theorem is newly exported.\n"
        "- The next admissible move must construct an actual new strict-sourced symmetry-breaking provider or introduce a different new typed object; otherwise preserve the P2697-P2702 no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    main()
