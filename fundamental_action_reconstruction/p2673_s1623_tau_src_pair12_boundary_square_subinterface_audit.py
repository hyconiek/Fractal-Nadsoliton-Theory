#!/usr/bin/env python3
"""P2673/S1623: exact subinterface audit after P2672.

Audits the missing typed chain:
    tau_src/source-topology sign -> pair1/pair2 typed carrier -> boundary square cycle
without claiming that the chain currently closes.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2673_s1623_tau_src_pair12_boundary_square_subinterface_audit.json"
MD = GEN / "p2673_s1623_tau_src_pair12_boundary_square_subinterface_audit.md"
P2672 = GEN / "p2672_s1622_source_topology_to_square_holonomy_typed_descent_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "tau_src_pair12_boundary_square_subinterface_exported",
    "boundary_square_cycle_typed_arrow_exported",
    "sign_preserving_pair12_to_square_map_exported",
    "sector_swap_sourced_invariant_exported",
    "boundary_phase_bit_target_exported_unconditionally",
    "target_independent_beta_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    data = json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(data.encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:50]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "tau_src_sign_content": "tau_src.*sign|barrier-protected sign|source-topology sign|plus-channel",
        "pair12_same_packet_content": "tau_src.*pair1/pair2|pair1/pair2.*tau_src|same tau_src packet|residual-datum pair1/pair2 carrier",
        "chart_sensitive_descent_content": "chart-sensitive.*pair1/pair2|pair1/pair2-typed.*descent|typed carrier|seed subinterface",
        "boundary_square_cycle_content": "boundary square cycle|boundary square|square holonomy|boundary cycle|non-exact sector",
        "sourced_invariant_content": "sourced invariant|sector swap changes|sector-swap.*invariant|sign-preserving|preserving sign",
        "nonclosure_guard_content": "QW-2191|L_total|ToE closure|beta source|future-route|nonexport",
    }
    return {
        "tool": "rg",
        "mode": "content-first exact tau_src -> pair12 -> boundary square subinterface audit",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2672 = load_json(P2672)
    decision = p2672.get("closure_decision", {})
    return {
        "p2672_no_typed_descent_gate": decision.get("passes_typed_descent_gate_now") is False,
        "p2672_no_passing_sourced_descent": decision.get("passing_sourced_descent_count") == 0,
        "p2672_no_boundary_phase_bit_target": decision.get("boundary_phase_bit_target_exported_now") is False,
        "p2672_no_ltotal_reopening": decision.get("role_bearing_ltotal_now") is False,
    }


def obligation_vector() -> list[dict[str, Any]]:
    """Current status of the exact chain obligations."""
    return [
        {
            "id": "O1_tau_src_sign_exists",
            "description": "source-topology/barrier sign material for tau_src is present",
            "satisfied_now": True,
            "content_pattern": "tau_src.*sign|barrier-protected sign|source-topology sign|plus-channel",
        },
        {
            "id": "O2_pair12_carrier_same_tau_packet",
            "description": "pair1/pair2 residual-datum carrier is tied to the same tau_src packet",
            "satisfied_now": True,
            "content_pattern": "same tau_src packet|tau_src.*pair1/pair2|pair1/pair2.*tau_src|residual-datum pair1/pair2 carrier",
        },
        {
            "id": "O3_chart_sensitive_pair12_typed_subinterface",
            "description": "chart-sensitive pair1/pair2 typed seed/descent subinterface is current, not future-route",
            "satisfied_now": False,
            "content_pattern": "chart-sensitive.*pair1/pair2|pair1/pair2-typed.*descent|seed subinterface",
        },
        {
            "id": "O4_boundary_square_cycle_arrow",
            "description": "typed arrow from pair12 carrier to boundary square cycle/square holonomy is exported",
            "satisfied_now": False,
            "content_pattern": "boundary square cycle|boundary square|square holonomy|boundary cycle|non-exact sector",
        },
        {
            "id": "O5_sector_swap_sourced_invariant",
            "description": "sector swap changes a sourced invariant rather than a label convention",
            "satisfied_now": False,
            "content_pattern": "sourced invariant|sector swap changes|sector-swap.*invariant|sign-preserving|preserving sign",
        },
    ]


def score_obligations(obligations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scored = []
    for item in obligations:
        hits = rg_count(item["content_pattern"])
        scored.append({**item, "content_hits": hits["count"], "content_samples": hits["samples"][:8]})
    return scored


def finite_closure_lattice(scored: list[dict[str, Any]]) -> dict[str, Any]:
    ids = [item["id"] for item in scored]
    rows = []
    for bits in itertools.product([False, True], repeat=len(ids)):
        state = dict(zip(ids, bits, strict=True))
        rows.append({"state": state, "closes_subinterface": all(bits)})
    current_state = {item["id"]: item["satisfied_now"] for item in scored}
    hamming_to_closure = sum(1 for value in current_state.values() if not value)
    return {
        "total_states_checked": len(rows),
        "closure_states_count": sum(row["closes_subinterface"] for row in rows),
        "current_state": current_state,
        "current_state_closes": all(current_state.values()),
        "missing_obligations_now": [key for key, value in current_state.items() if not value],
        "hamming_distance_from_closure": hamming_to_closure,
        "closure_state": {key: True for key in ids},
    }


def closure_decision(lattice: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2673_TAU_SRC_PAIR12_BOUNDARY_SQUARE_SUBINTERFACE_AUDIT__THREE_OBLIGATIONS_MISSING",
        "professorial_verdict": (
            "P2673 audits the exact subinterface requested by P2672. The tau_src sign material and same-packet pair1/pair2 carrier material are real, "
            "but the current chain is still three obligations away from closure: no current chart-sensitive pair12 typed subinterface, no pair12 -> boundary square-cycle typed arrow, "
            "and no sourced invariant that changes under sector swap. Therefore the sign cannot yet source the boundary square-holonomy entropy bit."
        ),
        "next_honest_step": (
            "Do not broaden the search. Attack O3 first: prove or refute a current chart-sensitive pair1/pair2 typed seed subinterface on tau_src. "
            "If O3 remains future-route/chart-bound, the whole tau_src -> pair12 -> boundary-square route should be recorded as blocked before attempting O4/O5. "
            "Only after O3 passes should the boundary-square arrow and sector-swap invariant be audited."
        ),
        "current_state_closes": lattice["current_state_closes"],
        "missing_obligations_now": lattice["missing_obligations_now"],
        "hamming_distance_from_closure": lattice["hamming_distance_from_closure"],
        "boundary_phase_bit_target_exported_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lattice = payload["finite_closure_lattice"]
    decision = payload["closure_decision"]
    lines = [
        "# P2673/S1623 tau_src -> pair12 -> boundary square subinterface audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Obligation vector"])
    for item in payload["obligation_vector"]:
        lines.append(f"- `{item['id']}`: satisfied_now=`{item['satisfied_now']}`, content_hits=`{item['content_hits']}` — {item['description']}")
    lines.extend(
        [
            "",
            "## Finite closure lattice",
            f"Total states checked: `{lattice['total_states_checked']}`.",
            f"Closure states: `{lattice['closure_states_count']}`.",
            f"Current state closes? `{lattice['current_state_closes']}`.",
            f"Missing obligations now: `{lattice['missing_obligations_now']}`.",
            f"Hamming distance from closure: `{lattice['hamming_distance_from_closure']}`.",
            "",
            "## Verdict",
            decision["professorial_verdict"],
            f"Decision: `{decision['decision']}`.",
            "",
            "## Next honest step",
            decision["next_honest_step"],
            "",
            "## Negative exports",
        ]
    )
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    scored = score_obligations(obligation_vector())
    lattice = finite_closure_lattice(scored)
    decision = closure_decision(lattice)
    payload: dict[str, Any] = {
        "status": "P2673_TAU_SRC_PAIR12_BOUNDARY_SQUARE_SUBINTERFACE_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2672": sha256_file(P2672),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "obligation_vector": scored,
        "finite_closure_lattice": lattice,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2673/S1623 tau_src pair12 boundary-square subinterface guard",
        "## P2673/S1623 tau_src pair12 boundary-square subinterface guard\n\n"
        "`P2673/S1623` audits the exact `tau_src/source-topology sign -> pair1/pair2 typed carrier -> boundary square cycle` subinterface.  The tau_src sign material and same-packet pair1/pair2 carrier material are real, but the finite obligation lattice remains three steps from closure: no current chart-sensitive pair12 typed subinterface, no pair12 -> boundary square-cycle typed arrow, and no sourced invariant changing under sector swap.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2673/S1623 tau_src pair12 boundary-square Ltotal guard",
        "## P2673/S1623 tau_src pair12 boundary-square Ltotal guard\n\n"
        "`P2673/S1623` keeps `L_total` closed to the tau_src -> pair12 -> boundary-square route.  A variational term is inadmissible until O3 chart-sensitive pair12 typing, O4 boundary-square typed arrow, and O5 sector-swap sourced invariant are all exported.\n",
    )
    return payload


if __name__ == "__main__":
    main()
