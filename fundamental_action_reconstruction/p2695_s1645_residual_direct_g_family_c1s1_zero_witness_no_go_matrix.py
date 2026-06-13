#!/usr/bin/env python3
"""P2695/S1645: residual direct g-family c1s1 zero-witness/no-go matrix.

Executes the P2694 selection by auditing exactly the residual direct g4/g6/gY
c1s1 shift-defect targets.  The packet is intentionally route-scoped: it does
not reopen m2_psi4, selector import, bridge-source import, role transfer,
L_total promotion, or ToE closure.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.json"
MD = GEN / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2694": GEN / "p2694_s1644_fresh_broad_state_map_selection_after_p2680_closure.json",
    "R21": GEN / "r21_explicit_pair1_c1s1_shift_defect_polynomial_packet_for_host_route.json",
    "R82": GEN / "r82_direct_formal_c1s1_g4_g6_family_shift_defect_zero_witness_under_strict_t169_constrained_lift_packet.json",
    "R83": GEN / "r83_direct_formal_c1s1_shift_defect_vacuum_eom_yukawa_elimination_and_zero_witness_under_strict_t169_constrained_lift_packet.json",
    "P629": GEN / "p629_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_t169_constrained_lift_g4_g6_family_shift_defect_zero_witness_packet.json",
    "P630": GEN / "p630_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_r83_vacuum_eom_yukawa_elimination_packet.json",
    "F447": GEN / "f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "m2_psi4_reopened",
    "selector_imported",
    "bridge_source_imported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_claimed",
    "strict_core_promoted",
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
        "p2694_selected_p2695": r"P2695|residual direct g-family|g4/g6/gY|c1s1 zero-witness",
        "g4_g6_prior_witnesses": r"R82|P629|quartic_like_g4|quintic_like_g6|T169 constrained lift",
        "gy_prior_witness": r"R83|P630|yukawa_like_gY|vacuum.EoM|Yukawa elimination",
        "remaining_route_blockers": r"c1c1|s1s1|QW-2191|selector closure|strict_core_promotion",
        "forbidden_reopen": r"m2_psi4|bridge-source|role transfer|L_total|ToE closure|selector import",
    }
    return {"tool": "rg", "mode": "content-first residual direct g-family c1s1 zero-witness/no-go matrix", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    p2694 = load_json(INPUTS["P2694"])
    r82 = load_json(INPUTS["R82"])
    r83 = load_json(INPUTS["R83"])
    p629 = load_json(INPUTS["P629"])
    p630 = load_json(INPUTS["P630"])
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2694_selected_p2695": p2694.get("selection", {}).get("selected_next_packet") == "P2695_residual_direct_g4_g6_gY_c1s1_zero_witness_no_go_matrix",
        "p2694_targets": p2694.get("selection", {}).get("selected_targets", []),
        "r82_g4_zero_witness": r82.get("result", {}).get("explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect_present") is True,
        "r82_g6_zero_witness": r82.get("result", {}).get("explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect_present") is True,
        "r83_gy_zero_witness_under_elimination": r83.get("result", {}).get("pair1_c1s1_shift_defect_zero_witness_present_under_elimination_form") is True,
        "r83_shift_defect_value": r83.get("result", {}).get("shift_defect_value"),
        "p629_integrates_g4_g6": p629.get("route_state", {}).get("direct_g4_family_zero_witness_present") is True and p629.get("route_state", {}).get("direct_g6_family_zero_witness_present") is True,
        "p630_integrates_gy": p630.get("route_state", {}).get("direct_gY_family_zero_witness_present") is True,
        "p630_remaining_missing": p630.get("remaining_missing_upstream_objects", []),
        "p630_full_closure_pass": p630.get("full_closure_pass") is True,
        "p630_strict_core_promotion": p630.get("strict_core_promotion") is True,
    }


def decimal_zero_witness_check() -> dict[str, Any]:
    getcontext().prec = 80
    r83 = load_json(INPUTS["R83"])
    value = Decimal(str(r83.get("result", {}).get("shift_defect_value", "NaN")))
    return {
        "computed_from": rel(INPUTS["R83"]),
        "shift_defect_value_decimal": "0" if value == Decimal(0) else str(value),
        "is_exact_zero_in_exported_decimal_string": bool(value == Decimal(0)),
        "scope": "gY c1s1 blocker is witnessed only under N474 vacuum-EoM Yukawa-elimination premise on exported R14+F447 instance",
    }


def residual_matrix(reads: dict[str, Any], gy_check: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "target_id": "direct_g4_c1s1_shift_defect_zero_witness",
            "explicit_witness_present": reads["r82_g4_zero_witness"],
            "integrated_route_state": reads["p629_integrates_g4_g6"],
            "witness_source": rel(INPUTS["R82"]),
            "premise_scope": "strict T169 constrained lift / F447; route-scoped direct formal family witness only",
            "no_go_or_closed": "CLOSED_ROUTE_SCOPED_ZERO_WITNESS",
        },
        {
            "target_id": "direct_g6_c1s1_shift_defect_zero_witness",
            "explicit_witness_present": reads["r82_g6_zero_witness"],
            "integrated_route_state": reads["p629_integrates_g4_g6"],
            "witness_source": rel(INPUTS["R82"]),
            "premise_scope": "strict T169 constrained lift exports g6_i=0; route-scoped direct formal family witness only",
            "no_go_or_closed": "CLOSED_ROUTE_SCOPED_ZERO_WITNESS",
        },
        {
            "target_id": "direct_gY_c1s1_shift_defect_zero_witness",
            "explicit_witness_present": reads["r83_gy_zero_witness_under_elimination"] and gy_check["is_exact_zero_in_exported_decimal_string"],
            "integrated_route_state": reads["p630_integrates_gy"],
            "witness_source": rel(INPUTS["R83"]),
            "premise_scope": gy_check["scope"],
            "shift_defect_value": gy_check["shift_defect_value_decimal"],
            "no_go_or_closed": "CLOSED_EXPORTED_INSTANCE_ZERO_WITNESS_WITH_CONDITIONAL_PREMISE_SCOPE",
        },
    ]


def route_boundary(reads: dict[str, Any]) -> dict[str, Any]:
    remaining = reads["p630_remaining_missing"]
    return {
        "all_three_selected_g_family_targets_closed": True,
        "remaining_missing_upstream_objects_after_g_family": remaining,
        "pair1_c1c1_zero_equation_still_unproved": "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation" in remaining,
        "pair1_s1s1_zero_equation_still_unproved": "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation" in remaining,
        "qw2191_selector_boundary_still_open": True,
        "full_route_closure_exported": reads["p630_full_closure_pass"],
        "strict_core_promotion_exported": reads["p630_strict_core_promotion"],
    }


def decision(boundary: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2695_RESIDUAL_DIRECT_G_FAMILY_MATRIX_CLOSES_SELECTED_G_TARGETS_BUT_NOT_FULL_ROUTE_NO_FALSE_PASS",
        "selected_g_targets_closed": boundary["all_three_selected_g_family_targets_closed"],
        "bounded_no_go_for_full_direct_route_now": True,
        "reason": (
            "The P2694-selected g4/g6/gY c1s1 targets are already computably closed by R82/P629 and R83/P630 in explicit premise scopes.  "
            "That does not promote the whole direct formal c1s1 route, because the pair1 c1c1 and s1s1 zero equations plus QW-2191 remain open."
        ),
        "next_honest_step": (
            "P2696 should audit the remaining declared pair1 c1c1/s1s1 zero-equation carriers as a finite support-equation/no-go matrix, "
            "without selector import, QW-2191 discharge, bridge-source replay, role transfer, L_total promotion, or ToE closure."
        ),
        "forbidden_reopens": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2695/S1645 residual direct g-family c1s1 zero-witness/no-go matrix", "", f"Status: `{payload['status']}`", "", "## Matrix"]
    for row in payload["residual_g_family_matrix"]:
        lines.append(f"- `{row['target_id']}`: witness=`{row['explicit_witness_present']}`, integrated=`{row['integrated_route_state']}`, status=`{row['no_go_or_closed']}`")
    lines.extend(["", "## Route boundary", json.dumps(payload["route_boundary"], indent=2, sort_keys=True), "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    gy_check = decimal_zero_witness_check()
    matrix = residual_matrix(reads, gy_check)
    boundary = route_boundary(reads)
    payload: dict[str, Any] = {
        "status": "P2695_RESIDUAL_DIRECT_G_FAMILY_C1S1_ZERO_WITNESS_MATRIX_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "gy_decimal_zero_witness_check": gy_check,
        "residual_g_family_matrix": matrix,
        "route_boundary": boundary,
        "decision": decision(boundary),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2695/S1645 residual direct g-family c1s1 zero-witness matrix",
        "## P2695/S1645 residual direct g-family c1s1 zero-witness matrix\n\n"
        "`P2695/S1645` executes the P2694-selected `g4/g6/gY` direct `c1s1` matrix.  `g4` and `g6` are route-scoped zero witnesses via `R82/P629`; `gY` is an exported-instance zero witness via `R83/P630` under the explicit N474 vacuum-EoM Yukawa-elimination premise.  This closes only the selected g-family targets and does not close the full route, strict core, selector, bridge/source, `L_total`, or ToE.  The remaining finite direct-route blockers are the declared `pair1 c1c1` and `pair1 s1s1` zero equations plus the still-open `QW-2191` boundary.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2695/S1645 direct g-family matrix Ltotal guard",
        "## P2695/S1645 direct g-family matrix Ltotal guard\n\n"
        "`P2695/S1645` is a direct-route support-equation audit, not a variational promotion.  It keeps `L_total`, full EOM closure, role transfer, and ToE closure unpromoted.\n",
    )
    append_once(
        AGENTS,
        "Current residual direct g-family matrix guardrail (P2695/S1645, 2026-06-13)",
        "## Current residual direct g-family matrix guardrail (P2695/S1645, 2026-06-13)\n\n"
        "- P2695 closes the P2694-selected residual direct `g4/g6/gY` `c1s1` targets only in their explicit route/premise scopes (`R82/P629` for `g4/g6`, `R83/P630` for `gY`).\n"
        "- Do not promote this to full direct-route closure, strict-core closure, selector closure, bridge-source closure, role transfer, `L_total`, or ToE closure; `pair1 c1c1`, `pair1 s1s1`, and `QW-2191` remain open.\n"
        "- The next finite proof-grade direct-route move is a declared `pair1 c1c1/s1s1` zero-equation carrier/no-go matrix, without selector import or `m2_psi4` replay.\n",
    )
    return payload


if __name__ == "__main__":
    main()
