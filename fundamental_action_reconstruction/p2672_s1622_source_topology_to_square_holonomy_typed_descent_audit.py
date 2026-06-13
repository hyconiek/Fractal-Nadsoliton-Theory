#!/usr/bin/env python3
"""P2672/S1622: typed-descent audit for source-topology sign -> square holonomy.

This is the concrete near-miss selected by P2671.  It checks whether the
source-topology plus/sign lane can be typed as a descent into the boundary
square-holonomy sector without importing a sector-label convention.
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
OUT = GEN / "p2672_s1622_source_topology_to_square_holonomy_typed_descent_audit.json"
MD = GEN / "p2672_s1622_source_topology_to_square_holonomy_typed_descent_audit.md"
P2671 = GEN / "p2671_s1621_boundary_observable_physical_origin_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "source_topology_to_square_holonomy_typed_descent_exported",
    "sign_preserving_boundary_sector_map_exported",
    "sector_swap_reversal_forbidden_internally",
    "boundary_phase_bit_target_exported_unconditionally",
    "intrinsic_entropy_level_exported",
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
            "rg",
            "-n",
            pattern,
            ".",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "*.json",
            "-g",
            "!fundamental_action_reconstruction/generated/**",
            "-g",
            "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "source_topology_sign_content": "source-topology selector|plus-channel|barrier-protected sign|tau_src|selector polarity",
        "square_holonomy_content": "square holonomy|non-exact sector|boundary sector|boundary-phase|sector swap",
        "typed_descent_content": "typed descent|pair1/pair2-typed|chart-sensitive|descent interface|typed carrier",
        "sign_preservation_content": "sign-preserving|preserving sign|sign flow|orientation preserving|holonomy sign",
        "reversal_guard_content": "sector-swap reversal|convention reversal|reverse convention|sector swap|coding convention",
        "nonclosure_guard_content": "QW-2191|L_total|ToE closure|beta source|nonexport|no false pass",
    }
    return {
        "tool": "rg",
        "mode": "content-first source-topology sign to boundary square-holonomy typed-descent audit",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2671 = load_json(P2671)
    decision = p2671.get("closure_decision", {})
    return {
        "p2671_has_near_miss_pair_sector_subsets": decision.get("near_miss_pair_and_sector_subset_count", 0) > 0,
        "p2671_no_passing_observable_origin": decision.get("passing_subset_count") == 0,
        "p2671_no_boundary_phase_bit_target": decision.get("boundary_phase_bit_target_exported_now") is False,
        "p2671_no_ltotal_reopening": decision.get("role_bearing_ltotal_now") is False,
    }


def typed_descent_gate() -> dict[str, Any]:
    source_hits = rg_count("source-topology selector|plus-channel|barrier-protected sign|tau_src|selector polarity")
    boundary_hits = rg_count("square holonomy|non-exact sector|boundary sector|boundary-phase|sector swap")
    typed_hits = rg_count("pair1/pair2-typed descent|chart-sensitive.*descent|typed carrier|descent interface")
    reversal_hits = rg_count("sector-swap reversal|convention reversal|reverse convention|sector swap|coding convention")
    criteria = {
        "source_topology_sign_material_present": source_hits["count"] > 0,
        "boundary_square_holonomy_material_present": boundary_hits["count"] > 0,
        "typed_descent_material_present": typed_hits["count"] > 0,
        "sector_swap_reversal_guard_material_present": reversal_hits["count"] > 0,
        "current_strict_typed_descent_exported": False,
        "current_sign_preservation_proof_exported": False,
        "current_sector_swap_reversal_forbidden": False,
        "bridge_completed_boundary_dynamics_provenance": False,
    }
    criteria["passes_typed_descent_gate_now"] = all(criteria.values())
    return {
        "criteria": criteria,
        "hit_counts": {
            "source_topology_sign": source_hits["count"],
            "boundary_square_holonomy": boundary_hits["count"],
            "typed_descent": typed_hits["count"],
            "reversal_guard": reversal_hits["count"],
        },
    }


def finite_map_witness(gate: dict[str, Any]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    signs = {"negative": 0, "positive": 1}
    for positive_sector in (0, 1):
        mapping = {"positive": positive_sector, "negative": 1 - positive_sector}
        for swap in (False, True):
            effective = {sign: (sector ^ int(swap)) for sign, sector in mapping.items()}
            positive_maps_to_sector1 = effective["positive"] == 1
            sign_distinguishes_sector = effective["positive"] != effective["negative"]
            convention_reversal_available = not gate["criteria"]["current_sector_swap_reversal_forbidden"]
            passes = all(
                [
                    gate["criteria"]["source_topology_sign_material_present"],
                    gate["criteria"]["boundary_square_holonomy_material_present"],
                    gate["criteria"]["current_strict_typed_descent_exported"],
                    gate["criteria"]["current_sign_preservation_proof_exported"],
                    not convention_reversal_available,
                    positive_maps_to_sector1,
                    sign_distinguishes_sector,
                ]
            )
            rows.append(
                {
                    "declared_positive_sector_convention": positive_sector,
                    "sector_swap_applied": swap,
                    "effective_sign_to_sector_map": effective,
                    "positive_maps_to_sector1_after_swap": positive_maps_to_sector1,
                    "sign_distinguishes_sector": sign_distinguishes_sector,
                    "convention_reversal_available": convention_reversal_available,
                    "passes_as_sourced_descent_now": passes,
                }
            )
    return {
        "statement": (
            "The finite witness enumerates the two sign-to-sector conventions and their sector-swap images. "
            "A convention can make positive sign land in sector 1, but the swapped convention is equally available unless an internal typed descent forbids it."
        ),
        "source_sign_encoding": signs,
        "rows": rows,
        "total_maps_checked": len(rows),
        "maps_with_positive_to_sector1": sum(row["positive_maps_to_sector1_after_swap"] for row in rows),
        "passing_sourced_descent_count": sum(row["passes_as_sourced_descent_now"] for row in rows),
    }


def closure_decision(gate: dict[str, Any], witness: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2672_SOURCE_TOPOLOGY_TO_SQUARE_HOLONOMY_TYPED_DESCENT_AUDIT__NEAR_MISS_NO_SOURCE",
        "professorial_verdict": (
            "P2672 audits the strongest P2671 near-miss directly: source-topology sign/plus-channel material plus boundary square-holonomy carrier material. "
            "Both content classes are present, and finite sign-to-sector maps can place the positive sign in sector 1. "
            "However the current repo still lacks a strict typed descent, a sign-preservation proof, bridge-completed boundary-dynamics provenance, and an internal ban on sector-swap reversal. "
            "Therefore this is a real near-miss, not a boundary-phase entropy-bit source theorem."
        ),
        "next_honest_step": (
            "Audit the exact missing subinterface: tau_src/source-topology sign -> pair1/pair2-typed chart-sensitive carrier -> boundary square cycle. "
            "The next proof must export a concrete typed arrow and prove that sector swap changes a sourced invariant, not merely a label. "
            "If the arrow remains future-route or chart-bound, record a no-go for source-topology sign descent into square-holonomy entropy-bit sourcing."
        ),
        "passes_typed_descent_gate_now": gate["criteria"]["passes_typed_descent_gate_now"],
        "passing_sourced_descent_count": witness["passing_sourced_descent_count"],
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    gate = payload["typed_descent_gate"]
    witness = payload["finite_map_witness"]
    decision = payload["closure_decision"]
    lines = [
        "# P2672/S1622 source-topology sign to square-holonomy typed-descent audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Typed descent gate"])
    for name, value in gate["criteria"].items():
        lines.append(f"- `{name}`: `{value}`")
    lines.extend(
        [
            "",
            "## Finite map witness",
            witness["statement"],
            f"Total maps checked: `{witness['total_maps_checked']}`.",
            f"Maps with positive -> sector 1: `{witness['maps_with_positive_to_sector1']}`.",
            f"Passing sourced descents: `{witness['passing_sourced_descent_count']}`.",
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
    gate = typed_descent_gate()
    witness = finite_map_witness(gate)
    decision = closure_decision(gate, witness)
    payload: dict[str, Any] = {
        "status": "P2672_SOURCE_TOPOLOGY_TO_SQUARE_HOLONOMY_TYPED_DESCENT_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2671": sha256_file(P2671),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_descent_gate": gate,
        "finite_map_witness": witness,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2672/S1622 source-topology to square-holonomy typed-descent guard",
        "## P2672/S1622 source-topology to square-holonomy typed-descent guard\n\n"
        "`P2672/S1622` audits the strongest P2671 near-miss: source-topology sign/plus-channel material mapped to boundary square-holonomy sector.  Content exists on both sides and finite sign-to-sector maps can place positive sign in sector `1`, but the current repo lacks a strict typed descent, sign-preservation proof, bridge-completed boundary-dynamics provenance, and internal sector-swap-reversal ban.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2672/S1622 source-topology to square-holonomy Ltotal guard",
        "## P2672/S1622 source-topology to square-holonomy Ltotal guard\n\n"
        "`P2672/S1622` keeps `L_total` closed to source-topology-sign-derived square-holonomy entropy terms.  A variational coefficient still requires a concrete typed arrow from `tau_src` sign data through a pair1/pair2 chart-sensitive carrier into the boundary square cycle, with sign preservation and internal sector-swap rejection.\n",
    )
    return payload


if __name__ == "__main__":
    main()
