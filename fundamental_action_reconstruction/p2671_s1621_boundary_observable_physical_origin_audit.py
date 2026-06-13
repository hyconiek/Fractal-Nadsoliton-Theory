#!/usr/bin/env python3
"""P2671/S1621: audit physical-origin observables for P2670 Boolean variables.

The audit is deliberately conservative: it does not introduce another formal
Boolean lift.  It checks whether existing repository content already supplies a
bridge-completed boundary-dynamics observable that can define the variables used
by P2669/P2670 (`pair2`, `sector1`, and any auxiliary lift) without imported
coding conventions.
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
OUT = GEN / "p2671_s1621_boundary_observable_physical_origin_audit.json"
MD = GEN / "p2671_s1621_boundary_observable_physical_origin_audit.md"
P2670 = GEN / "p2670_s1620_higher_order_boolean_cross_invariant_lift_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "boundary_observable_physical_origin_exported",
    "pair_variable_physical_origin_exported",
    "boundary_sector_variable_physical_origin_exported",
    "auxiliary_lift_physical_origin_exported",
    "boolean_cross_invariant_source_exported",
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
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


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
    return {"count": len(lines), "samples": lines[:120]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "observable_origin_content": "boundary-dynamics observable|boundary dynamics|physical origin|not from labels|coding convention",
        "pair_variable_content": "pair2_positive|pair1/pair2|pair12|w_break|pair variable",
        "boundary_sector_content": "boundary-sector variable|boundary sector|sector1|sector swap|square holonomy|non-exact sector",
        "auxiliary_lift_content": "auxiliary bit|auxiliary lift|higher-order Boolean|Boolean lift|GF\\(2\\)",
        "selector_content": "selector|source-topology selector|QW-2191|strict selector source",
        "nonclosure_guard_content": "L_total|ToE closure|beta source|no false pass|nonexport",
    }
    return {
        "tool": "rg",
        "mode": "content-first boundary-observable physical-origin audit, not number/name-only search",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2670 = load_json(P2670)
    decision = p2670.get("closure_decision", {})
    return {
        "p2670_higher_order_candidates_exist": decision.get("mathematical_higher_order_candidates_exist") is True,
        "p2670_no_passing_source": decision.get("passing_source_count") == 0,
        "p2670_no_boundary_phase_bit_target": decision.get("boundary_phase_bit_target_exported_now") is False,
        "p2670_no_ltotal_reopening": decision.get("role_bearing_ltotal_now") is False,
    }


def observable_candidates() -> list[dict[str, Any]]:
    """Finite acceptance matrix for existing candidate observable lanes."""
    return [
        {
            "id": "source_topology_selector_plus_channel",
            "content_pattern": "source-topology selector|plus-channel|selector polarity|tau_src",
            "defines_pair_variable": True,
            "defines_boundary_sector_variable": False,
            "defines_auxiliary_lift": False,
            "bridge_completed_boundary_dynamics": False,
            "convention_free_coding": False,
            "blocker": "preLM/source-topology sign support does not descend to boundary-sector square holonomy coding",
        },
        {
            "id": "pair12_w_break_branch_split",
            "content_pattern": "w_break|pair1/pair2|pair2_positive|pair12",
            "defines_pair_variable": True,
            "defines_boundary_sector_variable": False,
            "defines_auxiliary_lift": False,
            "bridge_completed_boundary_dynamics": False,
            "convention_free_coding": False,
            "blocker": "branch split is real but sector-1 assignment remains convention-reversal sensitive",
        },
        {
            "id": "boundary_square_holonomy_carrier",
            "content_pattern": "square holonomy|non-exact sector|boundary sector|sector swap",
            "defines_pair_variable": False,
            "defines_boundary_sector_variable": True,
            "defines_auxiliary_lift": False,
            "bridge_completed_boundary_dynamics": False,
            "convention_free_coding": False,
            "blocker": "carrier exists but selector/sign source for sector 1 is still missing",
        },
        {
            "id": "orientation_cross_invariant_material",
            "content_pattern": "orientation datum|cross-invariant|boundary-cycle|sector swap",
            "defines_pair_variable": False,
            "defines_boundary_sector_variable": False,
            "defines_auxiliary_lift": False,
            "bridge_completed_boundary_dynamics": False,
            "convention_free_coding": False,
            "blocker": "orientation exports are lane-scoped/axis/convention data, not a boundary-cycle physical origin",
        },
        {
            "id": "self_recorded_endpoint_anchor",
            "content_pattern": "self-recorded endpoint anchor|endpoint-valued ledger|arrow action",
            "defines_pair_variable": False,
            "defines_boundary_sector_variable": False,
            "defines_auxiliary_lift": True,
            "bridge_completed_boundary_dynamics": False,
            "convention_free_coding": False,
            "blocker": "conditional selector candidate lacks strict endpoint-source convention and boundary-sector typing",
        },
        {
            "id": "apd_boundary_nullspace_measure_lane",
            "content_pattern": "APD boundary|nullspace|grid/measure|gamma selector|quadrature",
            "defines_pair_variable": False,
            "defines_boundary_sector_variable": False,
            "defines_auxiliary_lift": True,
            "bridge_completed_boundary_dynamics": False,
            "convention_free_coding": False,
            "blocker": "grid/measure-dependent selector cannot serve as convention-free auxiliary origin",
        },
    ]


def score_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scored = []
    for candidate in candidates:
        hits = rg_count(candidate["content_pattern"])
        passes = all(
            [
                hits["count"] > 0,
                candidate["bridge_completed_boundary_dynamics"],
                candidate["convention_free_coding"],
                candidate["defines_pair_variable"],
                candidate["defines_boundary_sector_variable"],
            ]
        )
        scored.append({**candidate, "content_hits": hits["count"], "content_samples": hits["samples"][:8], "passes_single_observable_origin_now": passes})
    return scored


def subset_witness(scored: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for size in range(1, len(scored) + 1):
        for subset in itertools.combinations(scored, size):
            covers_pair = any(item["defines_pair_variable"] for item in subset)
            covers_sector = any(item["defines_boundary_sector_variable"] for item in subset)
            covers_aux = any(item["defines_auxiliary_lift"] for item in subset)
            all_content_found = all(item["content_hits"] > 0 for item in subset)
            all_bridge_completed = all(item["bridge_completed_boundary_dynamics"] for item in subset)
            all_convention_free = all(item["convention_free_coding"] for item in subset)
            passes = all_content_found and covers_pair and covers_sector and all_bridge_completed and all_convention_free
            rows.append(
                {
                    "ids": [item["id"] for item in subset],
                    "covers_pair_variable": covers_pair,
                    "covers_boundary_sector_variable": covers_sector,
                    "covers_auxiliary_lift": covers_aux,
                    "all_content_found": all_content_found,
                    "all_bridge_completed_boundary_dynamics": all_bridge_completed,
                    "all_convention_free_coding": all_convention_free,
                    "passes_combined_origin_now": passes,
                }
            )
    passing = [row for row in rows if row["passes_combined_origin_now"]]
    near_misses = [
        row
        for row in rows
        if row["all_content_found"] and row["covers_pair_variable"] and row["covers_boundary_sector_variable"]
    ]
    return {
        "total_subsets_checked": len(rows),
        "passing_subset_count": len(passing),
        "near_miss_pair_and_sector_subset_count": len(near_misses),
        "sample_near_misses": near_misses[:12],
        "passing_subsets": passing,
    }


def closure_decision(subsets: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2671_BOUNDARY_OBSERVABLE_PHYSICAL_ORIGIN_AUDIT__NO_CURRENT_ORIGIN_SOURCE",
        "professorial_verdict": (
            "P2671 follows P2670 by auditing actual candidate observable lanes instead of adding another Boolean formal lift. "
            "The repo contains relevant pair, sector, auxiliary, and selector material, and finite subset search finds near-miss combinations covering pair and boundary sector. "
            "However every near-miss still fails bridge-completed boundary-dynamics provenance and convention-free coding. "
            "Thus no current observable defines the P2669/P2670 Boolean variables as physical origins."
        ),
        "next_honest_step": (
            "Pick exactly one near-miss observable bridge, preferably source_topology_selector_plus_channel + boundary_square_holonomy_carrier, "
            "and try to construct a typed descent theorem preserving sign into boundary square holonomy. Acceptance requires current content hits, "
            "bridge-completed boundary-dynamics provenance, convention-free coding, and sector-swap reversal forbidden internally. If that fails, record a no-go for observable-origin transfer into Boolean entropy-bit sourcing."
        ),
        "passing_subset_count": subsets["passing_subset_count"],
        "near_miss_pair_and_sector_subset_count": subsets["near_miss_pair_and_sector_subset_count"],
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2671/S1621 boundary-observable physical-origin audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Candidate observable matrix"])
    for row in payload["observable_candidate_matrix"]:
        lines.append(
            f"- `{row['id']}`: hits=`{row['content_hits']}`, pair=`{row['defines_pair_variable']}`, "
            f"sector=`{row['defines_boundary_sector_variable']}`, aux=`{row['defines_auxiliary_lift']}`, "
            f"bridge_completed=`{row['bridge_completed_boundary_dynamics']}`, convention_free=`{row['convention_free_coding']}`, "
            f"passes=`{row['passes_single_observable_origin_now']}`; blocker: {row['blocker']}"
        )
    subsets = payload["subset_origin_witness"]
    decision = payload["closure_decision"]
    lines.extend(
        [
            "",
            "## Subset witness",
            f"Total subsets checked: `{subsets['total_subsets_checked']}`.",
            f"Near-miss pair+sector subsets: `{subsets['near_miss_pair_and_sector_subset_count']}`.",
            f"Passing origin subsets: `{subsets['passing_subset_count']}`.",
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
    candidates = score_candidates(observable_candidates())
    subsets = subset_witness(candidates)
    decision = closure_decision(subsets)
    payload: dict[str, Any] = {
        "status": "P2671_BOUNDARY_OBSERVABLE_PHYSICAL_ORIGIN_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2670": sha256_file(P2670),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "observable_candidate_matrix": candidates,
        "subset_origin_witness": subsets,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2671/S1621 boundary-observable physical-origin guard",
        "## P2671/S1621 boundary-observable physical-origin guard\n\n"
        "`P2671/S1621` audits actual candidate observable lanes for the P2669/P2670 Boolean variables instead of adding another formal lift.  Existing pair, boundary-sector, auxiliary, and selector materials give near-miss subsets, but none has both bridge-completed boundary-dynamics provenance and convention-free coding.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2671/S1621 boundary-observable physical-origin Ltotal guard",
        "## P2671/S1621 boundary-observable physical-origin Ltotal guard\n\n"
        "`P2671/S1621` keeps `L_total` closed to observable-origin Boolean entropy terms.  A future term must first provide a typed descent from a bridge-completed boundary-dynamics observable to the pair/sector/auxiliary variables with convention-free coding and internal sector-swap rejection.\n",
    )
    return payload


if __name__ == "__main__":
    main()
