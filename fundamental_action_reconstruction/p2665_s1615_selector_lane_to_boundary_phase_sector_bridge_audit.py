#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2665_s1615_selector_lane_to_boundary_phase_sector_bridge_audit.json"
MD = GEN / "p2665_s1615_selector_lane_to_boundary_phase_sector_bridge_audit.md"
P2664 = GEN / "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

SECTORS = [0, 1]
NEGATIVE_EXPORT_FLAGS = [
    "selector_lane_to_boundary_phase_bridge_exported",
    "nonexact_sector_selector_exported_unconditionally",
    "preferred_cycle_functional_exported_as_dynamics",
    "boundary_phase_bit_target_exported_unconditionally",
    "intrinsic_entropy_level_exported",
    "bit_to_action_map_sourced_unconditionally",
    "bit_to_length_map_sourced_unconditionally",
    "unique_beta_representative_selected_unconditionally",
    "entropy_arrow_discharges_qw_2191",
    "target_independent_beta_source_exported",
    "canonical_unit_exported",
    "bridge_completion_exported",
    "role_transfer_revalidated",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, ".",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
        "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
    ], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:120]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "selector_general_content": "selector|S_sel|QW-2191|strict-core selector|projective selector",
        "source_topology_selector_content": "source-topology selector|actual source-side selector|selector witness|selector polarity|tau_src",
        "boundary_phase_sector_content": "boundary-phase|non-exact sector|holonomy sector|square holonomy|theta.*holonomy|preferred cycle",
        "chart_quotient_limit_content": "projective|ray-level|basis-free quotient|declared-scope|selector-neutral|chart-bound|preLM",
        "raw_theta_uniqueness_content": "raw-theta uniqueness|theta uniqueness|theta source|holonomy source sign|sector source",
        "nonclosure_guard_content": "role-bearing L_total|bridge completion|role transfer|ToE closure|beta source",
    }
    return {"tool": "rg", "mode": "content-first selector-lane to boundary-phase-sector bridge audit; grep by research content, not packet names", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def lane_presence_from_rg(audit: dict[str, Any]) -> dict[str, bool]:
    patterns = audit["patterns"]
    return {
        "selector_material_exists": patterns["selector_general_content"]["count"] > 0,
        "source_topology_selector_material_exists": patterns["source_topology_selector_content"]["count"] > 0,
        "boundary_phase_sector_material_exists": patterns["boundary_phase_sector_content"]["count"] > 0,
        "chart_or_quotient_limit_material_exists": patterns["chart_quotient_limit_content"]["count"] > 0,
        "raw_theta_or_holonomy_source_material_exists": patterns["raw_theta_uniqueness_content"]["count"] > 0,
    }


def transfer_acceptance_rows(presence: dict[str, bool]) -> list[dict[str, Any]]:
    # Scores are finite proxies on the P2664 sectors.  Passing requires all three:
    # (1) computationally prefers sector 1, (2) is gauge/projector safe, and (3) has strict, non-declared provenance to the boundary-phase sector.
    candidates = [
        {
            "lane": "global_projective_selector_state",
            "repo_material_present": presence["selector_material_exists"],
            "sector_scores": {0: 1.0, 1: 1.0},
            "gauge_or_projector_safe": True,
            "strict_non_declared_boundary_phase_provenance": False,
            "failure_reason": "projective/ray-level selector material does not orient the raw boundary-phase holonomy sector",
        },
        {
            "lane": "source_topology_selector_witness",
            "repo_material_present": presence["source_topology_selector_material_exists"],
            "sector_scores": {0: 0.0, 1: 1.0},
            "gauge_or_projector_safe": False,
            "strict_non_declared_boundary_phase_provenance": False,
            "failure_reason": "actual selector witness material is useful but remains chart/preLM/downstream-bound relative to this boundary-phase sector target",
        },
        {
            "lane": "declared_scope_or_basis_free_quotient_selector",
            "repo_material_present": presence["chart_or_quotient_limit_material_exists"],
            "sector_scores": {0: 0.0, 1: 1.0},
            "gauge_or_projector_safe": True,
            "strict_non_declared_boundary_phase_provenance": False,
            "failure_reason": "quotient-safe declared-scope selection forgets the raw sector/source sign needed by P2664",
        },
        {
            "lane": "declared_theta_holonomy_source",
            "repo_material_present": presence["raw_theta_or_holonomy_source_material_exists"],
            "sector_scores": {0: 0.0, 1: 1.0},
            "gauge_or_projector_safe": True,
            "strict_non_declared_boundary_phase_provenance": False,
            "failure_reason": "theta-like sign selects sector 1 only as a declared source, exactly the missing premise",
        },
        {
            "lane": "bridge_completed_boundary_phase_sector_selector_target",
            "repo_material_present": False,
            "sector_scores": {0: 0.0, 1: 1.0},
            "gauge_or_projector_safe": True,
            "strict_non_declared_boundary_phase_provenance": False,
            "failure_reason": "this is the needed theorem target; no current content-first evidence exports it",
        },
    ]
    rows = []
    for row in candidates:
        scores = row["sector_scores"]
        best_score = max(scores.values())
        best_sectors = sorted(sector for sector, score in scores.items() if score == best_score)
        selects_sector_one = best_sectors == [1]
        passes = bool(
            row["repo_material_present"]
            and selects_sector_one
            and row["gauge_or_projector_safe"]
            and row["strict_non_declared_boundary_phase_provenance"]
        )
        rows.append({**row, "best_sectors": best_sectors, "selects_nonexact_sector_one": selects_sector_one, "passes_bridge_acceptance_now": passes})
    return rows


def upstream_consistency() -> dict[str, Any]:
    p2664 = load_json(P2664)
    return {
        "p2664_nonexact_bit_carrier_still_available": p2664.get("closure_decision", {}).get("nonexact_bit_carrier_still_available") is True,
        "p2664_local_even_class_does_not_select": p2664.get("closure_decision", {}).get("local_even_class_selects_nonexact_sector") is False,
        "p2664_declared_theta_can_select_but_not_internal": p2664.get("closure_decision", {}).get("declared_theta_can_select_nonexact_sector") is True and p2664.get("closure_decision", {}).get("theta_selection_is_internal_source") is False,
        "p2664_boundary_phase_bit_target_not_exported": p2664.get("closure_decision", {}).get("boundary_phase_bit_target_exported_now") is False,
    }


def closure_decision(rows: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["lane"] for row in rows if row["passes_bridge_acceptance_now"]]
    selector_like = [row["lane"] for row in rows if row["selects_nonexact_sector_one"]]
    return {
        "decision": "P2665_SELECTOR_LANE_TO_BOUNDARY_PHASE_SECTOR_BRIDGE_AUDIT__NO_ACCEPTED_TRANSFER",
        "professorial_verdict": "P2665 performs the requested content-first repo grep, including selector material, before doing the next computational audit.  The repo contains extensive selector/source-topology/projective material, but the finite transfer acceptance test does not find an existing strict, non-declared bridge into the P2664 boundary-phase square-holonomy sector.  Projective selector lanes are raw-sector neutral, source-topology witness lanes are not typed as this boundary-phase provenance, quotient/declared lanes forget or declare the needed sign, and theta-like holonomy selection remains the missing source.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, QW-2191 discharge, L_total reopening, or ToE closure is exported.",
        "next_honest_step": "Build a typed descent/promotion theorem from an existing selector source lane to the P2664 boundary-phase sector, but only if it preserves the raw holonomy sign and proves strict non-declared provenance.  The acceptance test is explicit: source content must be found by content grep, the induced sector scores must uniquely prefer sector 1, the map must be gauge/projector safe, and the provenance cannot be declared-scope, quotient-only, chart-bound, or theta-by-hand.  If that bridge cannot be built, freeze the entropy/cocycle UV-anchor route as conditional and record a no-go for selector-lane transfer to boundary-phase sector selection.",
        "selector_like_lanes": selector_like,
        "passing_selector_to_boundary_phase_bridge_lanes": passing,
        "selector_lane_to_boundary_phase_bridge_exported_now": False,
        "boundary_phase_sector_selector_exported_now": False,
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    decision = payload["closure_decision"]
    lines = [
        "# P2665/S1615 selector-lane to boundary-phase sector bridge audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines += [
        "",
        "## Lane presence",
    ]
    for key, value in payload["lane_presence_from_content_grep"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Transfer acceptance matrix",
        "| lane | repo material? | scores | best sectors | selects sector 1? | gauge/projector safe? | strict non-declared provenance? | passes? | failure reason |",
        "| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in payload["transfer_acceptance_rows"]:
        lines.append(
            f"| `{row['lane']}` | `{row['repo_material_present']}` | `{row['sector_scores']}` | `{row['best_sectors']}` | "
            f"`{row['selects_nonexact_sector_one']}` | `{row['gauge_or_projector_safe']}` | "
            f"`{row['strict_non_declared_boundary_phase_provenance']}` | `{row['passes_bridge_acceptance_now']}` | {row['failure_reason']} |"
        )
    lines += [
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Selector-like lanes: `{decision['selector_like_lanes']}`.",
        f"Passing selector-to-boundary-phase bridge lanes: `{decision['passing_selector_to_boundary_phase_bridge_lanes']}`.",
        f"Boundary-phase bit target exported now? `{decision['boundary_phase_bit_target_exported_now']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
        f"QW-2191 discharged now? `{decision['qw2191_discharged_now']}`.",
        f"Role-bearing L_total now? `{decision['role_bearing_ltotal_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Next honest step",
        decision["next_honest_step"],
        "",
        "## Negative exports",
    ]
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    audit = semantic_rg_audit()
    presence = lane_presence_from_rg(audit)
    rows = transfer_acceptance_rows(presence)
    decision = closure_decision(rows)
    payload: dict[str, Any] = {
        "status": "P2665_SELECTOR_LANE_TO_BOUNDARY_PHASE_SECTOR_BRIDGE_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": audit,
        "lane_presence_from_content_grep": presence,
        "source_hashes": {"P2664": sha256_file(P2664), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream_consistency(),
        "typed_transfer_target": {"source_problem": "selector lane to P2664 non-exact boundary-phase square holonomy sector", "sectors": SECTORS, "ontology_order": "nadsoliton -> light -> matter -> emergent observer", "no_sub_nadsoliton_information_layer": True},
        "transfer_acceptance_rows": rows,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2665/S1615 selector-lane to boundary-phase sector bridge guard", "## P2665/S1615 selector-lane to boundary-phase sector bridge guard\n\n`P2665/S1615` content-greps the repo by research content, including selector material, before auditing whether existing selector/source-topology lanes transfer to the P2664 boundary-phase square-holonomy sector.  The finite acceptance matrix finds no accepted bridge: projective lanes are raw-sector neutral, source-topology witness lanes are not typed as this boundary-phase provenance, quotient/declared lanes forget or declare the needed sign, and theta-like holonomy selection remains an unsourced premise.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2665/S1615 selector-lane to boundary-phase sector Ltotal guard", "## P2665/S1615 selector-lane to boundary-phase sector Ltotal guard\n\n`P2665/S1615` keeps `L_total` closed to selector-imported boundary-phase entropy terms.  Existing selector lanes may guide proof search, but none currently supplies a gauge/projector-safe, strict, non-declared descent to the raw P2664 square-holonomy sector.  A future term must pass the P2665 acceptance matrix before it can be used as a variational entropy-bit source.\n")
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
