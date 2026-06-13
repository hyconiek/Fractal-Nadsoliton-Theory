#!/usr/bin/env python3
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
OUT = GEN / "p2666_s1616_pair12_selector_witness_to_boundary_phase_sector_descent_audit.json"
MD = GEN / "p2666_s1616_pair12_selector_witness_to_boundary_phase_sector_descent_audit.md"
P2665 = GEN / "p2665_s1615_selector_lane_to_boundary_phase_sector_bridge_audit.json"
P731 = GEN / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
P741 = GEN / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe_summary.json"
P764 = GEN / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"
F947 = GEN / "f947_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_packet_summary.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

PAIR_BRANCHES = ["pair1", "pair2"]
BOUNDARY_SECTORS = [0, 1]
NEGATIVE_EXPORT_FLAGS = [
    "pair12_to_boundary_phase_sector_descent_exported",
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
        "pair12_witness_split_content": "pair1/pair2|pair12|w_break|witness split|orbit-direction",
        "typed_descent_promotion_content": "typed descent|promotion bridge|descent-bridge|source-side branch selection|chart-sensitive",
        "selector_source_content": "actual source-topology selector|selector witness|selector polarity|tau_src|S_sel",
        "boundary_phase_sector_content": "boundary-phase|square holonomy|non-exact sector|holonomy sector|preferred cycle",
        "nonclosure_guard_content": "nonexport|future-route|declared-scope|QW-2191|role-bearing L_total|ToE closure|beta source",
    }
    return {"tool": "rg", "mode": "content-first pair12 selector-witness to boundary-phase sector descent audit", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def prior_summary_audit() -> dict[str, Any]:
    p731 = load_json(P731)
    p741 = load_json(P741)
    p764 = load_json(P764)
    f947 = load_json(F947)
    p2665 = load_json(P2665)
    return {
        "p731_w_break_separates_pair12": p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches") is True,
        "p731_pair2_positive": p731.get("pair2_w_break_branch_score_sign") == "positive",
        "p731_promotion_bridge_exported": p731.get("t185_target_exported_on_current_repo_state") is True,
        "p741_actual_source_witness_exported": p741.get("current_actual_source_topology_selector_witness_exported") is True,
        "p741_witness_remains_prelm_not_pair12_typed": p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed") is True,
        "p741_promotion_bridge_exported": p741.get("t195_target_exported_on_current_repo_state") is True,
        "p764_typed_descent_target_exported": p764.get("t218_target_exported_on_current_repo_state") is True,
        "p764_typed_descent_target_future_route_only": p764.get("current_t218_target_is_future_route_only") is True,
        "f947_target_packet_no_false_pass": f947.get("no_false_pass") is True,
        "p2665_no_accepted_selector_bridge": p2665.get("closure_decision", {}).get("passing_selector_to_boundary_phase_bridge_lanes") == [],
    }


def descent_mapping_witness(summary: dict[str, Any]) -> dict[str, Any]:
    # Existing strict material gives a pair12 witness split (pair2 positive), but the P2664 target is sector {0,1}.
    # The only possible bijections are enumerated.  Each bijection is mathematically possible; none is exported as a typed descent map.
    rows = []
    for sector_perm in itertools.permutations(BOUNDARY_SECTORS):
        mapping = dict(zip(PAIR_BRANCHES, sector_perm, strict=True))
        positive_branch = "pair2"
        selected_sector = mapping[positive_branch]
        selects_nonexact_sector = selected_sector == 1
        preserves_raw_holonomy_sign = selects_nonexact_sector
        pair12_split_exported = summary["p731_w_break_separates_pair12"] and summary["p731_pair2_positive"]
        pair12_to_boundary_typed_map_exported = False
        current_source_witness_pair12_typed = not summary["p741_witness_remains_prelm_not_pair12_typed"]
        future_target_counts_as_current_source = summary["p764_typed_descent_target_exported"] and not summary["p764_typed_descent_target_future_route_only"]
        passes = bool(
            pair12_split_exported
            and current_source_witness_pair12_typed
            and pair12_to_boundary_typed_map_exported
            and future_target_counts_as_current_source
            and preserves_raw_holonomy_sign
        )
        rows.append({
            "candidate_mapping": mapping,
            "positive_pair12_branch": positive_branch,
            "selected_boundary_sector": selected_sector,
            "selects_nonexact_boundary_sector_one": selects_nonexact_sector,
            "preserves_raw_holonomy_sign": preserves_raw_holonomy_sign,
            "pair12_witness_split_exported": pair12_split_exported,
            "current_source_witness_pair12_typed": current_source_witness_pair12_typed,
            "pair12_to_boundary_typed_map_exported": pair12_to_boundary_typed_map_exported,
            "future_target_counts_as_current_source": future_target_counts_as_current_source,
            "passes_descent_acceptance_now": passes,
        })
    return {
        "statement": "P2666 enumerates the only two bijective descents from the existing pair1/pair2 witness split to the P2664 boundary sectors {0,1}.  One convention maps the positive pair2 witness branch to the non-exact boundary sector 1, but the opposite convention is equally mathematical.  Current summaries export the pair12 split and future typed-descent targets, not a current strict, non-declared pair12 -> boundary-phase sector map.",
        "rows": rows,
        "mathematically_possible_sector_one_mappings": [row for row in rows if row["selects_nonexact_boundary_sector_one"]],
        "passing_descent_mappings": [row for row in rows if row["passes_descent_acceptance_now"]],
        "orientation_convention_degeneracy_count": len(rows),
        "raw_holonomy_sign_preserved_by_at_least_one_mapping": any(row["preserves_raw_holonomy_sign"] for row in rows),
        "strict_current_descent_exported_now": any(row["passes_descent_acceptance_now"] for row in rows),
    }


def closure_decision(summary: dict[str, Any], witness: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2666_PAIR12_SELECTOR_WITNESS_TO_BOUNDARY_PHASE_SECTOR_DESCENT_AUDIT__ORIENTATION_MAP_NOT_EXPORTED",
        "professorial_verdict": "P2666 takes the next honest step after P2665 by using existing pair1/pair2 selector-witness packets instead of inventing a new selector.  The computation shows a real near-miss: the exported w_break lane separates pair1/pair2 and pair2 is the positive branch, so a convention could map pair2 to boundary sector 1.  But the opposite convention is equally available, the actual source-topology witness remains preLM/not pair12-typed, and the typed descent interface is future-route only.  Thus no current strict, non-declared pair12 -> boundary-phase holonomy-sector descent is exported, and no boundary-phase bit target, UV unit, beta source, QW-2191 discharge, L_total reopening, or ToE closure follows.",
        "next_honest_step": "Attempt one sharply scoped theorem: derive a canonical orientation map pair2_positive -> boundary_sector_1 from bridge-completed nadsoliton dynamics, not from naming convention.  The proof must cite the exported pair12 witness split, construct a typed pair12 -> boundary-cycle functor, prove convention reversal is forbidden internally, and then rerun P2662/P2665.  If convention reversal cannot be forbidden, record a no-go: selector witness splitting can at most provide a two-branch carrier, not the entropy-bit target source.",
        "pair12_split_exported": summary["p731_w_break_separates_pair12"],
        "pair2_positive_branch": summary["p731_pair2_positive"],
        "mathematically_possible_sector_one_mapping_exists": witness["raw_holonomy_sign_preserved_by_at_least_one_mapping"],
        "passing_descent_mapping_count": len(witness["passing_descent_mappings"]),
        "strict_current_descent_exported_now": witness["strict_current_descent_exported_now"],
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    decision = payload["closure_decision"]
    witness = payload["descent_mapping_witness"]
    lines = [
        "# P2666/S1616 pair12 selector-witness to boundary-phase sector descent audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines += [
        "",
        "## Prior-summary audit",
    ]
    for key, value in payload["prior_summary_audit"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Descent mapping witness",
        witness["statement"],
        f"Orientation convention degeneracy count: `{witness['orientation_convention_degeneracy_count']}`.",
        f"Raw holonomy sign preserved by at least one mapping? `{witness['raw_holonomy_sign_preserved_by_at_least_one_mapping']}`.",
        f"Strict current descent exported now? `{witness['strict_current_descent_exported_now']}`.",
        "",
        "| mapping | selected sector | selects sector 1? | pair12 split? | source pair12 typed? | typed boundary map? | future target current? | passes? |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in witness["rows"]:
        lines.append(
            f"| `{row['candidate_mapping']}` | `{row['selected_boundary_sector']}` | `{row['selects_nonexact_boundary_sector_one']}` | "
            f"`{row['pair12_witness_split_exported']}` | `{row['current_source_witness_pair12_typed']}` | "
            f"`{row['pair12_to_boundary_typed_map_exported']}` | `{row['future_target_counts_as_current_source']}` | `{row['passes_descent_acceptance_now']}` |"
        )
    lines += [
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Passing descent mapping count: `{decision['passing_descent_mapping_count']}`.",
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
    summary = prior_summary_audit()
    witness = descent_mapping_witness(summary)
    decision = closure_decision(summary, witness)
    payload: dict[str, Any] = {
        "status": "P2666_PAIR12_SELECTOR_WITNESS_TO_BOUNDARY_PHASE_SECTOR_DESCENT_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": audit,
        "source_hashes": {"P2665": sha256_file(P2665), "P731": sha256_file(P731), "P741": sha256_file(P741), "P764": sha256_file(P764), "F947": sha256_file(F947), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "prior_summary_audit": summary,
        "typed_descent_target": {"domain": PAIR_BRANCHES, "codomain": BOUNDARY_SECTORS, "ontology_order": "nadsoliton -> light -> matter -> emergent observer", "no_sub_nadsoliton_information_layer": True},
        "descent_mapping_witness": witness,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2666/S1616 pair12 selector-witness to boundary-phase sector descent guard", "## P2666/S1616 pair12 selector-witness to boundary-phase sector descent guard\n\n`P2666/S1616` reuses existing selector-source artifacts (`w_break`, pair1/pair2 witness split, actual source-topology witness, and typed-descent target summaries) to audit a concrete pair12 -> boundary-phase sector descent.  The near-miss is real: pair2 is the positive witness branch and can be conventionally mapped to boundary sector `1`; however the reverse convention is equally mathematical, the actual source witness remains preLM/not pair12-typed, and the typed descent interface is future-route only.  Therefore no strict current pair12 -> boundary-phase holonomy-sector map, boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2666/S1616 pair12 selector-witness to boundary-phase sector Ltotal guard", "## P2666/S1616 pair12 selector-witness to boundary-phase sector Ltotal guard\n\n`P2666/S1616` does not license a selector-derived boundary-phase entropy term in `L_total`.  A Lagrangian coefficient may use the pair12 witness split only after a typed pair12 -> boundary-cycle functor fixes the orientation convention internally and exports a current, non-declared descent to boundary sector `1`.\n")
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
