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
OUT = GEN / "p2667_s1617_pair12_boundary_orientation_reversal_no_go_audit.json"
MD = GEN / "p2667_s1617_pair12_boundary_orientation_reversal_no_go_audit.md"
P2666 = GEN / "p2666_s1616_pair12_selector_witness_to_boundary_phase_sector_descent_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

PAIR_BRANCHES = ("pair1", "pair2")
SECTORS = (0, 1)
NEGATIVE_EXPORT_FLAGS = [
    "canonical_pair12_boundary_orientation_map_exported",
    "orientation_reversal_forbidden_internally",
    "pair12_to_boundary_phase_sector_descent_exported",
    "boundary_phase_bit_target_exported_unconditionally",
    "intrinsic_entropy_level_exported",
    "bit_to_action_map_sourced_unconditionally",
    "bit_to_length_map_sourced_unconditionally",
    "unique_beta_representative_selected_unconditionally",
    "target_independent_beta_source_exported",
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
        "orientation_convention_content": "orientation convention|orientation map|orientation datum|cycle orientation|convention reversal",
        "pair12_positive_branch_content": "pair2.*positive|pair1/pair2|w_break|witness split|pair12",
        "boundary_cycle_functor_content": "boundary-cycle functor|boundary cycle|boundary-sector|boundary sector|square holonomy|holonomy sign",
        "symmetry_reversal_content": "reversal|swap|Z2|residual sign|sign-gauge|projective|ray-level",
        "nonclosure_guard_content": "QW-2191|role-bearing L_total|ToE closure|beta source|nonexport|future-route",
    }
    return {"tool": "rg", "mode": "content-first canonical orientation-map / convention-reversal audit", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def upstream_consistency() -> dict[str, Any]:
    p2666 = load_json(P2666)
    decision = p2666.get("closure_decision", {})
    return {
        "p2666_pair12_split_exported": decision.get("pair12_split_exported") is True,
        "p2666_pair2_positive": decision.get("pair2_positive_branch") is True,
        "p2666_possible_sector_one_mapping_exists": decision.get("mathematically_possible_sector_one_mapping_exists") is True,
        "p2666_no_passing_descent_mapping": decision.get("passing_descent_mapping_count") == 0,
        "p2666_strict_current_descent_not_exported": decision.get("strict_current_descent_exported_now") is False,
    }


def invariant_signature(mapping: dict[str, int]) -> dict[str, Any]:
    # Current exported invariants know: pair2 is the positive pair branch; sector 1 is the non-exact bit sector.
    # They do not export a cross-invariant tying pair2 to sector 1.  Thus both mappings have identical intrinsic data up to relabeling.
    positive_pair = "pair2"
    nonexact_sector = 1
    return {
        "positive_pair_branch": positive_pair,
        "nonexact_boundary_sector": nonexact_sector,
        "selected_sector_for_positive_pair": mapping[positive_pair],
        "selects_nonexact_sector": mapping[positive_pair] == nonexact_sector,
        "cross_invariant_pair2_equals_sector1_exported": False,
        "boundary_sector_labels_are_intrinsic_without_orientation": False,
    }


def orientation_reversal_witness() -> dict[str, Any]:
    mappings = []
    for perm in itertools.permutations(SECTORS):
        mapping = dict(zip(PAIR_BRANCHES, perm, strict=True))
        sig = invariant_signature(mapping)
        mappings.append({"mapping": mapping, **sig})
    reversal_pairs = []
    for a, b in itertools.combinations(mappings, 2):
        reversed_by_sector_swap = a["mapping"]["pair1"] == b["mapping"]["pair2"] and a["mapping"]["pair2"] == b["mapping"]["pair1"]
        same_exported_intrinsic_data = (
            a["positive_pair_branch"] == b["positive_pair_branch"]
            and a["nonexact_boundary_sector"] == b["nonexact_boundary_sector"]
            and not a["cross_invariant_pair2_equals_sector1_exported"]
            and not b["cross_invariant_pair2_equals_sector1_exported"]
        )
        reversal_pairs.append({
            "mapping_a": a["mapping"],
            "mapping_b": b["mapping"],
            "reversed_by_sector_swap": reversed_by_sector_swap,
            "same_exported_intrinsic_data_without_cross_invariant": same_exported_intrinsic_data,
            "orientation_reversal_forbidden": False,
        })
    return {
        "statement": "P2667 audits the exact missing theorem from P2666: can current data forbid convention reversal and canonically identify pair2_positive with boundary_sector_1?  Enumerating the two possible orientation maps shows that one selects sector 1, but the sector-swap reversal preserves all currently exported intrinsic data because no cross-invariant ties pair2 positivity to the boundary sector label.",
        "candidate_mappings": mappings,
        "reversal_pairs": reversal_pairs,
        "candidate_mapping_count": len(mappings),
        "sector_one_mapping_exists": any(row["selects_nonexact_sector"] for row in mappings),
        "all_reversal_pairs_unforbidden": all(not row["orientation_reversal_forbidden"] for row in reversal_pairs),
        "canonical_orientation_map_exported_now": False,
    }


def source_candidate_matrix(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"candidate": "pair2_positive_label_as_boundary_sector_1", "selects_sector_1": True, "forbids_reversal": False, "passes_now": False, "verdict": "false pass: it chooses the desired convention but does not source the cross-invariant"},
        {"candidate": "boundary_sector_label_1_as_intrinsic_orientation", "selects_sector_1": True, "forbids_reversal": False, "passes_now": False, "verdict": "blocked: sector labels are not physical without an internally oriented boundary-cycle functor"},
        {"candidate": "sector_swap_reversal_symmetry", "selects_sector_1": None, "forbids_reversal": False, "passes_now": False, "verdict": "negative control: reversal remains allowed by current exported data"},
        {"candidate": "future_cross_invariant_pair2_to_sector1", "selects_sector_1": True, "forbids_reversal": None, "passes_now": False, "verdict": "open theorem target: derive a bridge-completed cross-invariant that makes reversal impossible"},
    ]


def closure_decision(witness: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_now"]]
    return {
        "decision": "P2667_PAIR12_BOUNDARY_ORIENTATION_REVERSAL_NO_GO_AUDIT__NO_CANONICAL_ORIENTATION_MAP",
        "professorial_verdict": "P2667 is the sharper proof-grade obstruction requested by P2666.  The computation confirms that pair2-positive can be mapped to boundary sector 1, but the opposite sector convention remains equally consistent with all currently exported intrinsic data.  Since no cross-invariant or internally oriented boundary-cycle functor forbids sector swap, the desired orientation map would still be a convention.  Therefore no pair12-to-boundary source theorem, boundary-phase bit target, UV unit, beta source, QW-2191 discharge, L_total reopening, or ToE closure follows.",
        "next_honest_step": "Do not rerun P2662 as sourced yet.  The next honest proof target is a bridge-completed cross-invariant: an internally oriented boundary-cycle functor whose value changes under sector swap and is tied to the pair2-positive witness branch.  If such a cross-invariant cannot be produced, freeze the pair12 selector-witness route as a two-branch carrier no-go for entropy-bit sourcing.",
        "passing_orientation_source_candidates": passing,
        "sector_one_mapping_exists": witness["sector_one_mapping_exists"],
        "all_reversal_pairs_unforbidden": witness["all_reversal_pairs_unforbidden"],
        "canonical_orientation_map_exported_now": witness["canonical_orientation_map_exported_now"],
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["orientation_reversal_witness"]
    decision = payload["closure_decision"]
    lines = [
        "# P2667/S1617 pair12-boundary orientation reversal no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines += ["", "## Upstream consistency"]
    for key, value in payload["upstream_consistency"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Orientation reversal witness",
        witness["statement"],
        f"Candidate orientation maps: `{witness['candidate_mapping_count']}`.",
        f"A sector-one mapping exists? `{witness['sector_one_mapping_exists']}`.",
        f"All reversal pairs remain unforbidden? `{witness['all_reversal_pairs_unforbidden']}`.",
        f"Canonical orientation map exported now? `{witness['canonical_orientation_map_exported_now']}`.",
        "",
        "| mapping | selected sector for pair2 | selects sector 1? | cross-invariant exported? | sector labels intrinsic? |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in witness["candidate_mappings"]:
        lines.append(f"| `{row['mapping']}` | `{row['selected_sector_for_positive_pair']}` | `{row['selects_nonexact_sector']}` | `{row['cross_invariant_pair2_equals_sector1_exported']}` | `{row['boundary_sector_labels_are_intrinsic_without_orientation']}` |")
    lines += [
        "",
        "## Source candidate matrix",
        "| candidate | selects sector 1? | forbids reversal? | passes now? | verdict |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for row in payload["source_candidate_matrix"]:
        lines.append(f"| `{row['candidate']}` | `{row['selects_sector_1']}` | `{row['forbids_reversal']}` | `{row['passes_now']}` | {row['verdict']} |")
    lines += [
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Passing orientation-source candidates: `{decision['passing_orientation_source_candidates']}`.",
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
    upstream = upstream_consistency()
    witness = orientation_reversal_witness()
    matrix = source_candidate_matrix(witness)
    decision = closure_decision(witness, matrix)
    payload: dict[str, Any] = {
        "status": "P2667_PAIR12_BOUNDARY_ORIENTATION_REVERSAL_NO_GO_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": audit,
        "source_hashes": {"P2666": sha256_file(P2666), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream,
        "orientation_reversal_witness": witness,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2667/S1617 pair12-boundary orientation reversal no-go guard", "## P2667/S1617 pair12-boundary orientation reversal no-go guard\n\n`P2667/S1617` audits the exact orientation-map gap left by P2666.  Mapping `pair2_positive` to boundary sector `1` is mathematically possible, but sector-swap reversal preserves all currently exported intrinsic data because no cross-invariant or internally oriented boundary-cycle functor ties the positive pair12 branch to the boundary-sector label.  Therefore the orientation remains conventional, and no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2667/S1617 pair12-boundary orientation reversal Ltotal guard", "## P2667/S1617 pair12-boundary orientation reversal Ltotal guard\n\n`P2667/S1617` keeps `L_total` closed to pair12-derived boundary-phase entropy terms.  A future variational coefficient must first derive a cross-invariant that forbids the sector-swap convention and internally orients the boundary-cycle functor; naming `pair2_positive -> sector_1` is not sufficient.\n")
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
