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
OUT = GEN / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.json"
MD = GEN / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.md"
P2662 = GEN / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NODES = ["nadsoliton", "light", "matter", "observer", "compression"]
EDGES = [
    ("nadsoliton", "light"),
    ("light", "observer"),
    ("observer", "matter"),
    ("matter", "nadsoliton"),
    ("light", "compression"),
    ("compression", "observer"),
]
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
TRIANGLE = ("light", "observer", "compression")
CYCLE_SQUARE = [("nadsoliton", "light"), ("light", "observer"), ("observer", "matter"), ("matter", "nadsoliton")]
CYCLE_TRIANGLE = [("light", "observer"), ("light", "compression"), ("compression", "observer")]
CYCLE_TWO_STEP = [("light", "observer"), ("light", "compression"), ("compression", "observer")]
BIT_TARGETS = [0, 1, 2]
NEGATIVE_EXPORT_FLAGS = [
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
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "chain_boundary_phase_content": "chain-level|boundary-phase|boundary.*phase|cocycle|coboundary|phase law",
        "entropy_bit_target_content": "N_bits|bit target|log\\(2\\)|entropy target|entropy level",
        "unit_map_content": "unit map|bit-to-action|bit-to-length|action unit|length unit",
        "nonclosure_guard_content": "QW-2191|role-bearing L_total|bridge completion|role transfer|ToE closure|beta source",
    }
    return {"tool": "rg", "mode": "content-first boundary-phase bit-target audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def edge_value(cochain: tuple[int, ...], edge: tuple[str, str]) -> int:
    if edge in EDGE_INDEX:
        return cochain[EDGE_INDEX[edge]]
    rev = (edge[1], edge[0])
    return cochain[EDGE_INDEX[rev]]


def holonomy(cochain: tuple[int, ...], cycle: list[tuple[str, str]]) -> int:
    return sum(edge_value(cochain, edge) for edge in cycle) % 2


def triangle_flat(cochain: tuple[int, ...]) -> bool:
    return holonomy(cochain, CYCLE_TRIANGLE) == 0


def coboundary(vertex_bits: dict[str, int]) -> tuple[int, ...]:
    return tuple((vertex_bits[a] + vertex_bits[b]) % 2 for a, b in EDGES)


def chain_level_witness() -> dict[str, Any]:
    all_cochains = list(itertools.product([0, 1], repeat=len(EDGES)))
    flat = [c for c in all_cochains if triangle_flat(c)]
    exact = {coboundary(dict(zip(NODES, bits, strict=True))) for bits in itertools.product([0, 1], repeat=len(NODES))}
    rows = []
    for cochain in flat:
        square_h = holonomy(cochain, CYCLE_SQUARE)
        triangle_h = holonomy(cochain, CYCLE_TRIANGLE)
        candidate_n_bits = square_h + triangle_h
        rows.append({
            "cochain": list(cochain),
            "square_holonomy_bit": square_h,
            "filled_triangle_holonomy_bit": triangle_h,
            "candidate_n_bits_from_two_cycles": candidate_n_bits,
            "is_exact_coboundary": cochain in exact,
        })
    unique_targets = sorted({row["candidate_n_bits_from_two_cycles"] for row in rows})
    square_values = sorted({row["square_holonomy_bit"] for row in rows})
    exact_square_values = sorted({holonomy(c, CYCLE_SQUARE) for c in exact})
    return {
        "statement": "Enumerate finite Z2 edge phase cochains on the P2662 typed complex. Imposing flatness on the filled triangle preserves multiple square holonomy sectors. Exact coboundaries give zero holonomy, while non-exact cocycles can carry one bit only after a non-exact sector/cycle representative is selected.",
        "node_order": "nadsoliton -> light -> matter -> emergent observer",
        "no_sub_nadsoliton_information_layer": True,
        "edge_count": len(EDGES),
        "all_z2_cochain_count": len(all_cochains),
        "triangle_flat_cochain_count": len(flat),
        "exact_coboundary_count": len(exact),
        "rows": rows,
        "unique_candidate_n_bits_after_flatness": unique_targets,
        "square_holonomy_values_after_flatness": square_values,
        "exact_coboundary_square_holonomy_values": exact_square_values,
        "chain_law_derives_unique_n_bits_without_sector_choice": len(unique_targets) == 1,
        "nonexact_sector_needed_for_nonzero_bit": 1 in square_values and exact_square_values == [0],
        "boundary_of_boundary_zero_gives_only_flatness_not_entropy_level": True,
    }


def upstream_consistency() -> dict[str, Any]:
    p2662 = load_json(P2662)
    return {
        "p2662_conditional_unique_scale_verified": p2662.get("closure_decision", {}).get("conditional_unique_scale_verified") is True,
        "p2662_boundary_phase_bit_target_missing": "boundary_phase_integer_to_entropy_target" in p2662.get("closure_decision", {}).get("missing_premises", []),
        "p2662_unconditional_uv_unit_not_selected": p2662.get("closure_decision", {}).get("unconditional_uv_unit_selected_now") is False,
        "p2662_qw2191_not_discharged": p2662.get("closure_decision", {}).get("qw2191_discharged_now") is False,
    }


def source_candidate_matrix(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"candidate": "boundary_of_boundary_zero_as_bit_target", "computational_result": "flatness only", "passes_as_boundary_phase_bit_target_now": False, "verdict": "blocked: d^2=0 removes inconsistent filled-triangle phase but does not choose N_bits log(2)"},
        {"candidate": "exact_coboundary_phase_law", "computational_result": f"square holonomy values {witness['exact_coboundary_square_holonomy_values']}", "passes_as_boundary_phase_bit_target_now": False, "verdict": "blocked: gauge-exact phases carry no nonzero cycle bit"},
        {"candidate": "nonexact_square_holonomy_bit", "computational_result": f"available values {witness['square_holonomy_values_after_flatness']}", "passes_as_boundary_phase_bit_target_now": False, "verdict": "promising carrier only: one bit appears in a chosen non-exact sector, but the sector/cycle is not internally selected"},
        {"candidate": "declared_cycle_basis_entropy_target", "computational_result": "selects target only by declaration", "passes_as_boundary_phase_bit_target_now": False, "verdict": "false pass: choosing the cycle basis or entropy level by hand reimports the missing premise"},
    ]


def closure_decision(witness: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [r["candidate"] for r in matrix if r["passes_as_boundary_phase_bit_target_now"]]
    return {
        "decision": "P2663_CHAIN_LEVEL_BOUNDARY_PHASE_BIT_TARGET_AUDIT__NONEXACT_BIT_CARRIER_FOUND_NO_INTRINSIC_TARGET_SOURCE",
        "professorial_verdict": "P2663 makes the next missing P2662 premise computational: a finite chain-level Z2 boundary-phase model can distinguish exact zero-holonomy phases from non-exact one-bit cycle holonomy. This is genuine progress because it identifies the only viable carrier class for an entropy bit target. It is not a source theorem: flatness and d^2=0 do not select a non-exact sector, a preferred cycle, or an entropy level, and exact coboundaries give zero. Therefore no unconditional UV unit, beta source, QW-2191 discharge, L_total reopening, or ToE closure follows.",
        "next_honest_step": "Audit whether the bridge-completed nadsoliton dynamics supplies a non-exact boundary-phase sector selector and preferred cycle functional. If it does, rerun P2662 with the derived N_bits target. If it does not, promote P2663 into a no-go theorem: chain-level boundary/cocycle topology can carry an entropy bit, but cannot source the bit target without an additional internal sector selector.",
        "passing_boundary_phase_bit_target_candidates": passing,
        "nonexact_bit_carrier_found": witness["nonexact_sector_needed_for_nonzero_bit"],
        "unique_n_bits_derived_without_sector_choice": witness["chain_law_derives_unique_n_bits_without_sector_choice"],
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["chain_level_boundary_phase_witness"]
    d = payload["closure_decision"]
    lines = [
        "# P2663/S1613 chain-level boundary-phase bit-target audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first audit",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines += [
        "",
        "## Computational witness",
        w["statement"],
        f"All Z2 cochains: `{w['all_z2_cochain_count']}`.",
        f"Triangle-flat cochains: `{w['triangle_flat_cochain_count']}`.",
        f"Exact coboundaries: `{w['exact_coboundary_count']}`.",
        f"Candidate N_bits values after flatness: `{w['unique_candidate_n_bits_after_flatness']}`.",
        f"Exact coboundary square holonomy values: `{w['exact_coboundary_square_holonomy_values']}`.",
        f"Unique N_bits derived without sector choice? `{w['chain_law_derives_unique_n_bits_without_sector_choice']}`.",
        "",
        "## Source candidate matrix",
        "| candidate | computational result | source now? | verdict |",
        "| --- | --- | ---: | --- |",
    ]
    for row in payload["source_candidate_matrix"]:
        lines.append(f"| `{row['candidate']}` | {row['computational_result']} | `{row['passes_as_boundary_phase_bit_target_now']}` | {row['verdict']} |")
    lines += [
        "",
        "## Verdict",
        d["professorial_verdict"],
        f"Decision: `{d['decision']}`.",
        f"Passing boundary-phase bit-target candidates: `{d['passing_boundary_phase_bit_target_candidates']}`.",
        f"Nonexact bit carrier found? `{d['nonexact_bit_carrier_found']}`.",
        f"Boundary-phase bit target exported now? `{d['boundary_phase_bit_target_exported_now']}`.",
        f"Beta source exported now? `{d['beta_source_exported_now']}`.",
        f"QW-2191 discharged now? `{d['qw2191_discharged_now']}`.",
        f"Role-bearing L_total now? `{d['role_bearing_ltotal_now']}`.",
        f"ToE closure now? `{d['toe_closure_now']}`.",
        "",
        "## Next honest step",
        d["next_honest_step"],
        "",
        "## Negative exports",
    ]
    for k, v in payload["negative_export_flags"].items():
        lines.append(f"- `{k}`: `{v}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    witness = chain_level_witness()
    matrix = source_candidate_matrix(witness)
    decision = closure_decision(witness, matrix)
    payload: dict[str, Any] = {
        "status": "P2663_CHAIN_LEVEL_BOUNDARY_PHASE_BIT_TARGET_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {"P2662": sha256_file(P2662), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream_consistency(),
        "typed_chain_complex": {"nodes": NODES, "edges": EDGES, "filled_triangle": TRIANGLE, "ontology_order": "nadsoliton -> light -> matter -> emergent observer", "no_sub_nadsoliton_information_layer": True},
        "chain_level_boundary_phase_witness": witness,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {k: False for k in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2663/S1613 chain-level boundary-phase bit-target guard", "## P2663/S1613 chain-level boundary-phase bit-target guard\n\n`P2663/S1613` audits the first missing P2662 premise at chain level.  A finite `Z2` boundary-phase model finds a real non-exact one-bit holonomy carrier, while exact coboundaries give zero and filled-triangle flatness only enforces consistency.  The carrier is not yet a source theorem because the repo still lacks an internal selector for the non-exact sector, preferred cycle functional, and entropy target `N_bits log(2)`.  Therefore no intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2663/S1613 chain-level boundary-phase bit-target Ltotal guard", "## P2663/S1613 chain-level boundary-phase bit-target Ltotal guard\n\n`P2663/S1613` does not re-open `L_total`.  Non-exact boundary-phase holonomy may become a future entropy-bit source only if bridge-completed nadsoliton dynamics derives the sector selector and preferred cycle functional.  Until then, adding a boundary-phase entropy term to `L_total` would be a declared target, not an internally sourced variational coefficient.\n")
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
