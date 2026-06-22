#!/usr/bin/env python3
"""P3003/S1953: square-map basin/attractor localizer obstruction.

P3002 left three square-map-functional-graph routes. This audit attacks exactly
one: a nonpremise basin/attractor localizer for the Z/12Z squaring functional
graph. It does not replay strict provenance, source-atom coupling, action
installation, zero-divisor graph, CRT idempotent, nilradical, annihilator,
Fourier, selector, bridge, role-transfer, or L_total lanes.

The finite calculation constructs exact localizer receivers: basin rows,
attractor rows, fiber rows, unit-orbit signatures, translation-stability tests,
and singleton/uniqueness witnesses. These are useful graph-dynamical labels, but
they remain algebraic bookkeeping without a strict provenance theorem or a
nonpremise physical sector/localizer theorem.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3001_s1951_z12_square_map_functional_graph_source_candidate_obstruction import MODULUS, square_map_witness
from p3002_s1952_square_map_strict_provenance_obstruction import OUT as P3002, UNITS

OUT = GEN / "p3003_s1953_square_map_basin_attractor_localizer_obstruction.json"
MD = GEN / "p3003_s1953_square_map_basin_attractor_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def unit_orbit(x: int) -> list[int]:
    return sorted({(u * x) % MODULUS for u in UNITS})


def localizer_witness() -> dict[str, Any]:
    graph = square_map_witness()
    basins = {int(k): v for k, v in graph["basins"].items()}
    fibers = {int(k): v for k, v in graph["fibers"].items()}
    attractors = graph["fixed_points"]
    basin_rows = []
    for a in attractors:
        members = basins[a]
        rows = [r for r in graph["row_data"] if r["fixed_attractor"] == a]
        basin_rows.append({
            "attractor": a,
            "basin_members": members,
            "basin_size": len(members),
            "depth_signature": sorted([r["depth_to_attractor"] for r in rows]),
            "fiber_to_attractor": fibers[a],
            "unit_orbit_of_attractor": unit_orbit(a),
            "singleton_unit_orbit": len(unit_orbit(a)) == 1,
            "singleton_basin": len(members) == 1,
            "nonpremise_physical_sector_exported": False,
        })
    residue_rows = []
    for r in graph["row_data"]:
        x = r["residue"]
        residue_rows.append({
            "residue": x,
            "square_target": r["square_target"],
            "fixed_attractor": r["fixed_attractor"],
            "basin_size": len(basins[r["fixed_attractor"]]),
            "depth_to_attractor": r["depth_to_attractor"],
            "fiber_size_of_target": len(fibers[r["square_target"]]),
            "unit_orbit": unit_orbit(x),
            "unit_orbit_size": len(unit_orbit(x)),
            "localizer_signature": [r["fixed_attractor"], len(basins[r["fixed_attractor"]]), r["depth_to_attractor"], len(fibers[r["square_target"]]), len(unit_orbit(x))],
            "accepted_nonpremise_localizer": False,
        })
    signatures: dict[str, list[int]] = {}
    for row in residue_rows:
        signatures.setdefault(json.dumps(row["localizer_signature"]), []).append(row["residue"])
    translation_rows = []
    basin_sets = {a: set(v) for a, v in basins.items()}
    for t in range(MODULUS):
        translated = {a: {((x + t) % MODULUS) for x in members} for a, members in basin_sets.items()}
        preserves_partition = set(map(frozenset, translated.values())) == set(map(frozenset, basin_sets.values()))
        translation_rows.append({"translation": t, "preserves_basin_partition": preserves_partition})
    return {
        "basin_rows": basin_rows,
        "residue_rows": residue_rows,
        "signature_classes": {k: v for k, v in sorted(signatures.items(), key=lambda item: item[1])},
        "signature_class_count": len(signatures),
        "singleton_signature_classes": [v[0] for v in signatures.values() if len(v) == 1],
        "translation_rows": translation_rows,
        "basin_partition_preserving_translations": [r["translation"] for r in translation_rows if r["preserves_basin_partition"]],
        "strict_provenance_exported_by_P3002": False,
        "nonpremise_physical_sector_theorem_exported": False,
        "accepted_nonpremise_basin_attractor_localizers": [],
    }


def obligation_rows(w: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_basin_attractor_receivers", "satisfied": len(w["basin_rows"]) == 4 and len(w["residue_rows"]) == 12, "evidence": "four attractor/basin rows and twelve residue-localizer rows are computed"},
        {"obligation": "signature_partition_computed", "satisfied": w["signature_class_count"] >= 4, "evidence": f"computed {w['signature_class_count']} localizer signature classes"},
        {"obligation": "translation_nonlocalizer_witness", "satisfied": w["basin_partition_preserving_translations"] == [0, 3, 6, 9], "evidence": "four translations preserve the basin partition, so the partition is not an origin/localizer"},
        {"obligation": "strict_provenance_available", "satisfied": w["strict_provenance_exported_by_P3002"], "evidence": "P3002 found no strict nadsoliton source map for the square map"},
        {"obligation": "nonpremise_physical_sector_theorem", "satisfied": w["nonpremise_physical_sector_theorem_exported"], "evidence": "basin/attractor signatures are algebraic labels, not exported physical sectors"},
        {"obligation": "accepted_current_nonpremise_localizer", "satisfied": bool(w["accepted_nonpremise_basin_attractor_localizers"]), "evidence": "no row satisfies provenance plus a nonpremise physical sector theorem"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_basin_receivers", "signature_partition", "translation_witness", "strict_provenance", "physical_sector_theorem", "nonpremise_localizer_export"]
    return [{"present": dict(zip(names, bits)), "accepts_nonpremise_localizer": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3002_path: Any) -> dict[str, Any]:
    witness = localizer_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3003_SQUARE_MAP_BASIN_ATTRACTOR_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3002": hashlib.sha256(p3002_path.read_bytes()).hexdigest() if p3002_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "SquareMapBasinAttractorLocalizer_ObstructionMatrix",
            "localizer_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "localizer_certificate": {
            "basin_row_count": len(witness["basin_rows"]),
            "residue_row_count": len(witness["residue_rows"]),
            "signature_class_count": witness["signature_class_count"],
            "singleton_signature_classes": witness["singleton_signature_classes"],
            "basin_partition_preserving_translations": witness["basin_partition_preserving_translations"],
            "accepted_nonpremise_localizer_count": len(witness["accepted_nonpremise_basin_attractor_localizers"]),
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_nonpremise_localizer"]),
        },
        "decision": {
            "positive_progress": "P3003 attacks exactly one P3002 remaining square-map route: a nonpremise basin/attractor localizer for the Z/12Z square-map functional graph.",
            "breakthrough": "Bounded no-go: exact basin, attractor, fiber, unit-orbit, signature, and translation-stability receivers exist, but no strict provenance or nonpremise physical sector/localizer theorem is exported.",
            "negative_export_flags": {k: False for k in ["nonpremise_basin_attractor_localizer_exported", "square_map_strict_provenance_exported", "named_source_atom_coupling_exported", "unit_bearing_action_installation_exported", "selector_closure_exported", "bridge_closure_exported", "role_transfer_exported", "nonproxy_ltotal_exported", "toe_closure_exported"]},
            "next_honest_step": "A next admissible square-map-functional-graph move may attack exactly one remaining route: named source-atom coupling or unit-bearing action installation. Do not replay basin/attractor signatures, multiplicative compatibility, CRT idempotents, zero-divisor graph, selector, bridge, role-transfer, or L_total lanes as localization.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["localizer_certificate"]
    lines = [
        "# P3003/S1953 square-map basin/attractor localizer obstruction", "",
        f"Status: `{payload['status']}`", "", "## Localizer certificate",
        f"- basin rows: `{cert['basin_row_count']}`",
        f"- residue rows: `{cert['residue_row_count']}`",
        f"- signature classes: `{cert['signature_class_count']}`",
        f"- singleton signature classes: `{cert['singleton_signature_classes']}`",
        f"- basin-partition-preserving translations: `{cert['basin_partition_preserving_translations']}`",
        f"- accepted nonpremise localizers: `{cert['accepted_nonpremise_localizer_count']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3002)
    payload = build_payload(P3002)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3003/S1953 square-map basin/attractor localizer obstruction", "## P3003/S1953 square-map basin/attractor localizer obstruction\n\n`P3003/S1953` attacks exactly one remaining P3002 route: a nonpremise basin/attractor localizer for the `Z/12Z` square-map functional graph.  The finite localizer receivers are exact: four basin/attractor rows, twelve residue rows, computed fiber/unit-orbit signatures, and translations `0,3,6,9` preserve the basin partition.  These are algebraic graph-dynamics signatures only; P3002 exported no strict provenance and no nonpremise physical sector/localizer theorem is present.  No named source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3003/S1953 square-map basin/localizer `L_total` guard", "## P3003/S1953 square-map basin/localizer `L_total` guard\n\n`P3003/S1953` adds no square-map basin or attractor term to `L_total`.  Basin sizes, fibers, unit-orbit signatures, singleton signatures, and translation-stability witnesses do not supply strict field provenance, a named density theorem, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current square-map basin/attractor localizer obstruction guardrail (P3003/S1953, 2026-06-21)", "## Current square-map basin/attractor localizer obstruction guardrail (P3003/S1953, 2026-06-21)\n\n- P3003 attacks exactly one remaining P3002 route: a nonpremise basin/attractor localizer for the `Z/12Z` square-map functional graph.\n- Finite positives are exact but graph-dynamical only: four basin/attractor rows, twelve residue rows, computed fiber/unit-orbit signatures, singleton signatures, and translations `0,3,6,9` preserving the basin partition.\n- The route is bounded no-go because P3002 exported no strict provenance and no nonpremise physical sector/localizer theorem is present.\n- Do not promote basin/attractor signatures, singleton signatures, fibers, unit orbits, translation witnesses, CRT idempotent replay, zero-divisor graph replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  A next admissible square-map-functional-graph move may attack exactly one remaining route, or else introduce a genuinely new strict typed object/provider while preserving the P2929-P3003 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
