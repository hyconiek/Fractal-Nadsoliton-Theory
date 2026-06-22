#!/usr/bin/env python3
"""P3008/S1958: cube-map basin/fixed-sector localizer obstruction.

This is the next narrow cube-map follow-up after P3007. It attacks exactly one
remaining missing theorem: a nonpremise basin/fixed-sector localizer for the
cubing map x -> x^3 mod 12. It does not replay strict provenance, square-map,
zero-divisor graph, CRT-idempotent, nilradical, annihilator, Fourier, selector,
bridge, role-transfer, or L_total lanes.

The constructed object is a finite localizer receiver matrix: basin rows,
residue rows, fixed-sector signatures, fiber signatures, unit-orbit signatures,
and symmetry-stability rows. The exact computation shows that the basin
partition is preserved by every even translation and by every unit action, so
basin/fixed-sector bookkeeping cannot select a unique nonpremise origin or
physical sector on current artifacts.
"""
from __future__ import annotations

import hashlib, json, math
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction import cube_map
from p3007_s1957_cube_map_strict_provenance_obstruction import OUT as P3007

OUT = GEN / "p3008_s1958_cube_map_basin_fixed_sector_localizer_obstruction.json"
MD = GEN / "p3008_s1958_cube_map_basin_fixed_sector_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12
UNITS = [1, 5, 7, 11]


def orbit_path(x: int) -> list[int]:
    path = [x]
    y = cube_map(x)
    if y != x:
        path.append(y)
    return path


def localizer_witness() -> dict[str, Any]:
    residues = list(range(MODULUS))
    fixed = [x for x in residues if cube_map(x) == x]
    basins = {str(fp): sorted(x for x in residues if orbit_path(x)[-1] == fp) for fp in fixed}
    partition = {frozenset(v) for v in basins.values()}
    fibers = {str(y): sorted(x for x in residues if cube_map(x) == y) for y in fixed}
    basin_rows = []
    for fp, members in basins.items():
        m = [int(x) for x in members]
        basin_rows.append({
            "fixed_attractor": int(fp),
            "members": m,
            "size": len(m),
            "parities": sorted({x % 2 for x in m}),
            "gcd_profile": sorted({math.gcd(x, MODULUS) for x in m}),
            "singleton": len(m) == 1,
            "nonpremise_physical_sector": False,
        })
    residue_rows = []
    for x in residues:
        attractor = orbit_path(x)[-1]
        residue_rows.append({
            "residue": x,
            "cube_target": cube_map(x),
            "fixed_attractor": attractor,
            "path": orbit_path(x),
            "basin_size": len(basins[str(attractor)]),
            "fiber_size_at_target": len(fibers[str(cube_map(x))]),
            "unit_orbit": sorted({(u * x) % MODULUS for u in UNITS}),
            "localizer_signature": [len(basins[str(attractor)]), len(fibers[str(cube_map(x))]), x % 2, math.gcd(x, MODULUS)],
            "accepted_nonpremise_localizer": False,
        })
    translation_stability = []
    for t in residues:
        image = {frozenset((x + t) % MODULUS for x in block) for block in partition}
        translation_stability.append({"translation": t, "preserves_basin_partition": image == partition})
    unit_stability = []
    for u in UNITS:
        image = {frozenset((u * x) % MODULUS for x in block) for block in partition}
        unit_stability.append({"unit": u, "preserves_basin_partition": image == partition})
    return {
        "fixed_residues": fixed,
        "basins": basins,
        "fibers_over_fixed_targets": fibers,
        "basin_rows": basin_rows,
        "residue_rows": residue_rows,
        "localizer_signature_classes": {str(tuple(sig)): sorted(r["residue"] for r in residue_rows if r["localizer_signature"] == sig) for sig in sorted({tuple(r["localizer_signature"]) for r in residue_rows})},
        "basin_partition_preserving_translations": [r["translation"] for r in translation_stability if r["preserves_basin_partition"]],
        "basin_partition_preserving_units": [r["unit"] for r in unit_stability if r["preserves_basin_partition"]],
        "translation_stability_rows": translation_stability,
        "unit_stability_rows": unit_stability,
        "accepted_nonpremise_localizers": [],
    }


def proof_obligations(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"obligation": "finite_basin_fixed_sector_receivers", "satisfied": len(witness["basin_rows"]) == 9 and len(witness["residue_rows"]) == 12, "evidence": "nine fixed-attractor basin rows and twelve residue rows are built exactly"},
        {"obligation": "fiber_and_signature_partition", "satisfied": bool(witness["localizer_signature_classes"]), "evidence": "finite fiber sizes and localizer signatures are computed for every residue"},
        {"obligation": "symmetry_nonlocalizer_witness", "satisfied": witness["basin_partition_preserving_translations"] == [0, 2, 4, 6, 8, 10] and witness["basin_partition_preserving_units"] == [1, 5, 7, 11], "evidence": "even translations and all unit actions preserve the basin partition"},
        {"obligation": "strict_provenance_available", "satisfied": False, "evidence": "P3007 bounded strict provenance; this audit may not import it"},
        {"obligation": "nonpremise_physical_sector_theorem", "satisfied": False, "evidence": "basin/fixed-sector labels remain ring-dynamical labels, not accepted physical sectors"},
        {"obligation": "unique_origin_or_sector_localizer", "satisfied": False, "evidence": "even-translation symmetry leaves no unique origin selected by the basin partition"},
        {"obligation": "accepted_current_nonpremise_localizer", "satisfied": False, "evidence": "no row satisfies finite receiver, strict provenance, physical-sector theorem, and uniqueness premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_receivers", "signature_partition", "symmetry_witness", "strict_provenance", "physical_sector_theorem", "unique_origin_selector", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_nonpremise_localizer": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3007_path: Any) -> dict[str, Any]:
    witness = localizer_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3008_CUBE_MAP_BASIN_FIXED_SECTOR_LOCALIZER_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3007": hashlib.sha256(p3007_path.read_bytes()).hexdigest() if p3007_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "CubeMapBasinFixedSectorLocalizer_ObstructionMatrix",
            "localizer_witness": witness,
            "proof_obligation_rows": proof_obligations(witness),
            "finite_acceptance_matrix": matrix,
        },
        "localizer_certificate": {
            "fixed_residue_count": len(witness["fixed_residues"]),
            "basin_count": len(witness["basin_rows"]),
            "residue_row_count": len(witness["residue_rows"]),
            "signature_class_count": len(witness["localizer_signature_classes"]),
            "basin_partition_preserving_translations": witness["basin_partition_preserving_translations"],
            "basin_partition_preserving_units": witness["basin_partition_preserving_units"],
            "accepted_nonpremise_localizers": witness["accepted_nonpremise_localizers"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_nonpremise_localizer"]),
        },
        "decision": {
            "positive_progress": "P3008 attacks exactly one P3007 remaining theorem: a nonpremise basin/fixed-sector localizer for the cube-map x -> x^3 mod 12.",
            "breakthrough": "Bounded no-go: exact basin, fixed-sector, fiber, unit-orbit, and signature receivers exist, but even translations [0,2,4,6,8,10] and all unit actions preserve the basin partition; no strict provenance or nonpremise physical-sector theorem is exported.",
            "negative_export_flags": {k: False for k in ["cube_map_nonpremise_localizer_exported", "strict_provenance_exported", "named_source_atom_coupling_exported", "unit_bearing_action_installation_exported", "selector_closure_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next cube-map-functional-graph move may attack exactly one remaining missing theorem for this object: named source-atom coupling or unit-bearing action installation. Do not replay nonpremise localizer, strict provenance, square-map, zero-divisor graph, CRT idempotent, nilradical, annihilator, Fourier, selector, bridge, role-transfer, or L_total lanes through cube-map labels; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P3008 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["localizer_certificate"]
    lines = [
        "# P3008/S1958 cube-map basin/fixed-sector localizer obstruction", "",
        f"Status: `{payload['status']}`", "", "## Localizer certificate",
        f"- fixed residue count: `{cert['fixed_residue_count']}`",
        f"- basin/residue rows: `{cert['basin_count']}/{cert['residue_row_count']}`",
        f"- signature class count: `{cert['signature_class_count']}`",
        f"- basin-partition-preserving translations: `{cert['basin_partition_preserving_translations']}`",
        f"- basin-partition-preserving units: `{cert['basin_partition_preserving_units']}`",
        f"- accepted nonpremise localizers: `{cert['accepted_nonpremise_localizers']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3007)
    payload = build_payload(P3007)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3008/S1958 cube-map basin/fixed-sector localizer obstruction", "## P3008/S1958 cube-map basin/fixed-sector localizer obstruction\n\n`P3008/S1958` attacks exactly one P3007 remaining theorem: a nonpremise basin/fixed-sector localizer for the cube-map `x -> x^3 mod 12`.  The finite localizer receivers are exact: nine fixed-attractor basin rows, twelve residue rows, fiber and localizer signatures, unit-orbit signatures, and symmetry-stability rows.  The route remains bounded no-go: even translations `[0, 2, 4, 6, 8, 10]` and all unit actions `[1, 5, 7, 11]` preserve the basin partition, and no strict provenance or nonpremise physical-sector theorem is exported.  No named source-atom coupling, unit-bearing action installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3008/S1958 cube-map localizer `L_total` guard", "## P3008/S1958 cube-map localizer `L_total` guard\n\n`P3008/S1958` adds no cube-map basin/fixed-sector term to `L_total`.  Basin rows, fixed-sector rows, fiber signatures, unit-orbit signatures, and even-translation/unit-stability witnesses are finite localizer receivers only; they do not supply strict field provenance, a nonpremise physical-sector theorem, named density theorem, unit-bearing measure, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current cube-map basin/fixed-sector localizer obstruction guardrail (P3008/S1958, 2026-06-22)", "## Current cube-map basin/fixed-sector localizer obstruction guardrail (P3008/S1958, 2026-06-22)\n\n- P3008 attacks exactly one P3007 remaining theorem: a nonpremise basin/fixed-sector localizer for the cube-map `x -> x^3 mod 12`.\n- Exact finite positives are localizer receivers only: nine fixed-attractor basin rows, twelve residue rows, fiber/localizer signatures, unit-orbit signatures, and symmetry-stability rows are computed.\n- The current route is bounded no-go: even translations `[0, 2, 4, 6, 8, 10]` and all unit actions `[1, 5, 7, 11]` preserve the basin partition, and no strict provenance or nonpremise physical-sector theorem is exported.\n- Do not promote cube-map basin/fixed-sector signatures, fibers, unit orbits, even-translation witnesses, strict provenance replay, square-map replay, zero-divisor/CRT-idempotent/nilradical/annihilator/Fourier replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  A next admissible cube-map-functional-graph move may attack exactly one remaining theorem: named source-atom coupling or unit-bearing action installation; otherwise introduce a genuinely new strict typed object/provider while preserving the P2929-P3008 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
