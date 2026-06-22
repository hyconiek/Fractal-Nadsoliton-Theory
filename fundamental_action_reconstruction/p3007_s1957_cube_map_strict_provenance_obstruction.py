#!/usr/bin/env python3
"""P3007/S1957: cube-map strict provenance obstruction.

This is the first narrow follow-up allowed by P3006 for the new cube-map
functional graph.  It attacks exactly one missing theorem: strict provenance of
the map x -> x^3 mod 12.  The audit does not replay square-map, zero-divisor,
CRT-idempotent, nilradical, annihilator, Fourier, selector, bridge,
role-transfer, or L_total lanes as closure evidence.

The constructed object is a finite provenance receiver matrix.  It computes the
exact CRT defect of x^3 from identity on Z/12Z, tests unit/fixed-sector facts,
and then separates algebraic availability from strict-source export.  The
outcome is bounded no-go: the cube map has exact finite algebraic witnesses, but
current artifacts export no theorem that sources this cubing law from the
nadsoliton, APD/selector data, damping/compression, or variational dynamics.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3006_s1956_z12_cube_map_functional_graph_source_candidate_obstruction import OUT as P3006, cube_map

OUT = GEN / "p3007_s1957_cube_map_strict_provenance_obstruction.json"
MD = GEN / "p3007_s1957_cube_map_strict_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12


def crt_defect_rows() -> list[dict[str, Any]]:
    rows = []
    for x in range(MODULUS):
        y = cube_map(x)
        defect = (y - x) % MODULUS
        rows.append({
            "residue": x,
            "cube_target": y,
            "defect_mod12": defect,
            "defect_mod3": defect % 3,
            "defect_mod4": defect % 4,
            "fixed": y == x,
            "unit": __import__("math").gcd(x, MODULUS) == 1,
            "even_mod4_defect_carrier": defect % 4 != 0,
        })
    return rows


def provenance_receiver_matrix(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    positive_algebra = [
        ("finite_z12_power_map", True, "x -> x^3 mod 12 is exactly computable for all 12 residues"),
        ("crt_mod3_identity", all(r["defect_mod3"] == 0 for r in rows), "Fermat/finite check gives x^3 == x mod 3 for every residue"),
        ("unit_fixed_sector", all(r["fixed"] for r in rows if r["unit"]), "all four units 1,5,7,11 are fixed by cubing"),
        ("nonidentity_defect_witness", [r["residue"] for r in rows if not r["fixed"]] == [2, 6, 10], "the only nonfixed residues are 2,6,10"),
        ("crt_mod4_even_defect_witness", [r["residue"] for r in rows if r["even_mod4_defect_carrier"]] == [2, 6, 10], "the nonidentity part is carried by the mod-4 even sector"),
    ]
    missing_strict = [
        ("strict_nadsoliton_cube_law_export", False, "no current theorem derives x -> x^3 from nadsoliton source data"),
        ("apd_or_phase_provenance_export", False, "no APD/phase/topological-bit export selects cubing rather than another finite power map"),
        ("damping_compression_provenance_export", False, "strict damping/compression artifacts do not source this residue cubing law"),
        ("nonpremise_sector_localizer_export", False, "fixed/unit/even-sector labels are ring bookkeeping, not physical sectors"),
        ("unit_bearing_density_or_action_export", False, "no unit-bearing density, measure, or nonproxy action installation is exported"),
        ("accepted_current_strict_provenance", False, "algebraic receivers are insufficient without a strict source theorem"),
    ]
    return [{"receiver": name, "satisfied": ok, "evidence": ev} for name, ok, ev in positive_algebra + missing_strict]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = [
        "finite_power_map", "crt_identity_data", "nonidentity_defect_witness",
        "strict_nadsoliton_source", "apd_phase_source", "damping_or_variational_source",
        "nonpremise_localizer", "unit_bearing_export", "nonproxy_export",
    ]
    return [{"present": dict(zip(names, bits)), "accepts_strict_provenance": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3006_path: Any) -> dict[str, Any]:
    rows = crt_defect_rows()
    matrix = acceptance_matrix()
    moved = [r["residue"] for r in rows if not r["fixed"]]
    return {
        "status": "P3007_CUBE_MAP_STRICT_PROVENANCE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3006": hashlib.sha256(p3006_path.read_bytes()).hexdigest() if p3006_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "CubeMapStrictProvenance_ReceiverObstructionMatrix",
            "crt_defect_rows": rows,
            "provenance_receiver_matrix": provenance_receiver_matrix(rows),
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "residue_count": len(rows),
            "fixed_residues": [r["residue"] for r in rows if r["fixed"]],
            "moved_residues": moved,
            "unit_residues": [r["residue"] for r in rows if r["unit"]],
            "all_defects_zero_mod3": all(r["defect_mod3"] == 0 for r in rows),
            "mod4_defect_carriers": [r["residue"] for r in rows if r["even_mod4_defect_carrier"]],
            "strict_provenance_receivers_accepted": [],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_strict_provenance"]),
        },
        "decision": {
            "positive_progress": "P3007 attacks exactly one P3006 missing theorem: strict provenance for the cube-map x -> x^3 mod 12.",
            "breakthrough": "Bounded no-go: exact CRT/unit/defect receivers exist, but they remain finite algebraic bookkeeping; no current artifact exports a strict nadsoliton, APD/phase, damping, or variational source theorem for cubing.",
            "negative_export_flags": {k: False for k in ["cube_map_strict_provenance_exported", "nonpremise_basin_fixed_sector_localizer_exported", "named_source_atom_coupling_exported", "unit_bearing_action_installation_exported", "selector_closure_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "A next cube-map-functional-graph move may attack exactly one remaining missing theorem for this object: nonpremise basin/fixed-sector localizer, named source-atom coupling, or unit-bearing action installation. Do not replay strict provenance, square-map, zero-divisor graph, CRT idempotent, nilradical, annihilator, Fourier, selector, bridge, role-transfer, or L_total lanes through cube-map labels; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P3007 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P3007/S1957 cube-map strict provenance obstruction", "",
        f"Status: `{payload['status']}`", "", "## Provenance certificate",
        f"- residue count: `{cert['residue_count']}`",
        f"- fixed residues: `{cert['fixed_residues']}`",
        f"- moved residues: `{cert['moved_residues']}`",
        f"- unit residues: `{cert['unit_residues']}`",
        f"- all defects zero mod 3: `{cert['all_defects_zero_mod3']}`",
        f"- mod-4 defect carriers: `{cert['mod4_defect_carriers']}`",
        f"- strict provenance receivers accepted: `{cert['strict_provenance_receivers_accepted']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3006)
    payload = build_payload(P3006)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3007/S1957 cube-map strict provenance obstruction", "## P3007/S1957 cube-map strict provenance obstruction\n\n`P3007/S1957` attacks exactly one P3006 missing theorem: strict provenance for the cube-map `x -> x^3 mod 12`.  The exact finite receivers are positive as algebra: `x^3-x` vanishes mod `3` for all residues, all units `[1, 5, 7, 11]` are fixed, and the only nonidentity residues are `[2, 6, 10]`, carried by the mod-4 even-sector defect.  This remains bounded no-go for strict provenance: no current artifact exports a nadsoliton source law, APD/phase/topological-bit source, damping/compression source, nonpremise sector localizer, unit-bearing density, nonproxy `L_total`, bridge closure, role transfer, or ToE from these cube-map receivers.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3007/S1957 cube-map provenance `L_total` guard", "## P3007/S1957 cube-map provenance `L_total` guard\n\n`P3007/S1957` adds no cube-map provenance term to `L_total`.  The CRT defect rows, fixed-unit rows, mod-3 identity witness, and mod-4 even-sector defect witness are finite algebraic provenance receivers only; they do not supply a strict field source, named density theorem, unit-bearing measure, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current cube-map strict provenance obstruction guardrail (P3007/S1957, 2026-06-22)", "## Current cube-map strict provenance obstruction guardrail (P3007/S1957, 2026-06-22)\n\n- P3007 attacks exactly one P3006 missing theorem: strict provenance for the cube-map `x -> x^3 mod 12`.\n- Exact finite positives are algebraic only: `x^3-x` is zero mod `3` for every residue, all unit residues `[1, 5, 7, 11]` are fixed, and the nonidentity residues `[2, 6, 10]` are exactly the mod-4 even-sector defect carriers.\n- The current route is bounded no-go for strict provenance: no nadsoliton source law, APD/phase/topological-bit source, damping/compression source, nonpremise sector localizer, unit-bearing density, or nonproxy export is present.\n- Do not promote cube-map CRT defects, fixed-unit rows, mod-3 identity rows, mod-4 defect carriers, square-map replay, zero-divisor graph replay, CRT-idempotent replay, nilradical/annihilator/Fourier replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  A next admissible cube-map-functional-graph move may attack exactly one remaining missing theorem: nonpremise basin/fixed-sector localizer, named source-atom coupling, or unit-bearing action installation; otherwise introduce a genuinely new strict typed object/provider while preserving the P2929-P3007 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
