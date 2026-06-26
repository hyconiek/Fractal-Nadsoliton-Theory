#!/usr/bin/env python3
"""P3111/S2061: symplectic/action phase source-law audit.

P3110 isolated the strongest formal object, `symplectic_phase_area_candidate`,
but left its positive representative and calibration coupling unsourced.  This
step builds the missing source-law object directly: can the self-coupled
nadsoliton's four-bit Shannon cell alpha_geo=4 ln 2=ln(16) fix an internal
phase-area representative through phase periodicity, and does that representative
become a physical action/length/time unit without importing standard physics?
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit import ALPHA_GEO
from p3110_s2060_alpha_geo_dimension_comparison_standard_audit import OUT as P3110

OUT = GEN / "p3111_s2061_symplectic_action_phase_source_law_audit.json"
MD = GEN / "p3111_s2061_symplectic_action_phase_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TARGETS = ("action", "length", "time")
WINDINGS = tuple(range(1, 9))
SCALE_MULTIPLIERS = (Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(4, 1))
GATES = (
    "uses_alpha_geo_four_bit_cell",
    "uses_internal_phase_periodicity",
    "fixes_unique_positive_phase_area",
    "nonconventional_source_law_exported",
    "dimensionful_action_unit_exported",
    "length_time_calibration_exported",
    "standard_physics_import_free",
    "selector_bridge_ltotal_toe_free",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "symplectic_action_phase_source_law": r"symplectic|phase-area|phase area|action phase|source-law|Bohr|periodicity",
        "alpha_geo_shannon_cell": r"alpha_geo|4 ln\(2\)|ln\(16\)|four-bit|4-bit|Shannon|nadsoliton|self-coupled information",
        "unit_calibration_gap": r"physical action unit|length/time|calibration coupling|dimensionful|unit-source|scale representative",
        "blocked_imports_and_closures": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def phase_source_rows() -> list[dict[str, Any]]:
    rows = []
    alpha = ALPHA_GEO
    for winding in WINDINGS:
        representative = (2.0 * math.pi * winding) / alpha
        phase_residual = alpha * representative - 2.0 * math.pi * winding
        rows.append(
            {
                "winding_n": winding,
                "phase_equation": "alpha_geo * A_n = 2*pi*n",
                "phase_area_representative": round(representative, 12),
                "phase_residual_abs": abs(phase_residual),
                "positive": representative > 0,
                "minimal_positive_representative": winding == 1,
                "internal_phase_section_selected": winding == 1,
            }
        )
    return rows


def scale_orbit_rows(minimal_area: float) -> list[dict[str, Any]]:
    rows = []
    for target in TARGETS:
        for multiplier in SCALE_MULTIPLIERS:
            scaled_area = minimal_area * float(multiplier)
            phase_value = ALPHA_GEO * scaled_area
            rows.append(
                {
                    "target": target,
                    "scale_multiplier": f"{multiplier.numerator}/{multiplier.denominator}",
                    "scaled_phase_area": round(scaled_area, 12),
                    "phase_value_over_2pi": round(phase_value / (2.0 * math.pi), 12),
                    "preserves_internal_phase_periodicity": multiplier.denominator == 1,
                    "is_minimal_positive_internal_section": multiplier == Fraction(1, 1),
                    "physical_calibration_supplied": False,
                }
            )
    return rows


def candidate_law_rows() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "minimal_internal_phase_periodicity_law",
            "formula": "A_phi := min{A>0 : alpha_geo*A in 2*pi*Z} = 2*pi/alpha_geo",
            "uses_alpha_geo_four_bit_cell": True,
            "uses_internal_phase_periodicity": True,
            "fixes_unique_positive_phase_area": True,
            "nonconventional_source_law_exported": True,
            "dimensionful_action_unit_exported": False,
            "length_time_calibration_exported": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "exports a dimensionless internal phase-area section, not a physical action quantum or length/time calibration",
        },
        {
            "candidate": "integer_winding_family_law",
            "formula": "A_n = 2*pi*n/alpha_geo for n>=1",
            "uses_alpha_geo_four_bit_cell": True,
            "uses_internal_phase_periodicity": True,
            "fixes_unique_positive_phase_area": False,
            "nonconventional_source_law_exported": True,
            "dimensionful_action_unit_exported": False,
            "length_time_calibration_exported": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "without minimality it is a winding family rather than a section",
        },
        {
            "candidate": "imported_bohr_sommerfeld_hbar_law",
            "formula": "S/hbar in 2*pi*Z with hbar identified externally",
            "uses_alpha_geo_four_bit_cell": True,
            "uses_internal_phase_periodicity": True,
            "fixes_unique_positive_phase_area": True,
            "nonconventional_source_law_exported": False,
            "dimensionful_action_unit_exported": True,
            "length_time_calibration_exported": False,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "blocker": "dimensionful action appears only after importing hbar/standard quantum calibration",
        },
    ]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": row["candidate"],
            "required_gate": gate,
            "gate_passed": bool(row[gate]),
            "detail": "passed" if row[gate] else row["blocker"],
        }
        for row in candidates
        for gate in GATES
    ]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidates = sorted({row["candidate"] for row in gates})
    return [
        {
            "candidate": candidate,
            "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate),
            "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate),
            "accepted_physical_unit_source_law": all(row["gate_passed"] for row in gates if row["candidate"] == candidate),
        }
        for candidate in candidates
    ]


def build_payload() -> dict[str, Any]:
    p3110 = read_json(P3110)
    greps = content_grep()
    phases = phase_source_rows()
    minimal = next(row for row in phases if row["minimal_positive_representative"])
    orbits = scale_orbit_rows(float(minimal["phase_area_representative"]))
    candidates = candidate_law_rows()
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted_physical = [row for row in aggs if row["accepted_physical_unit_source_law"]]
    obligations = [
        {"obligation": "read_p3110_next_atom", "satisfied": True, "detail": "P3110 selected the symplectic/action phase source-law candidate"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by research content patterns"},
        {"obligation": "construct_phase_periodicity_source_law", "satisfied": any(row["minimal_positive_representative"] for row in phases), "detail": "minimal internal phase-area section 2*pi/alpha_geo was constructed"},
        {"obligation": "test_action_length_time_scale_orbits", "satisfied": len(orbits) == len(TARGETS) * len(SCALE_MULTIPLIERS), "detail": "phase section was tested against action/length/time scale targets"},
        {"obligation": "export_dimensionful_physical_unit", "satisfied": False, "detail": "0 candidates export physical action/length/time calibration without standard-physics import"},
    ]
    return {
        "status": "P3111_INTERNAL_PHASE_AREA_SECTION_POSITIVE_PHYSICAL_UNIT_NO_GO",
        "input_hashes": {"P3110": hashlib.sha256(P3110.read_bytes()).hexdigest() if P3110.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "audit_object": {
                "object": "AlphaGeoSymplecticActionPhaseSourceLawAudit",
                "ontology": "the nadsoliton is self-coupled primordial information; alpha_geo is a four-bit Shannon cell, and internal phase periodicity can define a dimensionless phase-area section before any physical unit calibration",
                "source_law_tested": "A_phi = 2*pi/alpha_geo as the minimal positive representative of alpha_geo*A in 2*pi*Z",
            },
            "phase_periodicity_rows": phases,
            "scale_orbit_rows": orbits,
            "candidate_source_law_rows": candidates,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3110_accepted_dimension_comparison_standards": p3110.get("finite_certificate", {}).get("accepted_dimension_comparison_standards"),
            "phase_periodicity_rows": len(phases),
            "minimal_positive_phase_sections": sum(row["minimal_positive_representative"] for row in phases),
            "scale_orbit_rows": len(orbits),
            "candidate_source_laws": len(candidates),
            "candidate_gate_rows": len(gates),
            "accepted_physical_unit_source_laws": len(accepted_physical),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3111 makes real internal progress: alpha_geo and phase periodicity export a unique minimal positive internal phase-area section A_phi=2*pi/alpha_geo.  This is not yet a physical unit source: it is dimensionless internal phase accounting and does not calibrate action, length, or time without importing hbar/Planck or external rods/clocks/light/apparatus.",
            "negative_export_flags": {
                key: False
                for key in [
                    "dimensionful_action_unit_exported",
                    "physical_length_time_unit_exported",
                    "detector_calibration_exported",
                    "spacetime_eom_exported",
                    "QW_2191_discharged",
                    "strict_selector_closure_exported",
                    "ltotal_exported",
                    "bridge_role_transfer_exported",
                    "toe_closure_exported",
                ]
            },
            "positive_scoped_flags": {
                "alpha_geo_four_bit_entropy_confirmed": True,
                "minimal_internal_phase_area_section_exported": True,
                "symplectic_phase_area_candidate_strengthened": True,
                "imported_bohr_sommerfeld_hbar_law_rejected": True,
            },
            "next_honest_step": "Construct exactly one internal calibration functional C_phi that maps the dimensionless phase-area section A_phi=2*pi/alpha_geo to a dimensionful action comparison without hbar/Planck, rods, clocks, observed light, apparatus, selector replay, L_total, bridge/role-transfer, or ToE.  It must specify its codomain, prove scale covariance is broken internally, and show how action calibration would induce length/time calibration; otherwise preserve the P3105-P3111 physical-unit no-go.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3111/S2061 symplectic/action phase source-law audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite certificate",
        f"- P3110 accepted dimension comparison standards: `{cert['p3110_accepted_dimension_comparison_standards']}`",
        f"- content grep lanes: `{cert['content_grep_lanes']}`",
        f"- phase periodicity rows: `{cert['phase_periodicity_rows']}`",
        f"- minimal positive phase sections: `{cert['minimal_positive_phase_sections']}`",
        f"- scale-orbit rows: `{cert['scale_orbit_rows']}`",
        f"- candidate source laws: `{cert['candidate_source_laws']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted physical unit source laws: `{cert['accepted_physical_unit_source_laws']}`",
        f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`",
        "",
        "## Decision",
        payload["decision"]["bounded_result"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P3111/S2061 symplectic/action phase source-law audit",
        "## P3111/S2061 symplectic/action phase source-law audit\n\n`P3111/S2061` executes the P3110-recommended source-law audit for `symplectic_phase_area_candidate`.  It constructs `8` phase-periodicity rows, `12` action/length/time scale-orbit rows, `3` candidate source laws, and a `3 x 8 = 24` gate matrix.  The positive scoped result is that `A_phi = 2*pi/alpha_geo` is a unique minimal positive internal phase-area section for `alpha_geo*A in 2*pi*Z`.  The bounded obstruction is that this section is still dimensionless internal phase accounting and exports no physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P3111/S2061 phase-area section remains nonphysical",
        "## P3111/S2061 phase-area section remains nonphysical\n\n`P3111/S2061` strengthens the P3110 symplectic shape by deriving the internal phase-area section `A_phi=2*pi/alpha_geo`, but it remains dimensionless until an internal calibration functional maps it to a dimensionful action comparison.  Therefore it is not yet a physical Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n",
    )
    append_once(
        AGENTS,
        "Current symplectic/action phase source-law guardrail (P3111/S2061, 2026-06-26)",
        "## Current symplectic/action phase source-law guardrail (P3111/S2061, 2026-06-26)\n\n- P3111 tests the P3110-selected `symplectic_phase_area_candidate` and derives the unique minimal positive internal phase-area section `A_phi=2*pi/alpha_geo` from `alpha_geo*A in 2*pi*Z`.\n- This is real internal progress, but it remains dimensionless phase accounting: `0` candidates export a physical action/length/time unit or calibration coupling without importing `hbar/Planck`, rods, clocks, observed light, or apparatus.\n- Do not promote `A_phi`, winding families, or imported Bohr-Sommerfeld/Planck calibration to detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one internal calibration functional `C_phi` mapping `A_phi` to a dimensionful action comparison while proving internal scale-covariance breaking; otherwise preserve the P3105-P3111 physical-unit no-go.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
