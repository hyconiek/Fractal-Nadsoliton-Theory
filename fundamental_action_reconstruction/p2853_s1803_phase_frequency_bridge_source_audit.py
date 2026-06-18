#!/usr/bin/env python3
"""P2853/S1803: phase/frequency bridge-source audit.

P2852 identifies the cleanest non-replay bridge obligation: the phase/frequency
row for omega/phi/topological data.  This packet makes that row explicit.  It
separates a positive finite transport witness from the stronger missing theorem:
a strict, non-premise source law for the selected phase/frequency data.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2852 = GEN / "p2852_s1802_kernel_bridge_obligation_reconciliation_matrix.json"
OUT = GEN / "p2853_s1803_phase_frequency_bridge_source_audit.json"
MD = GEN / "p2853_s1803_phase_frequency_bridge_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LEGACY_OMEGA = math.pi / 4.0
LEGACY_PHI = math.pi / 6.0
STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
DOMAIN = tuple(range(12))
Z12_UNITS = (1, 5, 7, 11)
TOL = 1e-12


def theta_legacy(x_value: float) -> float:
    return LEGACY_OMEGA * x_value + LEGACY_PHI


def theta_strict(d: int) -> float:
    return float(STRICT_OMEGA) * d + float(STRICT_PHI)


def affine_legacy_coordinate(d: int) -> float:
    return (theta_strict(d) - LEGACY_PHI) / LEGACY_OMEGA


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def bit_from_sign(value: int) -> int:
    if value == 1:
        return 0
    if value == -1:
        return 1
    raise ValueError("zero sign cannot be converted to a phase bit")


def best_scalar_fit(source: list[float], target: list[float]) -> dict[str, float]:
    denominator = sum(value * value for value in source)
    scalar = sum(value * wanted for value, wanted in zip(source, target)) / denominator
    residuals = [scalar * value - wanted for value, wanted in zip(source, target)]
    return {
        "best_scalar": scalar,
        "max_abs_residual": max(abs(value) for value in residuals),
        "l2_residual": math.sqrt(sum(value * value for value in residuals)),
    }


def phase_rows() -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN:
        legacy_at_integer = math.cos(theta_legacy(float(d)))
        strict = math.cos(theta_strict(d))
        x_affine = affine_legacy_coordinate(d)
        legacy_at_affine = math.cos(theta_legacy(x_affine))
        phase_factor = strict / legacy_at_integer
        phase_sign = sign(phase_factor)
        rows.append(
            {
                "d": d,
                "legacy_argument_at_integer_d": theta_legacy(float(d)),
                "strict_argument_at_d": theta_strict(d),
                "affine_legacy_coordinate_x_d": x_affine,
                "legacy_cos_at_integer_d": legacy_at_integer,
                "strict_cos_at_d": strict,
                "legacy_cos_at_affine_x_d": legacy_at_affine,
                "affine_transport_residual": legacy_at_affine - strict,
                "phase_frequency_transport_factor": phase_factor,
                "phase_factor_sign": phase_sign,
                "phase_factor_bit": bit_from_sign(phase_sign),
            }
        )
    return rows


def z12_automorphism_rows(legacy_signs: list[int], strict_signs: list[int]) -> list[dict[str, Any]]:
    rows = []
    for unit in Z12_UNITS:
        for offset in DOMAIN:
            mapped = [legacy_signs[(unit * d + offset) % 12] for d in DOMAIN]
            mismatches = [d for d, got, wanted in zip(DOMAIN, mapped, strict_signs) if got != wanted]
            rows.append(
                {
                    "unit": unit,
                    "offset": offset,
                    "mismatch_positions_against_strict_sign": mismatches,
                    "mismatch_count": len(mismatches),
                    "matches_strict_sign_pattern": len(mismatches) == 0,
                }
            )
    return rows


def edge_bit_rows(bits: list[int]) -> list[dict[str, Any]]:
    rows = []
    for left, right in zip(DOMAIN[:-1], DOMAIN[1:]):
        bit = bits[left] ^ bits[right]
        rows.append({"edge": f"{left}->{right}", "left_bit": bits[left], "right_bit": bits[right], "edge_flip_bit": bit})
    return rows


def phase_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    legacy_cos = [row["legacy_cos_at_integer_d"] for row in rows]
    strict_cos = [row["strict_cos_at_d"] for row in rows]
    legacy_signs = [sign(value) for value in legacy_cos]
    strict_signs = [sign(value) for value in strict_cos]
    bits = [row["phase_factor_bit"] for row in rows]
    automorphisms = z12_automorphism_rows(legacy_signs, strict_signs)
    scalar_fit = best_scalar_fit(legacy_cos, strict_cos)
    x_values = [row["affine_legacy_coordinate_x_d"] for row in rows]
    distances_to_integer = [abs(x_value - round(x_value)) for x_value in x_values]
    strict_minus_legacy_arguments = [theta_strict(d) - theta_legacy(float(d)) for d in DOMAIN]
    return {
        "continuous_affine_phase_transport_exact": max(abs(row["affine_transport_residual"]) for row in rows) <= TOL,
        "max_abs_affine_transport_residual": max(abs(row["affine_transport_residual"]) for row in rows),
        "affine_slope_omega_s_over_omega_l": float(STRICT_OMEGA) / LEGACY_OMEGA,
        "affine_intercept": (float(STRICT_PHI) - LEGACY_PHI) / LEGACY_OMEGA,
        "affine_coordinates_not_all_integers": any(distance > TOL for distance in distances_to_integer),
        "minimum_distance_to_integer_affine_coordinate": min(distances_to_integer),
        "maximum_distance_to_integer_affine_coordinate": max(distances_to_integer),
        "z12_unit_offset_automorphism_rows": automorphisms,
        "z12_unit_offset_automorphism_count_checked": len(automorphisms),
        "no_z12_unit_offset_reindex_matches_strict_sign_pattern": not any(row["matches_strict_sign_pattern"] for row in automorphisms),
        "best_z12_unit_offset_mismatch_count": min(row["mismatch_count"] for row in automorphisms),
        "legacy_sign_pattern": legacy_signs,
        "strict_sign_pattern": strict_signs,
        "phase_factor_bits": bits,
        "phase_edge_bit_rows": edge_bit_rows(bits),
        "phase_factor_bits_nonconstant": len(set(bits)) > 1,
        "scalar_phase_replacement_fails": scalar_fit["max_abs_residual"] > 1e-6,
        "scalar_phase_best_fit": scalar_fit,
        "same_integer_coordinate_phase_identity_fails": abs(LEGACY_OMEGA - float(STRICT_OMEGA)) > TOL
        or abs(LEGACY_PHI - float(STRICT_PHI)) > TOL,
        "constant_phase_shift_only_fails": max(strict_minus_legacy_arguments) - min(strict_minus_legacy_arguments) > TOL,
        "strict_omega_fraction": {"numerator": STRICT_OMEGA.numerator, "denominator": STRICT_OMEGA.denominator},
        "strict_phi_fraction": {"numerator": STRICT_PHI.numerator, "denominator": STRICT_PHI.denominator},
    }


def source_candidate_matrix(summary: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "same_integer_coordinate_phase_identity",
            "finite_witness_passes": not summary["same_integer_coordinate_phase_identity_fails"],
            "exports_strict_source_law": False,
            "verdict": "blocked: legacy and strict omega/phi are not identical on the same integer coordinate.",
        },
        {
            "candidate": "constant_phase_shift_only",
            "finite_witness_passes": not summary["constant_phase_shift_only_fails"],
            "exports_strict_source_law": False,
            "verdict": "blocked: phase differences vary with d, so a shift-only passage cannot source the strict frequency.",
        },
        {
            "candidate": "z12_unit_offset_reindexing",
            "finite_witness_passes": not summary["no_z12_unit_offset_reindex_matches_strict_sign_pattern"],
            "exports_strict_source_law": False,
            "verdict": "blocked: exhaustive Aut(Z12) unit+offset reindexing leaves residual sign-pattern mismatches.",
        },
        {
            "candidate": "continuous_affine_phase_coordinate_transport",
            "finite_witness_passes": summary["continuous_affine_phase_transport_exact"],
            "exports_strict_source_law": False,
            "verdict": "positive finite transport witness, but not a discrete topological automorphism or a strict source of omega/phi.",
        },
        {
            "candidate": "phase_sign_gf2_bookkeeping",
            "finite_witness_passes": summary["phase_factor_bits_nonconstant"],
            "exports_strict_source_law": False,
            "verdict": "positive sign/coboundary bookkeeping only; it records consequences of chosen phase parameters, not their source.",
        },
    ]


def professorial_analysis() -> dict[str, Any]:
    return {
        "foundation_diagnosis": [
            "The program has strong finite witnesses and careful provenance discipline, but the decisive missing object is sourcehood: where selected numerical/topological data come from.",
            "Kernel completion must be decomposed into amplitude, damping/compression, phase/frequency, selector/topological data, completion semantics, and then role transfer; skipping any layer creates a false theorem.",
            "Finite transport identities are useful bridge evidence only when kept below the source-law threshold.",
            "A role-bearing L_total cannot be assembled from kernel similarity; it needs unit-bearing localization, coupling coefficients, and variational chain rules after bridge/source closure.",
        ],
        "research_path_to_closure": [
            "Export one strict phase/frequency source theorem: an internal law selecting omega=743/4000 and phi=13/80, or an equivalent topological phase datum, without imported coordinate convention.",
            "Re-run the bridge-obligation matrix with the new source as a typed atom, checking whether it reduces full-bridge missing premises without touching role transfer.",
            "Only after phase, damping, and amplitude atoms are source-level, build the completion-map theorem and prove residual-zero transport on the audited domain plus stated analytic extension conditions.",
            "Run a separate role-transfer audit for legacy physical claims; do not fold it into the bridge theorem.",
            "After bridge and role-transfer boundaries are explicit, revisit L_total/EOM/Hamiltonian with unit-bearing source density, localization/pullback, coupling coefficient, and variational derivative requirements.",
        ],
    }


def build_payload(p2852: dict[str, Any]) -> dict[str, Any]:
    rows = phase_rows()
    summary = phase_summary(rows)
    candidates = source_candidate_matrix(summary)
    accepted_count = sum(1 for row in candidates if row["exports_strict_source_law"])
    facts = {
        "p2852_rechecked": p2852.get("status") == "P2852_KERNEL_BRIDGE_OBLIGATION_RECONCILIATION_MATRIX_NO_CLOSURE",
        "continuous_affine_transport_witness_positive": summary["continuous_affine_phase_transport_exact"],
        "z12_automorphism_source_rejected": summary["no_z12_unit_offset_reindex_matches_strict_sign_pattern"],
        "scalar_phase_replacement_rejected": summary["scalar_phase_replacement_fails"],
        "same_coordinate_identity_rejected": summary["same_integer_coordinate_phase_identity_fails"],
        "no_candidate_exports_strict_phase_frequency_source": accepted_count == 0,
    }
    return {
        "status": "P2853_PHASE_FREQUENCY_BRIDGE_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2852": sha(P2852)},
        "phase_frequency_bridge_source_audit": {
            "input_statuses_rechecked": {"P2852": p2852.get("status")},
            "parameters": {
                "domain": list(DOMAIN),
                "legacy_omega": "pi/4",
                "legacy_phi": "pi/6",
                "strict_omega": "743/4000",
                "strict_phi": "13/80",
                "z12_units_checked": list(Z12_UNITS),
            },
            "phase_rows": rows,
            "summary": summary,
            "source_candidate_matrix": candidates,
            "accepted_candidate_count": accepted_count,
            "professorial_analysis": professorial_analysis(),
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_phase_frequency_bridge_source_obstruction_audit": all(facts.values()),
            "exports_phase_frequency_transport_witness": summary["continuous_affine_phase_transport_exact"],
            "exports_phase_frequency_bridge_source_atom": False,
        },
        "decision": {
            "negative_export_flags": {
                "phase_frequency_bridge_source_atom_exported": False,
                "strict_omega_phi_source_law_exported": False,
                "selector_topological_source_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "selector_closure_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2853 tests exactly one remaining bridge-source atom: phase/frequency passage for omega/phi/topological data.  Continuous affine phase-coordinate transport is exact on the audited Z12 domain, and the phase-bit profile is a real finite witness.  But the affine map is not a Z12 automorphism, scalar phase replacement fails, same-coordinate identity fails, and no audited candidate exports a non-premise strict source law for omega/phi.  Therefore this is a witness/obstruction audit, not a phase/frequency bridge-source closure.",
            "next_honest_step": "Do not replay phase-sign bookkeeping, affine transport, damping, amplitude, EML syntax, density normalizers, role transfer, L_total, EOM, Hamiltonian, or ToE promotion.  The next proof-grade move must supply one genuine strict phase/frequency source law for omega/phi/topological data, or a genuinely new eta/beta source law.  Without such a new typed source premise, preserve a no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["phase_frequency_bridge_source_audit"]
    summary = audit["summary"]
    analysis = audit["professorial_analysis"]
    lines = [
        "# P2853/S1803 phase/frequency bridge-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Professorial foundation analysis",
    ]
    lines.extend(f"- {item}" for item in analysis["foundation_diagnosis"])
    lines.extend(
        [
            "",
            "## Computation summary",
            f"- continuous_affine_phase_transport_exact={summary['continuous_affine_phase_transport_exact']}",
            f"- max_abs_affine_transport_residual={summary['max_abs_affine_transport_residual']}",
            f"- no_z12_unit_offset_reindex_matches_strict_sign_pattern={summary['no_z12_unit_offset_reindex_matches_strict_sign_pattern']}",
            f"- best_z12_unit_offset_mismatch_count={summary['best_z12_unit_offset_mismatch_count']}",
            f"- scalar_phase_replacement_fails={summary['scalar_phase_replacement_fails']}",
            f"- scalar_phase_best_fit_max_abs_residual={summary['scalar_phase_best_fit']['max_abs_residual']}",
            f"- phase_factor_bits={summary['phase_factor_bits']}",
            "",
            "## Source candidate matrix",
        ]
    )
    for row in audit["source_candidate_matrix"]:
        lines.append(
            f"- {row['candidate']}: finite_witness_passes={row['finite_witness_passes']}; "
            f"exports_strict_source_law={row['exports_strict_source_law']}; {row['verdict']}"
        )
    lines.extend(["", "## Research path"])
    lines.extend(f"- {item}" for item in analysis["research_path_to_closure"])
    lines.extend(
        [
            "",
            "## Boundary",
            payload["decision"]["reason"],
            "",
            "## Recommendation",
            payload["decision"]["next_honest_step"],
        ]
    )
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2852))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2853/S1803 phase/frequency bridge-source audit",
        "## P2853/S1803 phase/frequency bridge-source audit\n\n"
        "`P2853/S1803` audits the phase/frequency bridge-source atom for `omega/phi` and topological data.  "
        "It finds an exact continuous affine phase-coordinate transport on the audited `Z12` domain and a real finite phase-bit witness, but rejects same-coordinate identity, scalar phase replacement, and exhaustive `Aut(Z12)` unit+offset reindexing as strict source laws.  "
        "No strict `omega/phi` source law, selector/topological source, full kernel bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2853/S1803 phase-frequency `L_total` guard",
        "## P2853/S1803 phase-frequency `L_total` guard\n\n"
        "`P2853/S1803` contributes no action term.  Finite affine phase transport and phase-bit bookkeeping do not provide a unit-bearing source density, coupling coefficient, localization/pullback map, variational chain rule, nonproxy `L_total` insertion, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current phase/frequency bridge-source audit guardrail (P2853/S1803, 2026-06-18)",
        "## Current phase/frequency bridge-source audit guardrail (P2853/S1803, 2026-06-18)\n\n"
        "- P2853 tests the phase/frequency bridge-source atom for `omega/phi` and topological data after P2852.\n"
        "- Continuous affine phase-coordinate transport is exact on the audited `Z12` domain and the phase-bit profile is a real finite witness, but the map is not a `Z12` automorphism, scalar phase replacement fails, same-coordinate identity fails, and no non-premise strict source law for `omega/phi` is exported.\n"
        "- Do not promote phase-sign bookkeeping, affine phase transport, or topological witness language to full kernel bridge, strict selector/topological source, `QW-2191` discharge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a genuine strict phase/frequency source law for `omega/phi`/topological data, a genuinely new `eta/beta` source law, or else preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
