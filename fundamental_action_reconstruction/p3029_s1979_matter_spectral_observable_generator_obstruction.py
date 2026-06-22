#!/usr/bin/env python3
"""P3029/S1979: matter-row spectral observable generator obstruction.

Attack exactly one P3028 foundation atom: an observer-independent observable
generator for one classical row.  We choose the matter_fields row and construct a
strict-kernel spectral-magnitude generator with explicit domain/codomain.
"""
from __future__ import annotations

import cmath, hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, UNITS, k_strict
from p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice import OUT as P3028

OUT = GEN / "p3029_s1979_matter_spectral_observable_generator_obstruction.json"
MD = GEN / "p3029_s1979_matter_spectral_observable_generator_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def kernel_vector(unit: int = 1) -> list[float]:
    # Unit relabeling on Z12 labels; index j=0 represents label 1.
    values = [0.0] * N
    for label in range(1, N + 1):
        image_label = ((unit * label - 1) % N) + 1
        values[image_label - 1] = k_strict(label)
    return values


def dft_magnitudes(values: list[float]) -> list[float]:
    mags = []
    for m in range(N):
        total = sum(values[j] * cmath.exp(-2j * math.pi * m * j / N) for j in range(N))
        mags.append(abs(total))
    return mags


def rounded_sorted(values: list[float]) -> list[float]:
    return sorted(round(v, 12) for v in values)


def build_generator_matrix() -> dict[str, Any]:
    base_values = kernel_vector()
    base_signature = rounded_sorted(dft_magnitudes(base_values))
    equivariance_rows = []
    for unit in UNITS:
        signature = rounded_sorted(dft_magnitudes(kernel_vector(unit)))
        equivariance_rows.append({
            "unit": unit,
            "sorted_spectral_signature_preserved": signature == base_signature,
            "signature_l1_delta": round(sum(abs(a - b) for a, b in zip(signature, base_signature)), 15),
        })
    obligations = [
        {"obligation": "attacks_single_P3028_foundation_atom", "satisfied": True, "detail": "only observer-independent observable generator for matter_fields is tested"},
        {"obligation": "explicit_domain_codomain", "satisfied": True, "detail": "domain: sampled strict-kernel vector on Z12; codomain: sorted DFT magnitude signature in R^12"},
        {"obligation": "observer_independent_formula", "satisfied": True, "detail": "computed directly from K_strict_gate values, no observer input"},
        {"obligation": "U12_unit_invariant_signature", "satisfied": all(row["sorted_spectral_signature_preserved"] for row in equivariance_rows), "detail": "sorted magnitude spectrum is preserved under U(12) relabeling"},
        {"obligation": "field_representation_localizer", "satisfied": False, "detail": "spectral magnitudes do not assign modes to physical fields or particles"},
        {"obligation": "mass_coupling_provenance", "satisfied": False, "detail": "no mass/coupling scale or unit theorem maps the signature to matter parameters"},
        {"obligation": "selector_or_sector_source", "satisfied": False, "detail": "sorted magnitudes forget directed/sector labels and do not select a physical branch"},
    ]
    return {
        "object": "MatterSpectralObservableGenerator_DFTMagnitudeInvariantObstructionMatrix",
        "classical_row": "matter_fields",
        "domain": "K_strict_gate samples on Z12 labels 1..12",
        "codomain": "sorted DFT magnitude signature in R^12",
        "base_signature": base_signature,
        "unit_equivariance_rows": equivariance_rows,
        "proof_obligations": obligations,
        "accepted_as_observer_independent_observable_generator": all(row["satisfied"] for row in obligations[:4]),
        "accepted_as_matter_sector_export": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3028_path: Any) -> dict[str, Any]:
    read_json(p3028_path)
    matrix = build_generator_matrix()
    return {
        "status": "P3029_MATTER_SPECTRAL_OBSERVABLE_GENERATOR_OBSTRUCTION_NO_CLASSICAL_EXPORT",
        "input_hashes": {"P3028": hashlib.sha256(p3028_path.read_bytes()).hexdigest() if p3028_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "signature_length": len(matrix["base_signature"]),
            "unit_rows": len(matrix["unit_equivariance_rows"]),
            "unit_invariant_rows": sum(1 for row in matrix["unit_equivariance_rows"] if row["sorted_spectral_signature_preserved"]),
            "observer_independent_generator_accepted": matrix["accepted_as_observer_independent_observable_generator"],
            "matter_sector_export_accepted": matrix["accepted_as_matter_sector_export"],
        },
        "decision": {
            "breakthrough": "A real observer-independent observable generator for the P3028 matter_fields row was constructed: the sorted DFT magnitude signature of K_strict_gate on Z12.  It has explicit types and is invariant under U(12) relabeling, but it is not a matter-sector export because it lacks a field representation/localizer, mass/coupling provenance, and selector/sector source.",
            "negative_export_flags": {k: False for k in ["matter_sector_exported", "field_representation_exported", "mass_coupling_provenance_exported", "selector_source_exported", "unit_bearing_action_eom_source_exported", "energy_hamiltonian_exported", "observed_physics_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote spectral magnitude signatures to matter physics.  The next proof-grade move may attack exactly one missing matter atom for this generator: a field-representation/localizer theorem or a mass/coupling provenance theorem; otherwise return to the P3028 lattice and choose another single foundation atom.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3029/S1979 matter spectral observable generator obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- signature length: `{c['signature_length']}`",
        f"- U(12) rows: `{c['unit_rows']}`",
        f"- invariant rows: `{c['unit_invariant_rows']}`",
        f"- observer-independent generator accepted: `{c['observer_independent_generator_accepted']}`",
        f"- matter-sector export accepted: `{c['matter_sector_export_accepted']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(P3028)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3029/S1979 matter spectral observable generator obstruction", "## P3029/S1979 matter spectral observable generator obstruction\n\n`P3029/S1979` attacks exactly one P3028 foundation atom: an observer-independent observable generator for the `matter_fields` row.  It constructs an explicit typed generator from sampled `K_strict_gate` on `Z12` to the sorted DFT magnitude signature in `R^12`.  The finite positive is real: the formula is observer-independent and all `4/4` `U(12)` relabeling rows preserve the sorted spectral signature.  The bounded obstruction is that spectral magnitudes alone do not export a field representation/localizer, mass/coupling provenance, or selector/sector source.  No matter-sector export, observed physics, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3029/S1979 matter spectral observable `L_total` guard", "## P3029/S1979 matter spectral observable `L_total` guard\n\n`P3029/S1979` adds no physical `L_total` term.  The sorted DFT magnitude signature is a real observer-independent strict-kernel observable generator for the matter row, but without field-representation, mass/coupling, selector/sector, and unit-bearing action the signature remains pre-physical readout data.\n")
    append_once(AGENTS, "Current matter spectral observable generator guardrail (P3029/S1979, 2026-06-22)", "## Current matter spectral observable generator guardrail (P3029/S1979, 2026-06-22)\n\n- P3029 attacks one P3028 foundation atom: an observer-independent observable generator for the `matter_fields` row.\n- The sorted DFT magnitude signature of sampled `K_strict_gate` has explicit domain/codomain and is invariant under all `4/4` `U(12)` relabeling rows.\n- Do not promote spectral magnitude signatures to matter sector, observed physics, unit-bearing action/EOM, `L_total`, selector, bridge/role-transfer, or ToE closure without a field-representation/localizer plus mass/coupling provenance theorem.\n- The next honest move may attack exactly one missing matter atom for this generator, preferably field-representation/localizer or mass/coupling provenance.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
