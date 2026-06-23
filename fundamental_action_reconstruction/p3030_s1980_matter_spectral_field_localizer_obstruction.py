#!/usr/bin/env python3
"""P3030/S1980: matter spectral field-representation/localizer obstruction.

Attack exactly one P3029 missing matter atom: can the sorted DFT magnitude
observable be upgraded to a field-representation/localizer theorem?  The finite
answer is no on current artifacts because magnitude-only data quotients out
translation/reflection phase information needed to localize a field/sector.
"""
from __future__ import annotations

import cmath, hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, UNITS, k_strict
from p3029_s1979_matter_spectral_observable_generator_obstruction import OUT as P3029, dft_magnitudes, rounded_sorted

OUT = GEN / "p3030_s1980_matter_spectral_field_localizer_obstruction.json"
MD = GEN / "p3030_s1980_matter_spectral_field_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def base_vector() -> list[float]:
    return [k_strict(label) for label in range(1, N + 1)]


def translate(values: list[float], shift: int) -> list[float]:
    return [values[(j - shift) % N] for j in range(N)]


def reflect(values: list[float]) -> list[float]:
    return [values[(-j) % N] for j in range(N)]


def dft(values: list[float]) -> list[complex]:
    return [sum(values[j] * cmath.exp(-2j * math.pi * m * j / N) for j in range(N)) for m in range(N)]


def autocorrelation_from_magnitudes(magnitudes: list[float]) -> list[float]:
    # Wiener-Khinchin inverse of |DFT|^2.  It is translation-invariant and even.
    power = [m * m for m in magnitudes]
    out = []
    for lag in range(N):
        total = sum(power[m] * cmath.exp(2j * math.pi * m * lag / N) for m in range(N)) / N
        out.append(round(total.real, 12))
    return out


def zero_phase_reconstruction(magnitudes: list[float]) -> list[float]:
    # Best possible canonical reconstruction without phases: set all phases to zero.
    out = []
    for j in range(N):
        total = sum(magnitudes[m] * cmath.exp(2j * math.pi * m * j / N) for m in range(N)) / N
        out.append(round(total.real, 12))
    return out


def build_matrix() -> dict[str, Any]:
    values = base_vector()
    base_signature = rounded_sorted(dft_magnitudes(values))

    orbit_rows = []
    for shift in range(N):
        shifted = translate(values, shift)
        reflected_shifted = translate(reflect(values), shift)
        for name, vec in [("translation", shifted), ("reflection_translation", reflected_shifted)]:
            orbit_rows.append({
                "operation": name,
                "shift": shift,
                "sorted_signature_preserved": rounded_sorted(dft_magnitudes(vec)) == base_signature,
            })

    magnitudes = [abs(z) for z in dft(values)]
    autocorr = autocorrelation_from_magnitudes(magnitudes)
    zero_phase = zero_phase_reconstruction(magnitudes)
    localizer_rows = [
        {
            "candidate": "inverse_dft_with_zero_phase",
            "domain": "unsorted DFT magnitudes with all phases set to zero",
            "candidate_index": max(range(N), key=lambda i: zero_phase[i]) + 1,
            "strict_localizer_exported": False,
            "failure": "phase choice is a convention; translations of the source have the same magnitudes but different field sites",
        },
        {
            "candidate": "autocorrelation_peak",
            "domain": "inverse DFT of |DFT|^2",
            "candidate_lag": max(range(1, N), key=lambda i: autocorr[i]) if N > 1 else 0,
            "strict_localizer_exported": False,
            "failure": "autocorrelation is translation-invariant/even and supplies only relative lags, not an absolute field representation",
        },
        {
            "candidate": "largest_spectral_bin",
            "domain": "sorted or unsorted magnitude spectrum",
            "candidate_bin": max(range(N), key=lambda i: magnitudes[i]),
            "strict_localizer_exported": False,
            "failure": "spectral bin is not a matter field/particle localizer and carries no sector theorem or mass/coupling provenance",
        },
    ]
    obligations = [
        {"obligation": "attacks_single_P3029_missing_atom", "satisfied": True, "detail": "only field-representation/localizer is tested"},
        {"obligation": "uses_existing_P3029_generator", "satisfied": True, "detail": "input is the DFT magnitude signature of sampled K_strict_gate"},
        {"obligation": "translation_orbit_degeneracy_absent", "satisfied": False, "detail": "all translated rows preserve the magnitude signature"},
        {"obligation": "reflection_orbit_degeneracy_absent", "satisfied": False, "detail": "all reflected-translated rows preserve the magnitude signature"},
        {"obligation": "nonpremise_phase_recovery", "satisfied": False, "detail": "magnitude data omit Fourier phases; zero phase is a convention"},
        {"obligation": "absolute_field_site_or_sector_localized", "satisfied": False, "detail": "candidates produce conventional indices/lags/bins, not a strict field representation"},
    ]
    return {
        "object": "MatterSpectralFieldLocalizer_PhaseRetrievalOrbitObstructionMatrix",
        "tested_missing_atom": "field-representation/localizer theorem for P3029 matter spectral observable",
        "orbit_rows": orbit_rows,
        "localizer_candidate_rows": localizer_rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "dihedral_rows": len(orbit_rows),
            "signature_preserving_rows": sum(1 for row in orbit_rows if row["sorted_signature_preserved"]),
            "translation_rows": sum(1 for row in orbit_rows if row["operation"] == "translation"),
            "reflection_translation_rows": sum(1 for row in orbit_rows if row["operation"] == "reflection_translation"),
            "localizer_candidates": len(localizer_rows),
            "accepted_localizers": sum(1 for row in localizer_rows if row["strict_localizer_exported"]),
            "field_representation_localizer_exported": all(row["satisfied"] for row in obligations),
        },
    }


def build_payload(p3029_path: Any) -> dict[str, Any]:
    read_json(p3029_path)
    matrix = build_matrix()
    return {
        "status": "P3030_MATTER_SPECTRAL_FIELD_LOCALIZER_OBSTRUCTION_NO_MATTER_EXPORT",
        "input_hashes": {"P3029": hashlib.sha256(p3029_path.read_bytes()).hexdigest() if p3029_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The P3029 spectral magnitude observable cannot currently be upgraded to a field-representation/localizer theorem.  The full 24-row translation/reflection orbit preserves the sorted magnitude signature, Fourier phases are missing, and the tested inverse/lag/bin localizers remain conventional rather than strict matter-field localizers.",
            "negative_export_flags": {k: False for k in ["field_representation_localizer_exported", "matter_sector_exported", "mass_coupling_provenance_exported", "selector_source_exported", "unit_bearing_action_eom_source_exported", "energy_hamiltonian_exported", "observed_physics_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay phase-retrieval or orbit-localizer variants for the P3029 magnitude signature.  The next proof-grade move may attack the remaining mass/coupling provenance atom for this generator, or pivot back to the P3028 foundation lattice and choose a different single atom with a genuinely new strict typed object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3030/S1980 matter spectral field-localizer obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- dihedral orbit rows: `{c['dihedral_rows']}`",
        f"- signature-preserving rows: `{c['signature_preserving_rows']}`",
        f"- localizer candidates: `{c['localizer_candidates']}`",
        f"- accepted localizers: `{c['accepted_localizers']}`",
        f"- field-representation/localizer exported: `{c['field_representation_localizer_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(P3029)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3030/S1980 matter spectral field-localizer obstruction", "## P3030/S1980 matter spectral field-localizer obstruction\n\n`P3030/S1980` attacks exactly one P3029 missing matter atom: the field-representation/localizer theorem for the spectral magnitude observable.  It tests translation/reflection phase-retrieval degeneracy for sampled `K_strict_gate` on `Z12` and three candidate localizers: zero-phase inverse DFT, autocorrelation peak, and largest spectral bin.  All `24/24` dihedral orbit rows preserve the sorted magnitude signature, and `0/3` candidate localizers export a strict matter-field representation.  Therefore no matter-sector export, observed physics, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3030/S1980 matter spectral field-localizer `L_total` guard", "## P3030/S1980 matter spectral field-localizer `L_total` guard\n\n`P3030/S1980` adds no physical `L_total` term.  Magnitude-only spectral data remain phase-retrieval/dihedral-orbit ambiguous; the zero-phase inverse, autocorrelation-lag, and largest-bin anchors are localizer conventions rather than unit-bearing matter-field/EOM sources.\n")
    append_once(AGENTS, "Current matter spectral field-localizer guardrail (P3030/S1980, 2026-06-22)", "## Current matter spectral field-localizer guardrail (P3030/S1980, 2026-06-22)\n\n- P3030 attacks exactly one P3029 missing matter atom: field-representation/localizer for the sorted DFT magnitude observable.\n- The finite phase-retrieval/orbit test finds all `24/24` translation/reflection rows preserve the sorted magnitude signature and `0/3` tested localizer anchors export a strict field representation.\n- Do not promote zero-phase inverse DFT, autocorrelation peaks, largest spectral bins, or magnitude-orbit anchors to matter sector, observed physics, unit-bearing action/EOM, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- The next honest move may attack the remaining mass/coupling provenance atom for the P3029 generator, or pivot to a different P3028 foundation atom with a genuinely new strict typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
