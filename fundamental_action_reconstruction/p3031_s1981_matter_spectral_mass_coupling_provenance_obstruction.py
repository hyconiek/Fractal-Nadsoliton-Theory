#!/usr/bin/env python3
"""P3031/S1981: matter spectral mass/coupling provenance obstruction.

Attack exactly one remaining P3029 matter atom after P3030: can the sorted DFT
magnitude observable be promoted to mass/coupling provenance?  The finite answer
is no on current artifacts: raw spectral quantities scale with kernel amplitude,
normalized quantities are dimensionless, and P3030 still blocks field/sector
localization.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3029_s1979_matter_spectral_observable_generator_obstruction import dft_magnitudes, rounded_sorted
from p3030_s1980_matter_spectral_field_localizer_obstruction import OUT as P3030, dft

OUT = GEN / "p3031_s1981_matter_spectral_mass_coupling_provenance_obstruction.json"
MD = GEN / "p3031_s1981_matter_spectral_mass_coupling_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

RESCALE_C = 3.0


def kernel_vector(scale: float = 1.0) -> list[float]:
    return [scale * k_strict(label) for label in range(1, N + 1)]


def unsorted_magnitudes(values: list[float]) -> list[float]:
    return [abs(z) for z in dft(values)]


def spectral_moments(magnitudes: list[float]) -> dict[str, float]:
    total = sum(magnitudes)
    if total == 0:
        return {"mean_bin": 0.0, "variance_bin": 0.0, "rms_magnitude": 0.0}
    mean = sum(i * v for i, v in enumerate(magnitudes)) / total
    variance = sum(((i - mean) ** 2) * v for i, v in enumerate(magnitudes)) / total
    rms = math.sqrt(sum(v * v for v in magnitudes) / len(magnitudes))
    return {"mean_bin": round(mean, 12), "variance_bin": round(variance, 12), "rms_magnitude": round(rms, 12)}


def safe_ratios(values: list[float]) -> list[float]:
    ratios = []
    for a, b in zip(values, values[1:] + values[:1]):
        ratios.append(round(a / b, 12) if abs(b) > 1e-15 else float("inf"))
    return ratios


def scaling_exponent(base: float, scaled: float, factor: float = RESCALE_C) -> float | None:
    if base == 0 or scaled == 0 or factor <= 0:
        return None
    return round(math.log(abs(scaled / base), factor), 12)


def build_matrix() -> dict[str, Any]:
    base_values = kernel_vector(1.0)
    scaled_values = kernel_vector(RESCALE_C)
    base_mags = unsorted_magnitudes(base_values)
    scaled_mags = unsorted_magnitudes(scaled_values)
    base_signature = rounded_sorted(base_mags)
    scaled_signature = rounded_sorted(scaled_mags)
    normalized_base = [round(v / sum(base_signature), 12) for v in base_signature]
    normalized_scaled = [round(v / sum(scaled_signature), 12) for v in scaled_signature]
    base_moments = spectral_moments(base_mags)
    scaled_moments = spectral_moments(scaled_mags)

    candidate_rows = [
        {
            "candidate": "raw_sorted_magnitude_as_mass_proxy",
            "rescaling_exponent": scaling_exponent(base_signature[-1], scaled_signature[-1]),
            "dimensionless_or_unitless": False,
            "accepted_as_mass_coupling_provenance": False,
            "failure": "raw magnitudes scale with arbitrary K amplitude and have no physical mass unit theorem",
        },
        {
            "candidate": "normalized_spectral_shape_as_coupling_proxy",
            "rescaling_exponent": scaling_exponent(normalized_base[-1], normalized_scaled[-1]),
            "dimensionless_or_unitless": True,
            "accepted_as_mass_coupling_provenance": False,
            "failure": "normalization removes amplitude but leaves a dimensionless shape without mass/coupling units or target map",
        },
        {
            "candidate": "spectral_moment_mass_proxy",
            "rescaling_exponent": scaling_exponent(base_moments["rms_magnitude"], scaled_moments["rms_magnitude"]),
            "dimensionless_or_unitless": False,
            "accepted_as_mass_coupling_provenance": False,
            "failure": "RMS magnitude scales with amplitude; bin moments are index statistics, not physical mass provenance",
        },
        {
            "candidate": "adjacent_magnitude_ratio_coupling_proxy",
            "rescaling_exponent": scaling_exponent(safe_ratios(base_signature)[0], safe_ratios(scaled_signature)[0]),
            "dimensionless_or_unitless": True,
            "accepted_as_mass_coupling_provenance": False,
            "failure": "ratios are scale-invariant but dimensionless and inherit P3030's no field/sector localizer",
        },
    ]
    obligations = [
        {"obligation": "attacks_single_P3029_missing_atom", "satisfied": True, "detail": "only mass/coupling provenance for the P3029 spectral generator is tested"},
        {"obligation": "uses_existing_P3029_generator", "satisfied": True, "detail": "input is the DFT magnitude signature of sampled K_strict_gate"},
        {"obligation": "field_representation_localizer_available", "satisfied": False, "detail": "P3030 exports no field/sector localizer for this magnitude observable"},
        {"obligation": "absolute_mass_or_coupling_unit_source", "satisfied": False, "detail": "raw candidates scale under K -> cK and normalized candidates are unitless"},
        {"obligation": "target_independent_mass_coupling_map", "satisfied": False, "detail": "no theorem maps signature components to named matter masses/couplings"},
        {"obligation": "unit_bearing_action_eom_insertion", "satisfied": False, "detail": "no unit-bearing action/EOM/Hamiltonian insertion is exported"},
    ]
    return {
        "object": "MatterSpectralMassCouplingProvenance_RescalingNormalizationObstructionMatrix",
        "tested_missing_atom": "mass/coupling provenance theorem for P3029 matter spectral observable",
        "rescaling_test": {"K_scale_factor": RESCALE_C, "raw_signature_l1_scale_ratio": round(sum(scaled_signature) / sum(base_signature), 12), "normalized_signature_preserved": normalized_base == normalized_scaled},
        "candidate_rows": candidate_rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "candidate_rows": len(candidate_rows),
            "accepted_candidates": sum(1 for row in candidate_rows if row["accepted_as_mass_coupling_provenance"]),
            "raw_scale_covariant_rows": sum(1 for row in candidate_rows if row["rescaling_exponent"] is not None and abs(row["rescaling_exponent"] - 1.0) <= 1e-9),
            "dimensionless_unitless_rows": sum(1 for row in candidate_rows if row["dimensionless_or_unitless"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "mass_coupling_provenance_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_mass_coupling_provenance"] for row in candidate_rows),
        },
    }


def build_payload(p3030_path: Any) -> dict[str, Any]:
    read_json(p3030_path)
    matrix = build_matrix()
    return {
        "status": "P3031_MATTER_SPECTRAL_MASS_COUPLING_PROVENANCE_OBSTRUCTION_NO_MATTER_EXPORT",
        "input_hashes": {"P3030": hashlib.sha256(p3030_path.read_bytes()).hexdigest() if p3030_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The P3029/P3030 spectral magnitude observable does not currently export mass/coupling provenance.  Raw spectral mass proxies scale with arbitrary kernel amplitude, normalized shape/ratio proxies are dimensionless, P3030 supplies no field/sector localizer, and no target-independent unit-bearing mass/coupling map or action/EOM insertion theorem is exported.",
            "negative_export_flags": {k: False for k in ["mass_coupling_provenance_exported", "matter_sector_exported", "field_representation_localizer_exported", "selector_source_exported", "unit_bearing_action_eom_source_exported", "energy_hamiltonian_exported", "observed_physics_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay spectral mass/coupling proxies for the P3029 magnitude generator.  The honest next move is a P3029-P3031 matter spectral lane reconciliation/no-new-live-frontier certificate, or a pivot to a different P3028 foundation atom with a genuinely new strict typed object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3031/S1981 matter spectral mass/coupling provenance obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- candidate rows: `{c['candidate_rows']}`",
        f"- accepted candidates: `{c['accepted_candidates']}`",
        f"- raw scale-covariant rows: `{c['raw_scale_covariant_rows']}`",
        f"- dimensionless/unitless rows: `{c['dimensionless_unitless_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- mass/coupling provenance exported: `{c['mass_coupling_provenance_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(P3030)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3031/S1981 matter spectral mass/coupling provenance obstruction", "## P3031/S1981 matter spectral mass/coupling provenance obstruction\n\n`P3031/S1981` attacks exactly one remaining P3029/P3030 matter atom: mass/coupling provenance for the spectral magnitude observable.  It tests raw sorted magnitudes, normalized spectral shape, spectral moments, and adjacent magnitude ratios under `K -> 3K`.  Raw mass proxies scale with arbitrary kernel amplitude, normalized shape/ratio proxies are dimensionless, P3030 still exports no field/sector localizer, and `0/4` candidate rows export a target-independent unit-bearing mass/coupling theorem.  Therefore no matter-sector export, observed physics, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3031/S1981 matter spectral mass/coupling `L_total` guard", "## P3031/S1981 matter spectral mass/coupling `L_total` guard\n\n`P3031/S1981` adds no physical `L_total` term.  The spectral mass/coupling proxy family either scales with arbitrary kernel amplitude or becomes dimensionless after normalization, while the P3030 field/sector localizer remains absent.  No unit-bearing action/EOM/Hamiltonian insertion is exported.\n")
    append_once(AGENTS, "Current matter spectral mass/coupling provenance guardrail (P3031/S1981, 2026-06-22)", "## Current matter spectral mass/coupling provenance guardrail (P3031/S1981, 2026-06-22)\n\n- P3031 attacks exactly one remaining P3029/P3030 matter atom: mass/coupling provenance for the spectral magnitude observable.\n- The finite `K -> 3K` rescaling/normalization test finds `0/4` accepted mass/coupling provenance candidates: raw spectral proxies scale with arbitrary amplitude, normalized ratios are dimensionless, and P3030 still blocks field/sector localization.\n- Do not promote raw magnitudes, normalized spectral shapes, spectral moments, or adjacent magnitude ratios to matter sector, observed physics, unit-bearing action/EOM, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- The next honest move is a P3029-P3031 matter spectral lane reconciliation/no-new-live-frontier certificate, or a pivot to a different P3028 foundation atom with a genuinely new strict typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
