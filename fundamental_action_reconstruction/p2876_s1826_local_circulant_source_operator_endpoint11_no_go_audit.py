#!/usr/bin/env python3
"""P2876/S1826: local circulant source-operator endpoint-11 no-go audit.

P2875 showed that two-dimensional D12 irrep/Fourier endpoint waves are global and
that delta_11 appears only after importing target phase coefficients.  This
packet tests the next concrete candidate: a local translation-compatible source
operator/density on Z12.

The finite audit enumerates all ternary local circulant stencils of radius two
(coefficients in {-1,0,1}) and the reflection-symmetric subfamily.  It applies
those operators to the intrinsic global Fourier inputs k=0..6 (constant,
irrep modes, and alternating parity).  A circulant local operator only rescales
or kills a Fourier input, so every nonzero output remains global and every zero
output is empty; no stencil computes singleton endpoint 11 without importing a
delta/phase seed.
"""
from __future__ import annotations

import cmath
import json
import math
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2875 = GEN / "p2875_s1825_d12_two_dimensional_irrep_endpoint_localization_no_go_audit.json"
OUT = GEN / "p2876_s1826_local_circulant_source_operator_endpoint11_no_go_audit.json"
MD = GEN / "p2876_s1826_local_circulant_source_operator_endpoint11_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
ENDPOINTS = tuple(range(MODULUS))
TARGET_ENDPOINT = 11
RADIUS = 2
OFFSETS = tuple(range(-RADIUS, RADIUS + 1))
COEFFICIENTS = (-1, 0, 1)
INPUT_MODES = tuple(range(0, 7))
TOL = 1e-9


def all_stencils() -> list[tuple[int, ...]]:
    return list(product(COEFFICIENTS, repeat=len(OFFSETS)))


def is_reflection_symmetric(stencil: tuple[int, ...]) -> bool:
    lookup = dict(zip(OFFSETS, stencil))
    return all(lookup[offset] == lookup[-offset] for offset in OFFSETS)


def fourier_input(mode: int) -> list[complex]:
    return [cmath.exp(2j * math.pi * mode * endpoint / MODULUS) for endpoint in ENDPOINTS]


def apply_circulant_stencil(stencil: tuple[int, ...], values: list[complex]) -> list[complex]:
    output: list[complex] = []
    for endpoint in ENDPOINTS:
        total = 0j
        for coefficient, offset in zip(stencil, OFFSETS):
            total += coefficient * values[(endpoint + offset) % MODULUS]
        output.append(total)
    return output


def support(values: list[complex]) -> list[int]:
    return [endpoint for endpoint, value in enumerate(values) if abs(value) > TOL]


def multiplier(stencil: tuple[int, ...], mode: int) -> complex:
    return sum(coefficient * cmath.exp(2j * math.pi * mode * offset / MODULUS) for coefficient, offset in zip(stencil, OFFSETS))


def summarize_stencil_family(stencils: list[tuple[int, ...]]) -> dict[str, Any]:
    singleton_hits: list[dict[str, Any]] = []
    nonzero_global_outputs = 0
    zero_outputs = 0
    partial_non_singleton_outputs = 0
    killed_by_mode = {str(mode): 0 for mode in INPUT_MODES}
    for stencil in stencils:
        for mode in INPUT_MODES:
            output = apply_circulant_stencil(stencil, fourier_input(mode))
            out_support = support(output)
            if out_support == [TARGET_ENDPOINT]:
                singleton_hits.append({"stencil": list(stencil), "mode": mode})
            if out_support == []:
                zero_outputs += 1
                killed_by_mode[str(mode)] += 1
            elif out_support == list(ENDPOINTS):
                nonzero_global_outputs += 1
            else:
                partial_non_singleton_outputs += 1
    return {
        "stencil_count": len(stencils),
        "input_mode_count": len(INPUT_MODES),
        "operator_input_pair_count": len(stencils) * len(INPUT_MODES),
        "singleton_11_hits": singleton_hits,
        "singleton_11_hit_count": len(singleton_hits),
        "nonzero_global_output_count": nonzero_global_outputs,
        "zero_output_count": zero_outputs,
        "partial_non_singleton_output_count": partial_non_singleton_outputs,
        "killed_by_mode": killed_by_mode,
    }


def sample_records(stencils: list[tuple[int, ...]]) -> list[dict[str, Any]]:
    samples = [(-1, 0, 1, 0, -1), (1, -1, 0, -1, 1), (0, 1, -1, 1, 0), (1, 0, 0, 0, -1)]
    records: list[dict[str, Any]] = []
    for stencil in samples:
        if stencil not in stencils:
            continue
        mode_supports = {str(mode): support(apply_circulant_stencil(stencil, fourier_input(mode))) for mode in INPUT_MODES}
        mode_multipliers = {
            str(mode): [round(multiplier(stencil, mode).real, 12), round(multiplier(stencil, mode).imag, 12)]
            for mode in INPUT_MODES
        }
        records.append(
            {
                "stencil_offsets": list(OFFSETS),
                "stencil": list(stencil),
                "reflection_symmetric": is_reflection_symmetric(stencil),
                "mode_supports": mode_supports,
                "mode_multipliers": mode_multipliers,
            }
        )
    return records


def build_payload(p2875: dict[str, Any]) -> dict[str, Any]:
    stencils = all_stencils()
    reflection_symmetric_stencils = [stencil for stencil in stencils if is_reflection_symmetric(stencil)]
    full_summary = summarize_stencil_family(stencils)
    reflection_summary = summarize_stencil_family(reflection_symmetric_stencils)
    facts = {
        "p2875_rechecked": p2875.get("status") == "P2875_D12_TWO_DIMENSIONAL_IRREP_ENDPOINT_LOCALIZATION_NO_GO_AUDIT_NO_CLOSURE",
        "all_ternary_radius2_circulant_stencils_enumerated": len(stencils) == 3**len(OFFSETS),
        "reflection_symmetric_subfamily_enumerated": len(reflection_symmetric_stencils) == 3 ** (RADIUS + 1),
        "input_modes_k0_through_k6_checked": INPUT_MODES == tuple(range(7)),
        "no_full_family_singleton_11_hits": full_summary["singleton_11_hit_count"] == 0,
        "no_reflection_symmetric_singleton_11_hits": reflection_summary["singleton_11_hit_count"] == 0,
        "outputs_are_only_zero_or_global_for_fourier_inputs": full_summary["partial_non_singleton_output_count"] == 0,
    }
    return {
        "status": "P2876_LOCAL_CIRCULANT_SOURCE_OPERATOR_ENDPOINT11_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2875": sha(P2875)},
        "local_circulant_source_operator_endpoint11_no_go_audit": {
            "input_status_rechecked": p2875.get("status"),
            "candidate_class": "ternary radius-2 local circulant source operators on Z12 applied to intrinsic Fourier inputs k=0..6",
            "modulus": MODULUS,
            "offsets": list(OFFSETS),
            "coefficient_alphabet": list(COEFFICIENTS),
            "input_modes": list(INPUT_MODES),
            "full_family_summary": full_summary,
            "reflection_symmetric_summary": reflection_summary,
            "sample_records": sample_records(stencils),
            "proof_certificate": {
                "operator_enumeration_step": "All 3^5 ternary radius-2 circulant stencils are enumerated, with the 3^3 reflection-symmetric subfamily separated.",
                "fourier_diagonalization_step": "Circular local operators are diagonal on Fourier inputs, so each output is a scalar multiple of the same global mode.",
                "support_step": "Every checked output is either zero or full-support on Z12; no output has singleton endpoint-11 support.",
                "sourcehood_step": "A local circulant operator can produce delta_11 only from an imported delta/phase-localized seed, not from intrinsic global inputs; no unit-bearing 9/5 coefficient theorem is exported.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_strict_local_endpoint_11_source_operator": False,
            "exports_boundary_source_law": False,
            "exports_orientation_source_law": False,
            "exports_chiral_endpoint_11_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "boundary_source_law_exported": False,
                "orientation_source_law_exported": False,
                "chiral_endpoint_11_source_law_exported": False,
                "selector_or_localizer_source_exported": False,
                "unit_bearing_coupling_localization_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2876 enumerates all ternary radius-2 local circulant source operators and their reflection-symmetric subfamily on Fourier inputs k=0..6.  Every output is zero or full-support; no stencil computes singleton endpoint 11 without an imported localized seed.  Therefore no strict local endpoint-11 source operator or unit-bearing 9/5 coupling theorem is exported.",
            "next_honest_step": "Do not replay local circulant stencils, reflection-symmetric local stencils, Fourier-diagonal operator action, D12 irrep waves, or full-DFT delta reconstruction as sourcehood.  A next proof-grade move must provide a non-circulant strict local defect/source with independently exported origin/support at endpoint 11 and a unit-bearing 9/5 coefficient theorem, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["local_circulant_source_operator_endpoint11_no_go_audit"]
    lines = [
        "# P2876/S1826 local circulant source-operator endpoint-11 no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Local operator audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- full stencil count: `{audit['full_family_summary']['stencil_count']}`",
        f"- reflection-symmetric stencil count: `{audit['reflection_symmetric_summary']['stencil_count']}`",
        f"- full singleton-11 hits: `{audit['full_family_summary']['singleton_11_hit_count']}`",
        f"- reflection-symmetric singleton-11 hits: `{audit['reflection_symmetric_summary']['singleton_11_hit_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2875))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2876/S1826 local circulant source-operator endpoint-11 no-go audit",
        "## P2876/S1826 local circulant source-operator endpoint-11 no-go audit\n\n"
        "`P2876/S1826` enumerates all `3^5=243` ternary radius-2 local circulant stencils on `Z12`, including the `3^3=27` reflection-symmetric subfamily, and applies them to intrinsic Fourier inputs `k=0..6`.  Because circulant local operators diagonalize on Fourier modes, every checked output is either zero or full-support; no stencil computes singleton endpoint `11` without importing a localized delta/phase seed.  No strict local endpoint-11 source operator, unit-bearing `9/5` coupling theorem, boundary/orientation source law, selector/localizer source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2876/S1826 local circulant source operator `L_total` guard",
        "## P2876/S1826 local circulant source operator `L_total` guard\n\n"
        "`P2876/S1826` adds no strict action term.  Local circulant stencils and their reflection-symmetric subfamily do not export a localized unit-bearing boundary/source density at `11`, non-premise support/origin law, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current local circulant source-operator endpoint-11 no-go guardrail (P2876/S1826, 2026-06-18)",
        "## Current local circulant source-operator endpoint-11 no-go guardrail (P2876/S1826, 2026-06-18)\n\n"
        "- P2876 enumerates all ternary radius-2 local circulant stencils on `Z12`, plus the reflection-symmetric subfamily, after P2875.\n"
        "- On intrinsic Fourier inputs `k=0..6`, every output is zero or full-support; singleton endpoint `11` appears only if a localized delta/phase seed is imported.\n"
        "- Do not promote local circulant stencils, reflection-symmetric local stencils, Fourier-diagonal operator action, two-dimensional `D12` irrep waves, full-DFT delta reconstruction, imported phase-11 coefficients, or earlier endpoint-label imports to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export a non-circulant strict local defect/source with independently exported origin/support at endpoint `11` and a unit-bearing `9/5` coefficient theorem, or use a genuinely different typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
