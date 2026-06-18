#!/usr/bin/env python3
"""P2858/S1808: phase-bit cell continuum no-source audit.

P2857 showed that observer readout cannot turn P2856 phase-bit ambiguity into a
pre-observer source.  This packet strengthens that boundary from sampled finite
witnesses to a direct open-cell certificate: the strict phase-bit profile is
locally constant on a positive-radius neighborhood of (omega, phi).

Result: the audited Z12 phase-bit/topological readout is a robust finite
witness, but it cannot select exact omega=743/4000 and phi=13/80.  A continuum
of nearby phase/frequency pairs has the same phase-bit profile, so any source
law must add information beyond phase bits, observer frames, and sign-cell
bookkeeping.
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
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import (
    GEN,
    LEGACY_OMEGA,
    LEGACY_PHI,
    STRICT_OMEGA,
    STRICT_PHI,
    fraction_payload,
    phase_bits_for,
)

P2857 = GEN / "p2857_s1807_observer_readout_phase_source_effect_audit.json"
OUT = GEN / "p2858_s1808_phase_bit_cell_continuum_no_source_audit.json"
MD = GEN / "p2858_s1808_phase_bit_cell_continuum_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

DOMAIN = tuple(range(12))
ZERO_HALF_PERIOD = math.pi / 2.0
PI = math.pi


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    raise ValueError("zero sign has no stable phase-bit cell")


def legacy_required_sign(bit: int, legacy_cos_sign: int) -> int:
    return legacy_cos_sign if bit == 0 else -legacy_cos_sign


def distance_to_cos_zero(theta: float) -> float:
    nearest_index = round((theta - ZERO_HALF_PERIOD) / PI)
    nearest_zero = ZERO_HALF_PERIOD + nearest_index * PI
    return abs(theta - nearest_zero)


def phase_cell_rows(bits: list[int]) -> list[dict[str, Any]]:
    rows = []
    for d, bit in enumerate(bits):
        legacy_cos = math.cos(LEGACY_OMEGA * d + LEGACY_PHI)
        strict_theta = float(STRICT_OMEGA) * d + float(STRICT_PHI)
        strict_cos = math.cos(strict_theta)
        required = legacy_required_sign(bit, sign(legacy_cos))
        observed = sign(strict_cos)
        margin = distance_to_cos_zero(strict_theta)
        rows.append(
            {
                "d": d,
                "bit": bit,
                "legacy_cos_sign": sign(legacy_cos),
                "required_candidate_cos_sign": required,
                "strict_candidate_cos_sign": observed,
                "sign_constraint_satisfied": required == observed,
                "theta": strict_theta,
                "distance_to_nearest_cos_zero": margin,
                "safe_common_epsilon_bound": margin / (d + 1),
            }
        )
    return rows


def common_open_box(rows: list[dict[str, Any]]) -> dict[str, Any]:
    limiting = min(rows, key=lambda row: row["safe_common_epsilon_bound"])
    epsilon = limiting["safe_common_epsilon_bound"] / 2.0
    rational_probe_delta = Fraction(1, 1_000_000)
    probes = [
        (STRICT_OMEGA + rational_probe_delta, STRICT_PHI),
        (STRICT_OMEGA - rational_probe_delta, STRICT_PHI),
        (STRICT_OMEGA, STRICT_PHI + rational_probe_delta),
        (STRICT_OMEGA, STRICT_PHI - rational_probe_delta),
        (STRICT_OMEGA + rational_probe_delta, STRICT_PHI - rational_probe_delta),
    ]
    return {
        "limiting_domain_point": limiting["d"],
        "min_safe_common_epsilon_bound": limiting["safe_common_epsilon_bound"],
        "certified_open_box_half_width": epsilon,
        "rational_probe_delta": f"{rational_probe_delta.numerator}/{rational_probe_delta.denominator}",
        "rational_probe_delta_inside_box": float(rational_probe_delta) < epsilon,
        "probe_phase_bits": [phase_bits_for(omega, phi) for omega, phi in probes],
        "probe_parameters": [
            {"omega": fraction_payload(omega), "phi": fraction_payload(phi)} for omega, phi in probes
        ],
    }


def build_payload(p2857: dict[str, Any]) -> dict[str, Any]:
    bits = p2857["observer_readout_phase_source_effect_audit"]["target_phase_bits"]
    rows = phase_cell_rows(bits)
    box = common_open_box(rows)
    strict_bits = phase_bits_for(STRICT_OMEGA, STRICT_PHI)
    probe_bits_match = all(probe == bits for probe in box["probe_phase_bits"])
    all_sign_constraints = all(row["sign_constraint_satisfied"] for row in rows)
    facts = {
        "p2857_rechecked": p2857.get("status") == "P2857_OBSERVER_READOUT_PHASE_SOURCE_EFFECT_AUDIT_NO_CLOSURE",
        "strict_bits_match_input": strict_bits == bits,
        "all_phase_cell_sign_constraints_satisfied": all_sign_constraints,
        "positive_open_box_certified": box["certified_open_box_half_width"] > 0.0,
        "rational_probes_inside_box": box["rational_probe_delta_inside_box"],
        "rational_probes_preserve_phase_bits": probe_bits_match,
    }
    return {
        "status": "P2858_PHASE_BIT_CELL_CONTINUUM_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2857": sha(P2857)},
        "phase_bit_cell_continuum_no_source_audit": {
            "input_status_rechecked": p2857.get("status"),
            "strict_omega": fraction_payload(STRICT_OMEGA),
            "strict_phi": fraction_payload(STRICT_PHI),
            "target_phase_bits": bits,
            "cell_rows": rows,
            "open_box_certificate": box,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "sign_cell_step": "For each d in Z12, the strict theta=omega*d+phi has positive distance to the nearest zero of cos(theta).",
                "open_box_step": "If |delta_omega| and |delta_phi| are smaller than the certified common half-width, every audited sign inequality remains unchanged.",
                "continuum_step": "Therefore a continuum of nearby omega/phi pairs shares the same phase-bit profile; phase-bit/topological readout is not an exact source-selection law.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_phase_bit_open_cell_no_source_audit": all(facts.values()),
            "exports_strict_phase_frequency_source_law": False,
        },
        "decision": {
            "negative_export_flags": {
                "phase_bit_source_law_exported": False,
                "observer_source_law_exported": False,
                "prime5_phase_unit_source_exported": False,
                "strict_phase_frequency_source_law_exported": False,
                "selector_closure_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2858 proves a positive-radius phase-bit cell around the strict omega/phi tuple.  The audited Z12 phase-bit profile is stable under a continuum of perturbations, including explicit rational probes, so phase bits, observer readout, and sign-cell bookkeeping cannot select the exact strict tuple without a new pre-observer source law.",
            "next_honest_step": "Do not replay phase-bit, observer, affine-transport, rational-lattice, or sign-cell evidence as a source.  The next proof-grade move must add a genuinely new pre-observer source-selection law for the exact phase/frequency tuple, or pivot to a genuinely new eta/beta source law; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["phase_bit_cell_continuum_no_source_audit"]
    box = audit["open_box_certificate"]
    lines = [
        "# P2858/S1808 phase-bit cell continuum no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Open-cell certificate",
        f"- limiting domain point: `{box['limiting_domain_point']}`",
        f"- min safe common epsilon bound: `{box['min_safe_common_epsilon_bound']}`",
        f"- certified open-box half-width: `{box['certified_open_box_half_width']}`",
        f"- rational probe delta: `{box['rational_probe_delta']}`",
        f"- rational probe delta inside box: `{box['rational_probe_delta_inside_box']}`",
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
    payload = build_payload(read_json(P2857))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2858/S1808 phase-bit cell continuum no-source audit",
        "## P2858/S1808 phase-bit cell continuum no-source audit\n\n"
        "`P2858/S1808` strengthens the P2856/P2857 ambiguity result from sampled witnesses to an open-cell certificate.  Each audited strict `theta_d=omega*d+phi` has positive distance from the nearest cosine zero, yielding a positive-radius neighborhood of `omega=743/4000`, `phi=13/80` with the same `Z12` phase-bit profile; explicit rational probes inside the box preserve the bits.  Thus phase bits, observer readout, affine/sign-cell bookkeeping, and rational-lattice representability cannot select the exact tuple without a new pre-observer source law.  No strict phase/frequency source law, full kernel bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2858/S1808 phase-bit open cell `L_total` guard",
        "## P2858/S1808 phase-bit open cell `L_total` guard\n\n"
        "`P2858/S1808` adds no action term.  A positive-radius phase-bit cell and rational perturbation witnesses do not provide a pre-observer unit-bearing source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current phase-bit open-cell continuum guardrail (P2858/S1808, 2026-06-18)",
        "## Current phase-bit open-cell continuum guardrail (P2858/S1808, 2026-06-18)\n\n"
        "- P2858 proves a positive-radius open phase-bit cell around strict `omega=743/4000`, `phi=13/80`; explicit rational perturbations preserve the same `Z12` phase-bit profile.\n"
        "- Phase-bit/topological readout is therefore a robust finite witness but not an exact source-selection law for the strict tuple.\n"
        "- Do not promote phase bits, observer readout, affine transport, sign-cell stability, rational probes, or prime-5 representability to strict phase/frequency source law, selector closure, full bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a genuinely new pre-observer source-selection law for exact `omega/phi`, provide a genuinely new `eta/beta` source law, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
