#!/usr/bin/env python3
"""Scratch probe: stage-by-stage completion ladder from legacy carrier to strict kernel.

The previous factorization certificate proved the one-line product identity.  This
probe makes the same bridge candidate more transparent and checkable by applying
completion in explicit stages on the finite Z12 distance domain:

    K0(d) = alpha_geo*cos(omega_L*d+phi_L)/(1+beta_tors*d)
    K1(d) = alpha_geo^{-1} K0(d)
    K2(d) = [cos(strict phase)/cos(legacy phase)] K1(d)
    K3(d) = [(1+beta_tors*d)/(1+beta*d^eta)] K2(d) = K_strict(d)

This is the requested "how legacy is completed into strict" computation.  It is
still not a derivation of the stages from nadsoliton dynamics, nor a proof of
beta_tors -> chi_11.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_kernel_completion_ladder_report.json"
OUT_MD = HERE / "bridge_strict_kernel_completion_ladder_report.md"
FACTORIZATION = HERE / "bridge_completed_strict_kernel_factorization_certificate_report.json"
BLOCKER_LATTICE = HERE / "bridge_completed_kernel_blocker_dependency_lattice_report.json"

ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
DOMAIN = list(range(12))
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def cos_legacy(d: float) -> float:
    return math.cos(LEGACY["omega"] * d + LEGACY["phi"])


def cos_strict(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"])


def legacy_full(d: float) -> float:
    return ALPHA_GEO * cos_legacy(d) / (1.0 + LEGACY["beta_tors"] * d)


def strict_kernel(d: float) -> float:
    return cos_strict(d) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def completion_ladder_row(d: int) -> dict[str, Any]:
    df = float(d)
    amplitude_factor = 1.0 / ALPHA_GEO
    phase_factor = cos_strict(df) / cos_legacy(df)
    damping_factor = (1.0 + LEGACY["beta_tors"] * df) / (1.0 + STRICT["beta"] * df ** STRICT["eta"])

    stage0_legacy_full = legacy_full(df)
    stage1_alpha_removed = amplitude_factor * stage0_legacy_full
    stage2_phase_transported = phase_factor * stage1_alpha_removed
    stage3_damping_compressed = damping_factor * stage2_phase_transported
    target = strict_kernel(df)

    return {
        "d": d,
        "stage0_legacy_full": stage0_legacy_full,
        "stage1_alpha_removed_legacy_carrier": stage1_alpha_removed,
        "stage2_phase_frequency_transported": stage2_phase_transported,
        "stage3_damping_compressed_strict_kernel": stage3_damping_compressed,
        "strict_kernel_target": target,
        "amplitude_factor_alpha_removal": amplitude_factor,
        "phase_frequency_transport_factor": phase_factor,
        "damping_compression_factor": damping_factor,
        "delta_alpha_step": stage1_alpha_removed - stage0_legacy_full,
        "delta_phase_step": stage2_phase_transported - stage1_alpha_removed,
        "delta_damping_step": stage3_damping_compressed - stage2_phase_transported,
        "strict_residual_after_all_steps": stage3_damping_compressed - target,
    }


def l1(values: list[float]) -> float:
    return sum(abs(value) for value in values)


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def sign_pattern(values: list[float]) -> list[int]:
    return [1 if value > 0 else -1 if value < 0 else 0 for value in values]


def monotone_decreasing(values: list[float]) -> bool:
    return all(left > right for left, right in zip(values, values[1:]))


def build_payload() -> dict[str, Any]:
    factorization = load_json(FACTORIZATION)
    blocker = load_json(BLOCKER_LATTICE)
    rows = [completion_ladder_row(d) for d in DOMAIN]
    alpha_deltas = [row["delta_alpha_step"] for row in rows]
    phase_deltas = [row["delta_phase_step"] for row in rows]
    damping_deltas = [row["delta_damping_step"] for row in rows]
    residuals = [row["strict_residual_after_all_steps"] for row in rows]
    damping_factors = [row["damping_compression_factor"] for row in rows]
    phase_factors = [row["phase_frequency_transport_factor"] for row in rows]

    stage_norms = {
        "stage0_legacy_full_l2": l2([row["stage0_legacy_full"] for row in rows]),
        "stage1_alpha_removed_l2": l2([row["stage1_alpha_removed_legacy_carrier"] for row in rows]),
        "stage2_phase_transported_l2": l2([row["stage2_phase_frequency_transported"] for row in rows]),
        "stage3_strict_l2": l2([row["stage3_damping_compressed_strict_kernel"] for row in rows]),
    }

    return {
        "result_kind": "SCRATCH_STRICT_KERNEL_COMPLETION_LADDER__STAGEWISE_CERTIFICATE_NOT_DERIVATION",
        "status": "legacy-carrier-completed-to-current-strict-kernel-by-explicit-stage-ladder",
        "source_reports": {
            "factorization_certificate": str(FACTORIZATION.relative_to(HERE.parents[1])),
            "blocker_lattice": str(BLOCKER_LATTICE.relative_to(HERE.parents[1])),
        },
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "legacy": LEGACY,
            "strict": STRICT,
            "domain_d_values": DOMAIN,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "stage_definitions": [
            {
                "stage": "stage0_legacy_full",
                "formula": "alpha_geo*cos(omega_L*d+phi_L)/(1+beta_tors*d)",
                "meaning": "historical legacy/nadsoliton-characteristic carrier with explicit alpha_geo amplitude",
            },
            {
                "stage": "stage1_alpha_removed_legacy_carrier",
                "formula": "alpha_geo^{-1}*stage0",
                "meaning": "remove explicit legacy amplitude so strict gate amplitude convention can be used",
            },
            {
                "stage": "stage2_phase_frequency_transported",
                "formula": "[cos(omega_S*d+phi_S)/cos(omega_L*d+phi_L)]*stage1",
                "meaning": "transport legacy torsion/resonance phase to strict gate phase",
            },
            {
                "stage": "stage3_damping_compressed_strict_kernel",
                "formula": "[(1+beta_tors*d)/(1+beta*d^eta)]*stage2",
                "meaning": "replace legacy hyperbolic beta_tors*d damping by strict beta*d^eta compression",
            },
        ],
        "completion_summary": {
            "max_abs_final_residual": max(abs(value) for value in residuals),
            "final_residual_tolerance_pass": max(abs(value) for value in residuals) < 1e-15,
            "amplitude_factor": 1.0 / ALPHA_GEO,
            "damping_factor_positive": all(value > 0 for value in damping_factors),
            "damping_factor_nonincreasing": monotone_decreasing(damping_factors),
            "damping_factor_d0": damping_factors[0],
            "damping_factor_d11": damping_factors[-1],
            "phase_factor_sign_pattern": sign_pattern(phase_factors),
            "stage_l1_delta_alpha": l1(alpha_deltas),
            "stage_l1_delta_phase": l1(phase_deltas),
            "stage_l1_delta_damping": l1(damping_deltas),
            "stage_norms": stage_norms,
        },
        "completion_ladder_rows": rows,
        "blocker_context": {
            "factorization_status": factorization["status"],
            "blocker_lattice_status": blocker["status"],
            "still_open_for_theorem_level_bridge": blocker["current_frontier_summary"]["remaining_for_theorem_level_bridge"],
            "main_one_bit_blocker": blocker["current_frontier_summary"]["main_one_bit_blocker"],
            "transport_blocker": blocker["current_frontier_summary"]["transport_blocker"],
        },
        "exact_proof_certificate": {
            "stage_identity": "For every d in 0..11, stage3_damping_compressed_strict_kernel equals K_strict_gate(d) to floating precision.",
            "how_legacy_is_completed": "Legacy carrier is completed by alpha removal, strict phase/frequency transport, and strict eta damping compression; orientation/chi_11 remains a separate source theorem-target.",
            "nonduplication": "This is a stage-ladder explanation built on the prior factorization certificate, not another orbit/Reynolds/Puiseux enumeration.",
            "theoretical_limit": "The ladder shows what factors complete the carrier, not why nadsoliton dynamics must generate those factors.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state; legacy kernel data is read as a historical internal characteristic carrier completed into the current strict kernel.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate is the current live/full kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is asserted; the completion ladder uses explicit factors.",
            "No proof derives the completion factors from strict nadsoliton dynamics yet.",
            "No beta_tors -> chi_11 theorem is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["completion_summary"]
    blockers = payload["blocker_context"]
    lines = [
        "# Strict kernel completion ladder",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Stages",
        "",
    ]
    for stage in payload["stage_definitions"]:
        lines.extend(
            [
                f"### {stage['stage']}",
                f"- Formula: `{stage['formula']}`",
                f"- Meaning: {stage['meaning']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Completion summary",
            "",
            f"- Max final residual: `{summary['max_abs_final_residual']:.3e}`",
            f"- Residual tolerance pass: `{summary['final_residual_tolerance_pass']}`",
            f"- Amplitude factor: `{summary['amplitude_factor']:.12f}`",
            f"- Damping factor positive/nonincreasing: `{summary['damping_factor_positive']}` / `{summary['damping_factor_nonincreasing']}`",
            f"- Damping factor d=0 to d=11: `{summary['damping_factor_d0']:.12f}` -> `{summary['damping_factor_d11']:.12f}`",
            f"- Phase factor sign pattern: `{summary['phase_factor_sign_pattern']}`",
            f"- L1 deltas alpha/phase/damping: `{summary['stage_l1_delta_alpha']:.6f}` / `{summary['stage_l1_delta_phase']:.6f}` / `{summary['stage_l1_delta_damping']:.6f}`",
            "",
            "## Remaining blockers",
            "",
            f"- Still open for theorem-level bridge: `{blockers['still_open_for_theorem_level_bridge']}`",
            f"- Main one-bit blocker: {blockers['main_one_bit_blocker']}",
            f"- Transport blocker: {blockers['transport_blocker']}",
            "",
            "## Proof certificate",
            "",
        ]
    )
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
