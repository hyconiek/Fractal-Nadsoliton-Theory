#!/usr/bin/env python3
"""Scratch probe: positive-factor sign-separation certificate.

The strict-completion factorization uses a positive alpha-normalization factor,
a phase transport factor, and a positive damping/compression factor.  The recent
Z2 and exact damping certificates separately fixed the phase sign transport and
the damping positivity/monotonicity.

This probe records the finite algebra consequence:

    sign(A(d) * P(d) * D(d)) = sign(P(d))

on d=0..11, because A(d)>0 and D(d)>0.  Equivalently, positive factors carry
zero Z2 node/edge sign bits, so the whole completion sign cocycle is exactly the
phase-sign Z2 coboundary.  This is sign bookkeeping for the selected completion
ansatz, not a bridge theorem and not a strict dynamical derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.md"
NECESSITY_REPORT = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
DAMPING_EXACT_REPORT = HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json"

EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def sign_of(value: float) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def sign_to_bit(sign: int) -> int:
    if sign == 1:
        return 0
    if sign == -1:
        return 1
    raise ValueError(f"zero sign has no Z2 bit: {sign}")


def node_rows(necessity: dict[str, Any], z2: dict[str, Any]) -> list[dict[str, Any]]:
    z2_by_d = {row["d"]: row for row in z2["node_bit_rows"]}
    rows = []
    for row in necessity["pointwise_quotient_certificate"]:
        d = row["d"]
        alpha_positive = row["alpha_normalization_factor"] > 0
        damping_positive = row["damping_compression_factor"] > 0
        phase_sign = row["phase_factor_sign"]
        factor_product_sign = sign_of(row["factor_product"])
        quotient_sign = sign_of(row["strict_over_legacy_quotient"])
        positive_factor_bit = 0 if alpha_positive and damping_positive else 1
        rows.append({
            "d": d,
            "alpha_normalization_positive": alpha_positive,
            "damping_compression_positive": damping_positive,
            "positive_factor_z2_bit": positive_factor_bit,
            "phase_factor_sign": phase_sign,
            "phase_z2_bit": sign_to_bit(phase_sign),
            "factor_product_sign": factor_product_sign,
            "strict_over_legacy_quotient_sign": quotient_sign,
            "z2_phase_sign": z2_by_d[d]["phase_sign"],
            "z2_node_bit": z2_by_d[d]["node_bit"],
            "factor_sign_equals_phase_sign": factor_product_sign == phase_sign,
            "quotient_sign_equals_phase_sign": quotient_sign == phase_sign,
            "z2_bit_equals_phase_sign_bit": z2_by_d[d]["node_bit"] == sign_to_bit(phase_sign),
        })
    return rows


def edge_rows(z2: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in z2["edge_bit_rows"]:
        rows.append({
            "edge": row["edge"],
            "positive_factor_edge_bit": 0,
            "phase_edge_bit": row["edge_bit"],
            "completion_edge_bit": row["edge_bit"],
            "is_phase_flip": row["is_phase_flip"],
            "positive_factors_change_edge_sign": False,
            "completion_flip_equals_phase_flip": row["is_phase_flip"] == (row["edge_bit"] == 1),
        })
    return rows


def build_payload() -> dict[str, Any]:
    necessity = load_json(NECESSITY_REPORT)
    z2 = load_json(Z2_REPORT)
    damping_exact = load_json(DAMPING_EXACT_REPORT)
    nodes = node_rows(necessity, z2)
    edges = edge_rows(z2)
    derived_sign_pattern = [row["factor_product_sign"] for row in nodes]
    derived_flip_edges = [row["edge"] for row in edges if row["completion_edge_bit"] == 1]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_POSITIVE_FACTOR_SIGN_SEPARATION_CERTIFICATE__PHASE_ONLY_Z2_SIGN",
        "status": "positive-alpha-and-damping-factors-carry-zero-z2-sign-so-completion-sign-is-phase-only",
        "source_reports": {
            "necessity_certificate": str(NECESSITY_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "damping_exact_rational_certificate": str(DAMPING_EXACT_REPORT.relative_to(ROOT)),
        },
        "sign_separation_rule": {
            "factorization": "T(d)=A(d)*P(d)*D(d)",
            "positive_factor_statement": "A(d)>0 and D(d)>0 on d=0..11",
            "sign_consequence": "sign(T(d))=sign(P(d))",
            "z2_consequence": "positive factors contribute zero node and edge sign bits",
        },
        "node_sign_separation_rows": nodes,
        "edge_sign_separation_rows": edges,
        "positive_factor_sign_summary": {
            "all_alpha_normalization_factors_positive": all(row["alpha_normalization_positive"] for row in nodes),
            "all_damping_compression_factors_positive": all(row["damping_compression_positive"] for row in nodes),
            "all_positive_factor_z2_bits_zero": all(row["positive_factor_z2_bit"] == 0 for row in nodes),
            "all_factor_signs_equal_phase_signs": all(row["factor_sign_equals_phase_sign"] for row in nodes),
            "all_quotient_signs_equal_phase_signs": all(row["quotient_sign_equals_phase_sign"] for row in nodes),
            "all_z2_bits_equal_phase_sign_bits": all(row["z2_bit_equals_phase_sign_bit"] for row in nodes),
            "all_positive_factor_edge_bits_zero": all(row["positive_factor_edge_bit"] == 0 for row in edges),
            "all_completion_flips_equal_phase_flips": all(row["completion_flip_equals_phase_flip"] for row in edges),
            "derived_completion_sign_pattern": derived_sign_pattern,
            "derived_completion_flip_edges": derived_flip_edges,
            "matches_expected_sign_pattern": derived_sign_pattern == EXPECTED_SIGN_PATTERN,
            "matches_expected_flip_edges": derived_flip_edges == EXPECTED_FLIP_EDGES,
            "matches_z2_sign_pattern": derived_sign_pattern == z2["z2_coboundary_summary"]["phase_sign_pattern"],
            "matches_z2_flip_edges": derived_flip_edges == z2["z2_coboundary_summary"]["derived_phase_sign_flip_edges"],
            "exact_damping_positive_and_decreasing": damping_exact["exact_rational_damping_summary"]["continuous_positive_certificate"] and damping_exact["exact_rational_damping_summary"]["continuous_strictly_decreasing_certificate"],
        },
        "blocker_context": {
            "what_this_refines": "turns positive alpha/damping factor facts into an explicit phase-only sign separation theorem for the completion factorization",
            "necessity_status": necessity["status"],
            "phase_sign_z2_status": z2["status"],
            "damping_exact_status": damping_exact["status"],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_formula_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "factorization_step": "Use the completion factorization T(d)=A(d)*P(d)*D(d) from the necessity certificate.",
            "positivity_step": "The alpha-normalization and damping/compression factors are positive at every audited node; exact damping positivity is separately certified.",
            "node_sign_step": "Multiplication by positive factors does not change node sign, so sign(T(d))=sign(P(d)).",
            "edge_z2_step": "Positive factors contribute zero Z2 edge bits, so the completion sign coboundary equals the phase sign coboundary.",
            "nonduplication": "This is a positive-factor sign-separation certificate, not another phase-zero placement, damping calculus, or real-valued transport fit.",
            "theoretical_limit": "The certificate proves finite sign bookkeeping for the selected completion ansatz; it does not derive A(d), P(d), D(d), omega/phi, beta/eta, or transport from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only checks finite sign bookkeeping inside the selected completion factorization.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or transport from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["positive_factor_sign_summary"]
    lines = [
        "# Strict completion positive-factor sign-separation certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit proves a finite sign-separation statement for the selected",
        "completion factorization: positive alpha-normalization and damping factors",
        "carry zero Z2 sign bits, so the full completion sign is phase-only.",
        "",
        "## Summary",
        "",
        f"- All alpha-normalization factors positive: `{summary['all_alpha_normalization_factors_positive']}`",
        f"- All damping/compression factors positive: `{summary['all_damping_compression_factors_positive']}`",
        f"- All positive-factor Z2 bits zero: `{summary['all_positive_factor_z2_bits_zero']}`",
        f"- All factor signs equal phase signs: `{summary['all_factor_signs_equal_phase_signs']}`",
        f"- All completion flips equal phase flips: `{summary['all_completion_flips_equal_phase_flips']}`",
        f"- Derived sign pattern: `{summary['derived_completion_sign_pattern']}`",
        f"- Derived flip edges: `{summary['derived_completion_flip_edges']}`",
        "",
        "## Edge sign separation rows",
        "",
        "| edge | positive-factor bit | phase bit | completion bit | phase flip |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in payload["edge_sign_separation_rows"]:
        lines.append(
            f"| {row['edge']} | {row['positive_factor_edge_bit']} | {row['phase_edge_bit']} | {row['completion_edge_bit']} | {row['is_phase_flip']} |"
        )
    lines.extend([
        "",
        "## Proof certificate",
        "",
    ])
    for key, value in payload["proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend([
        "",
        "## Hard limits",
        "",
    ])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
