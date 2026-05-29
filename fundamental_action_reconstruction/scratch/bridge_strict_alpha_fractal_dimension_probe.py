#!/usr/bin/env python3
"""Scratch probe: strict-alpha fractal-dimension candidate for eta.

This is the next narrow bridge-prep step after the compression ontology audit. It
asks whether the strict-side Shannon source already exported in the repo,

    alpha_geo_strict_derived_v1 = ln(16) = 4 ln 2,

can supply a compression scale for the strict exponent eta=9/5 without importing
legacy physical roles.  The smallest honest calculation is a Hausdorff-style
count/scale dimension using the strict nad12 support size N=12 and the Shannon
microstate square-root scale sqrt(exp(alpha_geo))=4:

    eta_alpha_z12 = ln(12) / ln(4) = 2 ln(12) / alpha_geo.

The result is close to eta=1.8 but not exact.  This packet therefore exports a
candidate discriminator, not a theorem: strict alpha plus nad12 count plausibly
explains why eta is near 1.8, while an independent correction/rationalization is
still required to prove the exact 9/5 gate exponent.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ALPHA_PACKET = ROOT / "fundamental_action_reconstruction" / "generated" / "alpha_geo_strict_derived_v1.json"
COMPRESSION_ONTOLOGY = HERE / "bridge_compression_ontology_audit_report.json"
MEASURE_TRANSPORT = HERE / "bridge_phase_normalized_measure_transport_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_fractal_dimension_report.json"
OUT_MD = HERE / "bridge_strict_alpha_fractal_dimension_report.md"

NAD12_SUPPORT_SIZE = 12
STRICT_TARGET_ETA = 9.0 / 5.0
INTEGER_SCALE_CANDIDATES = list(range(2, 13))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def dimension_from_scale(n: int, scale: float) -> float:
    return math.log(n) / math.log(scale)


def candidate_rows() -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for scale in INTEGER_SCALE_CANDIDATES:
        eta = dimension_from_scale(NAD12_SUPPORT_SIZE, float(scale))
        rows.append(
            {
                "label": f"ln(12)/ln({scale})",
                "scale": scale,
                "eta_candidate": eta,
                "signed_residual_vs_9_5": eta - STRICT_TARGET_ETA,
                "abs_residual_vs_9_5": abs(eta - STRICT_TARGET_ETA),
            }
        )
    rows.sort(key=lambda row: float(row["abs_residual_vs_9_5"]))
    return rows


def main() -> None:
    alpha_packet = load_json(ALPHA_PACKET)
    compression_report = load_json(COMPRESSION_ONTOLOGY)
    measure_report = load_json(MEASURE_TRANSPORT)

    alpha_geo = 4.0 * math.log(2.0)
    alpha_scale = math.sqrt(math.exp(alpha_geo))
    eta_alpha_z12 = dimension_from_scale(NAD12_SUPPORT_SIZE, alpha_scale)
    eta_residual = eta_alpha_z12 - STRICT_TARGET_ETA
    exact_scale_for_eta = NAD12_SUPPORT_SIZE ** (1.0 / STRICT_TARGET_ETA)
    effective_count_at_scale4_eta = alpha_scale**STRICT_TARGET_ETA
    integer_rows = candidate_rows()
    best_integer = integer_rows[0]

    report = {
        "status": "OPEN_STRICT_ALPHA_FRACTAL_DIMENSION_CANDIDATE_NO_ETA_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_FRACTAL_DIMENSION_PROBE__NOT_A_THEOREM",
        "source_packets": {
            "alpha_geo_strict_derived_v1": str(ALPHA_PACKET.relative_to(ROOT)),
            "compression_ontology": str(COMPRESSION_ONTOLOGY.relative_to(ROOT)),
            "measure_transport": str(MEASURE_TRANSPORT.relative_to(ROOT)),
        },
        "strict_sources_used": {
            "alpha_geo_strict_derived_v1": alpha_packet["value"],
            "alpha_geo_numeric": alpha_geo,
            "nad12_support_size": NAD12_SUPPORT_SIZE,
            "scale_from_alpha": "sqrt(exp(alpha_geo_strict_derived_v1)) = 4",
            "scale_from_alpha_numeric": alpha_scale,
            "guardrail": "The calculation treats 4 ln 2 as strict-derived Shannon source data, not as a legacy kernel role transfer.",
        },
        "fractal_dimension_candidate": {
            "formula": "eta_alpha_z12 = ln(12)/ln(sqrt(exp(alpha_geo_strict_derived_v1))) = 2 ln(12)/alpha_geo",
            "eta_alpha_z12": eta_alpha_z12,
            "strict_target_eta_9_5": STRICT_TARGET_ETA,
            "signed_residual_vs_9_5": eta_residual,
            "abs_residual_vs_9_5": abs(eta_residual),
            "relative_abs_residual_vs_9_5": abs(eta_residual) / STRICT_TARGET_ETA,
            "within_one_percent_of_9_5": abs(eta_residual) / STRICT_TARGET_ETA < 0.01,
        },
        "inverse_checks": {
            "exact_scale_required_for_eta_9_5_with_count_12": exact_scale_for_eta,
            "scale_residual_exact_minus_alpha_scale": exact_scale_for_eta - alpha_scale,
            "relative_scale_residual": abs(exact_scale_for_eta - alpha_scale) / alpha_scale,
            "effective_count_if_scale4_and_eta9_5": effective_count_at_scale4_eta,
            "effective_count_residual_vs_12": effective_count_at_scale4_eta - NAD12_SUPPORT_SIZE,
            "relative_effective_count_residual": abs(effective_count_at_scale4_eta - NAD12_SUPPORT_SIZE) / NAD12_SUPPORT_SIZE,
        },
        "integer_scale_discriminator": {
            "best_integer_scale_by_abs_residual": best_integer,
            "all_integer_scale_rows_sorted": integer_rows,
            "interpretation": "Among integer scales 2..12, the strict-alpha scale 4 is the closest to eta=9/5 for a 12-state support.",
        },
        "upstream_replay": {
            "compression_audit_says_strict_terminal_candidate": compression_report["strict_final_form_assessment"]["is_strict_plausible_terminal_shape_candidate"],
            "compression_audit_says_not_final_formula": compression_report["strict_final_form_assessment"]["is_strict_proven_final_nadsoliton_formula"],
            "measure_transport_candidate_supported": measure_report["candidate_ontological_reading"]["supported_by_this_probe"],
            "measure_transport_balance_residual": measure_report["measure_balance"]["max_abs_balance_residual"],
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                best_integer["scale"] == 4
                and abs(eta_residual) / STRICT_TARGET_ETA < 0.01
                and abs(exact_scale_for_eta - alpha_scale) / alpha_scale < 0.01
            ),
            "content": "Strict alpha supplies a natural scale 4; nad12 supplies count 12; their Hausdorff-style dimension is 1.79248, close to eta=1.8.",
            "why_this_is_more_proof_like": "It uses an already exported strict Shannon value and a fixed nad12 count, then checks inverse scale/count residuals and competing integer scales.",
            "why_this_is_not_enough": "The calculation is close but not exact; it does not derive the rational correction from 1.79248 to 9/5 or prove that this dimension is the nadsoliton compression law.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives eta=9/5 exactly from alpha_geo_strict_derived_v1 and nad12 support.",
            "No theorem proves that ln(12)/ln(4) is the intrinsic nadsoliton Hausdorff dimension.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Search for a strict-side finite-size/rationalization correction that moves ln(12)/ln(4)=1.79248 to eta=9/5 without importing legacy roles; if none is found, record eta=9/5 as gate-selected with a close strict-alpha fractal-dimension shadow.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha fractal-dimension probe\n\n"
        "Status: strict alpha + nad12 count gives a close eta shadow; no eta theorem.\n\n"
        f"- Candidate: `ln(12)/ln(4) = {eta_alpha_z12:.12f}` vs strict `eta=9/5={STRICT_TARGET_ETA:.12f}`; residual `{eta_residual:.3e}`.\n"
        f"- Exact scale for `eta=9/5` and count 12: `{exact_scale_for_eta:.12f}` vs alpha scale `4`; relative scale residual `{abs(exact_scale_for_eta-alpha_scale)/alpha_scale:.3e}`.\n"
        f"- Integer-scale discriminator: best scale `{best_integer['scale']}` with abs residual `{best_integer['abs_residual_vs_9_5']:.3e}`.\n"
        "- Honest read: strict likely has a real compression-dimension shadow, but exact eta=9/5 still needs a strict-side correction or derivation.\n"
        "- No false pass: no kernel identity, no legacy role transfer, no exact eta theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
