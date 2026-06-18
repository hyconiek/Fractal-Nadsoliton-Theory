#!/usr/bin/env python3
"""P2851/S1801: amplitude-normalization kernel bridge-atom audit.

P2850 leaves the next proof-grade move as one typed bridge-source atom.  P2851
chooses the alternative atom not covered by P2849: the amplitude/normalization
passage from the legacy prefactor alpha_geo = 4 ln 2 to the strict gate kernel,
which has no explicit amplitude prefactor.

The audit is intentionally narrow.  It asks whether a single target-independent
constant amplitude A can map the strict gate shape to the legacy kernel values
on audited distances, or whether the legacy alpha_geo prefactor can simply be
removed/absorbed without the remaining phase/frequency/damping atoms.  It does
not claim a full bridge theorem or role transfer.
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
P2850 = GEN / "p2850_s1800_eml_single_operator_kernel_bridge_impact_audit.json"
OUT = GEN / "p2851_s1801_amplitude_normalization_kernel_bridge_atom_audit.json"
MD = GEN / "p2851_s1801_amplitude_normalization_kernel_bridge_atom_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY_OMEGA = math.pi / 4.0
LEGACY_PHI = math.pi / 6.0
LEGACY_BETA_TORS = 0.01
STRICT_OMEGA = 0.18575
STRICT_PHI = 0.16250
STRICT_BETA = 1.0
STRICT_ETA = 1.8
DISTANCES = tuple(range(0, 13))
REQUIRED_AMPLITUDE_BRIDGE_PREMISES = (
    "legacy_alpha_geo_defined",
    "strict_unit_amplitude_defined",
    "target_independent_constant_amplitude_map",
    "zero_residual_after_constant_amplitude_fit",
    "phase_frequency_compatibility",
    "damping_compression_compatibility",
    "alpha_geo_source_law_safe_for_strict_kernel",
    "role_transfer_theorem_available",
)


def legacy_shape(d: int) -> float:
    return math.cos(LEGACY_OMEGA * d + LEGACY_PHI) / (1.0 + LEGACY_BETA_TORS * d)


def legacy_kernel(d: int) -> float:
    return ALPHA_GEO * legacy_shape(d)


def strict_kernel(d: int) -> float:
    return math.cos(STRICT_OMEGA * d + STRICT_PHI) / (1.0 + STRICT_BETA * (d ** STRICT_ETA))


def pointwise_amplitude_ratio(d: int) -> float | None:
    s = strict_kernel(d)
    if abs(s) < 1e-15:
        return None
    return legacy_kernel(d) / s


def least_squares_amplitude() -> dict[str, float]:
    xs = [strict_kernel(d) for d in DISTANCES]
    ys = [legacy_kernel(d) for d in DISTANCES]
    numerator = sum(x * y for x, y in zip(xs, ys))
    denominator = sum(x * x for x in xs)
    a = numerator / denominator
    residuals = [a * x - y for x, y in zip(xs, ys)]
    sse = sum(r * r for r in residuals)
    max_abs = max(abs(r) for r in residuals)
    return {"best_fit_amplitude": a, "sse": sse, "max_abs_residual": max_abs}


def amplitude_rows() -> list[dict[str, object]]:
    rows = []
    for d in DISTANCES:
        strict = strict_kernel(d)
        legacy = legacy_kernel(d)
        ratio = pointwise_amplitude_ratio(d)
        rows.append(
            {
                "d": d,
                "legacy_kernel": legacy,
                "strict_kernel": strict,
                "pointwise_legacy_over_strict_amplitude": ratio,
                "legacy_without_alpha_shape": legacy_shape(d),
                "strict_minus_legacy_without_alpha_shape": strict - legacy_shape(d),
            }
        )
    return rows


def ratio_stats(rows: list[dict[str, object]]) -> dict[str, object]:
    ratios = [row["pointwise_legacy_over_strict_amplitude"] for row in rows if row["pointwise_legacy_over_strict_amplitude"] is not None]
    assert ratios
    return {
        "min_ratio": min(ratios),
        "max_ratio": max(ratios),
        "max_minus_min": max(ratios) - min(ratios),
        "sign_changes_present": any(r < 0 for r in ratios) and any(r > 0 for r in ratios),
        "constant_ratio": math.isclose(max(ratios), min(ratios), rel_tol=0.0, abs_tol=1e-12),
    }


def two_point_constant_amplitude_matrix(rows: list[dict[str, object]]) -> dict[str, object]:
    pairs = [(0, 1), (0, 2), (1, 2), (2, 4), (4, 8), (8, 12)]
    by_d = {row["d"]: row for row in rows}
    out = []
    for a, b in pairs:
        ra = by_d[a]["pointwise_legacy_over_strict_amplitude"]
        rb = by_d[b]["pointwise_legacy_over_strict_amplitude"]
        same = (ra is not None and rb is not None and math.isclose(float(ra), float(rb), rel_tol=0.0, abs_tol=1e-12))
        out.append({"pair": [a, b], "ratio_a": ra, "ratio_b": rb, "same_constant_amplitude": same})
    return {
        "proof_statement": "A target-independent amplitude-only bridge A*K_strict(d)=K_legacy(d) requires identical pointwise ratios K_legacy(d)/K_strict(d) at every nonzero strict point; two distinct ratios refute an amplitude-only bridge.",
        "rows": out,
        "all_pairs_reject_constant_amplitude": all(not row["same_constant_amplitude"] for row in out),
    }


def premise_matrix(rstats: dict[str, object], fit: dict[str, float], two_point: dict[str, object]) -> dict[str, object]:
    premises = {
        "legacy_alpha_geo_defined": True,
        "strict_unit_amplitude_defined": True,
        "target_independent_constant_amplitude_map": not rstats["sign_changes_present"] and not two_point["all_pairs_reject_constant_amplitude"],
        "zero_residual_after_constant_amplitude_fit": fit["max_abs_residual"] < 1e-12,
        "phase_frequency_compatibility": False,
        "damping_compression_compatibility": False,
        "alpha_geo_source_law_safe_for_strict_kernel": False,
        "role_transfer_theorem_available": False,
    }
    return {
        "premises": premises,
        "missing_premises": [key for key in REQUIRED_AMPLITUDE_BRIDGE_PREMISES if not premises[key]],
        "accepted_as_amplitude_normalization_bridge_atom": all(premises.values()),
    }


def build_payload(p2850: dict[str, object]) -> dict[str, object]:
    rows = amplitude_rows()
    rstats = ratio_stats(rows)
    fit = least_squares_amplitude()
    two_point = two_point_constant_amplitude_matrix(rows)
    matrix = premise_matrix(rstats, fit, two_point)
    facts = {
        "p2850_rechecked": p2850.get("status") == "P2850_EML_SINGLE_OPERATOR_KERNEL_BRIDGE_IMPACT_AUDIT_NO_CLOSURE",
        "legacy_alpha_and_strict_unit_amplitude_defined": matrix["premises"]["legacy_alpha_geo_defined"] and matrix["premises"]["strict_unit_amplitude_defined"],
        "two_point_ratios_reject_amplitude_only_bridge": two_point["all_pairs_reject_constant_amplitude"],
        "least_squares_residual_nonzero": fit["max_abs_residual"] > 1e-12,
        "ratio_not_constant": not rstats["constant_ratio"],
        "no_alpha_geo_source_law_safe_for_strict_kernel": not matrix["premises"]["alpha_geo_source_law_safe_for_strict_kernel"],
        "no_role_transfer_exported": not matrix["premises"]["role_transfer_theorem_available"],
    }
    return {
        "status": "P2851_AMPLITUDE_NORMALIZATION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2850": sha(P2850)},
        "amplitude_normalization_bridge_atom_audit": {
            "input_statuses_rechecked": {"P2850": p2850.get("status")},
            "parameters": {
                "alpha_geo_float": ALPHA_GEO,
                "alpha_geo_symbolic": "4*ln(2)",
                "legacy_omega": "pi/4",
                "legacy_phi": "pi/6",
                "legacy_beta_tors": LEGACY_BETA_TORS,
                "strict_omega": STRICT_OMEGA,
                "strict_phi": STRICT_PHI,
                "strict_beta": STRICT_BETA,
                "strict_eta": STRICT_ETA,
                "distances": list(DISTANCES),
            },
            "amplitude_rows": rows,
            "ratio_stats": rstats,
            "least_squares_constant_amplitude_fit": fit,
            "two_point_constant_amplitude_matrix": two_point,
            "required_amplitude_bridge_premises": list(REQUIRED_AMPLITUDE_BRIDGE_PREMISES),
            "premise_matrix": matrix,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_amplitude_normalization_bridge_obstruction_audit": all(facts.values()),
            "exports_amplitude_normalization_bridge_atom": False,
        },
        "decision": {
            "negative_export_flags": {
                "amplitude_normalization_bridge_atom_exported": False,
                "alpha_geo_strict_source_law_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "selector_closure_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2851 tests exactly one bridge atom: amplitude/normalization passage from legacy alpha_geo=4 ln 2 to the unit-amplitude strict gate kernel.  Pointwise legacy/strict amplitude ratios are not constant, two-point ratio checks reject an amplitude-only bridge, and the best constant-amplitude least-squares fit leaves a nonzero residual.  Thus alpha_geo cannot be silently absorbed into K_strict_gate without the missing phase/frequency and damping/compression bridge atoms plus a role-safe alpha_geo source law.",
            "next_honest_step": "Do not claim amplitude passage, full bridge, role transfer, L_total, EOM, Hamiltonian, or ToE closure from alpha_geo rescaling.  The next admissible move is a combined bridge-obligation reconciliation matrix over the now-tested atoms (damping/compression, amplitude, EML syntax) to identify whether any single remaining typed source atom is still live; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, object]) -> None:
    audit = payload["amplitude_normalization_bridge_atom_audit"]
    fit = audit["least_squares_constant_amplitude_fit"]
    rstats = audit["ratio_stats"]
    lines = [
        "# P2851/S1801 amplitude-normalization kernel bridge-atom audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Ratio stats",
        f"- min_ratio={rstats['min_ratio']}",
        f"- max_ratio={rstats['max_ratio']}",
        f"- max_minus_min={rstats['max_minus_min']}",
        f"- sign_changes_present={rstats['sign_changes_present']}",
        f"- constant_ratio={rstats['constant_ratio']}",
        "",
        "## Least-squares constant-amplitude fit",
        f"- best_fit_amplitude={fit['best_fit_amplitude']}",
        f"- sse={fit['sse']}",
        f"- max_abs_residual={fit['max_abs_residual']}",
        "",
        "## Two-point constant-amplitude matrix",
        audit["two_point_constant_amplitude_matrix"]["proof_statement"],
    ]
    for row in audit["two_point_constant_amplitude_matrix"]["rows"]:
        lines.append(f"- pair={row['pair']}: ratio_a={row['ratio_a']}; ratio_b={row['ratio_b']}; same={row['same_constant_amplitude']}")
    lines.extend([
        "",
        "## Premise matrix",
        f"- accepted_as_amplitude_normalization_bridge_atom={audit['premise_matrix']['accepted_as_amplitude_normalization_bridge_atom']}",
        f"- missing_premises={audit['premise_matrix']['missing_premises']}",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, object]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2850))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2851/S1801 amplitude-normalization kernel bridge-atom audit", "## P2851/S1801 amplitude-normalization kernel bridge-atom audit\n\n`P2851/S1801` tests the amplitude/normalization bridge atom from legacy `alpha_geo=4 ln 2` to the unit-amplitude `K_strict_gate`.  Pointwise legacy/strict amplitude ratios are not constant, audited two-point ratio checks reject an amplitude-only bridge, and the best constant-amplitude least-squares fit leaves a nonzero residual.  Therefore `alpha_geo` cannot be silently absorbed into `K_strict_gate` without missing phase/frequency and damping/compression bridge atoms plus a role-safe strict source law.  No amplitude bridge, full kernel bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2851/S1801 amplitude bridge Ltotal guard", "## P2851/S1801 amplitude bridge `L_total` guard\n\n`P2851/S1801` adds no action term.  Constant amplitude rescaling of `K_strict_gate` does not provide a role-safe `alpha_geo` source law, unit-bearing source density, nonproxy `L_total` insertion, EOM, Hamiltonian, or role-bearing term.\n")
    append_once(AGENTS, "Current amplitude-normalization kernel bridge-atom guardrail (P2851/S1801, 2026-06-18)", "## Current amplitude-normalization kernel bridge-atom guardrail (P2851/S1801, 2026-06-18)\n\n- P2851 tests exactly one bridge atom: legacy `alpha_geo=4 ln 2` amplitude passage into the unit-amplitude strict gate kernel.\n- Pointwise legacy/strict amplitude ratios are not constant; two-point ratio checks reject an amplitude-only bridge, and the best constant-amplitude fit has nonzero residual.\n- Do not promote `alpha_geo` rescaling to amplitude bridge, full kernel bridge, role transfer, selector closure, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n- A next admissible move is a combined bridge-obligation reconciliation matrix over damping/compression, amplitude, and EML syntax to identify whether any single remaining typed source atom is still live; otherwise preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
