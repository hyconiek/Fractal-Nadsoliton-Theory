#!/usr/bin/env python3
"""P2849/S1799: damping/compression kernel bridge-atom audit.

After P2848, the finite vertex-density route is bounded no-go.  P2849 pivots to
one concrete kernel bridge atom, as recommended: the denominator passage from
legacy linear torsion damping

    1 + beta_tors * d

to strict nonlinear compression

    1 + beta * d**eta.

The audit is deliberately narrow.  It does not attempt a full bridge theorem or
role transfer.  It proves and computes only whether the legacy linear damping
law itself can source the strict nonlinear beta/eta atom without an additional
strict source premise.
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
P2848 = GEN / "p2848_s1798_coupling_coefficient_unit_source_law_audit.json"
OUT = GEN / "p2849_s1799_damping_compression_kernel_bridge_atom_audit.json"
MD = GEN / "p2849_s1799_damping_compression_kernel_bridge_atom_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LEGACY_BETA_TORS = Fraction(1, 100)
STRICT_BETA = Fraction(1, 1)
STRICT_ETA = Fraction(9, 5)
DISTANCES = tuple(range(1, 13))
REQUIRED_DAMPING_BRIDGE_PREMISES = (
    "legacy_linear_damping_defined",
    "strict_nonlinear_compression_defined",
    "two_point_exact_completion_map",
    "target_independent_positive_beta_source",
    "eta_source_law",
    "amplitude_normalization_compatibility",
    "phase_frequency_compatibility",
    "selector_or_topological_source_compatibility",
    "role_transfer_theorem_available",
)


def legacy_denominator(d: int) -> Fraction:
    return Fraction(1, 1) + LEGACY_BETA_TORS * d


def strict_denominator(d: int) -> float:
    return 1.0 + float(STRICT_BETA) * (d ** float(STRICT_ETA))


def beta_eff_for_eta(d: int, eta: Fraction = STRICT_ETA) -> float:
    """Beta that would make beta*d**eta equal beta_tors*d at this d."""
    return float(LEGACY_BETA_TORS) * (d ** (1.0 - float(eta)))


def exact_two_point_eta_from_legacy(a: int, b: int) -> Fraction | None:
    """For beta*d**eta = beta_tors*d at a,b, exact equality implies eta=1.

    The proof is algebraic: dividing equations gives (a/b)**eta = a/b.  For
    distinct positive a,b this forces eta=1 over positive real powers.
    """
    if a <= 0 or b <= 0 or a == b:
        return None
    return Fraction(1, 1)


def compression_table() -> list[dict[str, Any]]:
    rows = []
    for d in DISTANCES:
        legacy = legacy_denominator(d)
        strict = strict_denominator(d)
        beta_eff = beta_eff_for_eta(d)
        rows.append(
            {
                "d": d,
                "legacy_denominator": str(legacy),
                "strict_denominator_float": strict,
                "strict_over_legacy_float": strict / float(legacy),
                "beta_eff_for_strict_eta_to_match_legacy": beta_eff,
            }
        )
    return rows


def beta_eff_stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    vals = [row["beta_eff_for_strict_eta_to_match_legacy"] for row in rows]
    mean = sum(vals) / len(vals)
    variance = sum((v - mean) ** 2 for v in vals) / len(vals)
    return {
        "min": min(vals),
        "max": max(vals),
        "max_over_min": max(vals) / min(vals),
        "mean": mean,
        "variance": variance,
        "constant_across_distances": math.isclose(max(vals), min(vals), rel_tol=0.0, abs_tol=0.0),
    }


def two_point_matrix() -> dict[str, Any]:
    pairs = [(1, 2), (1, 3), (2, 3), (3, 6), (6, 12)]
    rows = []
    for a, b in pairs:
        eta = exact_two_point_eta_from_legacy(a, b)
        rows.append(
            {
                "pair": [a, b],
                "forced_eta_for_exact_legacy_linear_match": str(eta),
                "strict_eta": str(STRICT_ETA),
                "strict_eta_matches_forced_eta": eta == STRICT_ETA,
            }
        )
    return {
        "proof_statement": "For two distinct positive distances a,b, exact equality beta*d^eta = beta_tors*d at both points forces eta=1 and beta=beta_tors; therefore the strict eta=9/5 nonlinear compression cannot be sourced by legacy linear damping alone.",
        "rows": rows,
        "all_pairs_reject_strict_eta_as_legacy_linear_exact_match": all(not row["strict_eta_matches_forced_eta"] for row in rows),
    }


def premise_matrix(two_point: dict[str, Any]) -> dict[str, Any]:
    premises = {
        "legacy_linear_damping_defined": True,
        "strict_nonlinear_compression_defined": True,
        "two_point_exact_completion_map": not two_point["all_pairs_reject_strict_eta_as_legacy_linear_exact_match"],
        "target_independent_positive_beta_source": False,
        "eta_source_law": False,
        "amplitude_normalization_compatibility": False,
        "phase_frequency_compatibility": False,
        "selector_or_topological_source_compatibility": False,
        "role_transfer_theorem_available": False,
    }
    return {
        "premises": premises,
        "missing_premises": [key for key in REQUIRED_DAMPING_BRIDGE_PREMISES if not premises[key]],
        "accepted_as_damping_compression_bridge_atom": all(premises.values()),
    }


def build_payload(p2848: dict[str, Any]) -> dict[str, Any]:
    rows = compression_table()
    stats = beta_eff_stats(rows)
    two_point = two_point_matrix()
    matrix = premise_matrix(two_point)
    facts = {
        "p2848_rechecked": p2848.get("status") == "P2848_COUPLING_COEFFICIENT_UNIT_SOURCE_LAW_AUDIT_NO_CLOSURE",
        "legacy_and_strict_denominators_defined": matrix["premises"]["legacy_linear_damping_defined"] and matrix["premises"]["strict_nonlinear_compression_defined"],
        "two_point_exact_match_rejects_strict_eta": two_point["all_pairs_reject_strict_eta_as_legacy_linear_exact_match"],
        "effective_beta_varies_for_strict_eta": not stats["constant_across_distances"],
        "no_eta_source_law_exported": not matrix["premises"]["eta_source_law"],
        "no_target_independent_beta_source_exported": not matrix["premises"]["target_independent_positive_beta_source"],
        "no_role_transfer_exported": not matrix["premises"]["role_transfer_theorem_available"],
    }
    return {
        "status": "P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2848": sha(P2848)},
        "damping_compression_bridge_atom_audit": {
            "input_statuses_rechecked": {"P2848": p2848.get("status")},
            "parameters": {
                "legacy_beta_tors": str(LEGACY_BETA_TORS),
                "strict_beta": str(STRICT_BETA),
                "strict_eta": str(STRICT_ETA),
                "distances": list(DISTANCES),
            },
            "compression_rows": rows,
            "beta_eff_stats": stats,
            "two_point_exact_match_matrix": two_point,
            "required_damping_bridge_premises": list(REQUIRED_DAMPING_BRIDGE_PREMISES),
            "premise_matrix": matrix,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_damping_compression_bridge_obstruction_audit": all(facts.values()),
            "exports_damping_compression_bridge_atom": False,
        },
        "decision": {
            "negative_export_flags": {
                "damping_compression_bridge_atom_exported": False,
                "target_independent_positive_beta_source_exported": False,
                "eta_source_law_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "selector_closure_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2849 tests exactly one kernel bridge atom: legacy linear torsion damping to strict nonlinear beta/eta compression.  A two-point exact-match theorem forces eta=1 and beta=beta_tors for any exact legacy-linear source, while the strict value eta=9/5 fails every audited two-point row.  Matching legacy damping at strict eta would require beta_eff(d)=beta_tors*d^(1-eta), which varies across distances.  Current artifacts therefore do not export a target-independent beta source or eta source law.",
            "next_honest_step": "Do not claim kernel bridge, role transfer, L_total, EOM, Hamiltonian, or ToE closure from damping/compression similarity.  The next admissible move is either one new strict source law for eta and target-independent beta, or a different single kernel bridge atom such as amplitude normalization passage.  Without a new source premise, preserve the no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["damping_compression_bridge_atom_audit"]
    lines = [
        "# P2849/S1799 damping/compression kernel bridge-atom audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Parameters",
        f"- legacy_beta_tors={audit['parameters']['legacy_beta_tors']}",
        f"- strict_beta={audit['parameters']['strict_beta']}",
        f"- strict_eta={audit['parameters']['strict_eta']}",
        "",
        "## Two-point exact-match matrix",
        audit["two_point_exact_match_matrix"]["proof_statement"],
    ]
    for row in audit["two_point_exact_match_matrix"]["rows"]:
        lines.append(f"- pair={row['pair']}: forced_eta={row['forced_eta_for_exact_legacy_linear_match']}; strict_eta_matches={row['strict_eta_matches_forced_eta']}")
    stats = audit["beta_eff_stats"]
    lines.extend([
        "",
        "## beta_eff stats for strict eta",
        f"- min={stats['min']}",
        f"- max={stats['max']}",
        f"- max_over_min={stats['max_over_min']}",
        f"- constant_across_distances={stats['constant_across_distances']}",
        "",
        "## Premise matrix",
        f"- accepted_as_damping_compression_bridge_atom={audit['premise_matrix']['accepted_as_damping_compression_bridge_atom']}",
        f"- missing_premises={audit['premise_matrix']['missing_premises']}",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2848))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2849/S1799 damping/compression kernel bridge-atom audit", "## P2849/S1799 damping/compression kernel bridge-atom audit\n\n`P2849/S1799` pivots after the bounded finite-density route to one concrete kernel bridge atom: legacy linear torsion damping `1 + beta_tors*d` versus strict nonlinear compression `1 + beta*d^eta`.  The two-point exact-match theorem shows that sourcing strict compression from the legacy linear law alone forces `eta=1` and `beta=beta_tors`; the strict `eta=9/5` fails this exact bridge condition, and the effective beta needed to match legacy damping varies with distance.  No target-independent beta source, eta source law, damping/compression bridge atom, full kernel bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2849/S1799 damping bridge Ltotal guard", "## P2849/S1799 damping bridge `L_total` guard\n\n`P2849/S1799` adds no action term.  The damping/compression denominator comparison is a kernel-bridge obstruction audit only; without a target-independent beta source and eta source law, it does not provide a unit-bearing source density, nonproxy `L_total` insertion, EOM, Hamiltonian, or role-bearing term.\n")
    append_once(AGENTS, "Current damping/compression kernel bridge-atom guardrail (P2849/S1799, 2026-06-18)", "## Current damping/compression kernel bridge-atom guardrail (P2849/S1799, 2026-06-18)\n\n- P2849 tests exactly one concrete kernel bridge atom: legacy linear torsion damping `1 + beta_tors*d` versus strict nonlinear compression `1 + beta*d^eta`.\n- Exact two-point matching of the legacy linear source forces `eta=1` and `beta=beta_tors`; strict `eta=9/5` is not sourced by that law, and the effective beta required at strict eta varies with distance.\n- Do not promote damping/compression similarity to full kernel bridge, role transfer, selector closure, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n- A next admissible move requires one new strict source law for `eta` and target-independent positive `beta`, or a different single bridge atom such as amplitude normalization passage; otherwise preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
