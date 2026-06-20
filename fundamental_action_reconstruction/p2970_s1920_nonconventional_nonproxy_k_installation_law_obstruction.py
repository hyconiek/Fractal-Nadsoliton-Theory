#!/usr/bin/env python3
"""P2970/S1920: nonconventional nonproxy k-installation law obstruction.

P2969 left one ratio-package move that is not a replay of symbolic U or unit
conventions: construct a nonconventional nonproxy installation law that fixes the
exponent k in c_k=(9/5)U^k internally.  This audit attacks that exact object.

The finite computation reuses the P2968 exponent set and tests candidate
installation laws as score/selection rules on k.  Some rules select a unique k,
but only by mathematical convention or by replaying the dimensionless k=0
section; nonproxy variation remains k-blind without an independent unit/source
coupling.  Thus no strict installation law is exported.
"""
from __future__ import annotations

import hashlib, json
from fractions import Fraction
from typing import Any, Callable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2968_s1918_coefficient_source_law_exponent_blind_obstruction import OUT as P2968, exponent_rows
from p2967_s1917_dimension_field_selector_source_obstruction import parse_frac
from p2969_s1919_strict_unit_U_source_candidate_obstruction import OUT as P2969

OUT = GEN / "p2970_s1920_nonconventional_nonproxy_k_installation_law_obstruction.json"
MD = GEN / "p2970_s1920_nonconventional_nonproxy_k_installation_law_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def p2968_exponents() -> list[Fraction]:
    return [parse_frac(row["exponent_k"]) for row in exponent_rows()]


def selection_summary(name: str, exponents: list[Fraction], score: Callable[[Fraction], tuple[Any, ...]]) -> dict[str, Any]:
    scored = [(k, score(k)) for k in exponents]
    best_score = min(v for _, v in scored)
    selected = [k for k, v in scored if v == best_score]
    return {
        "selector": name,
        "best_score": [str(x) for x in best_score],
        "selected_k": [f"{k.numerator}/{k.denominator}" for k in selected],
        "selected_count": len(selected),
    }


def candidate_installation_rows(exponents: list[Fraction]) -> list[dict[str, Any]]:
    selectors = {
        "minimal_absolute_exponent_installation": lambda k: (abs(k),),
        "dimensionless_k_zero_installation": lambda k: (0 if k == 0 else 1, abs(k)),
        "minimal_positive_exponent_installation": lambda k: (0 if k > 0 else 1, abs(k)),
        "maximal_negative_exponent_installation": lambda k: (0 if k < 0 else 1, abs(k)),
        "integer_exponent_subfamily_installation": lambda k: (0 if k.denominator == 1 else 1, abs(k)),
        "nonproxy_euler_stationarity_placeholder": lambda k: (0,),
    }
    summaries = {name: selection_summary(name, exponents, fn) for name, fn in selectors.items()}
    rows = [
        {
            "candidate": "minimal_absolute_exponent_installation",
            "selected_k": summaries["minimal_absolute_exponent_installation"]["selected_k"],
            "unique_k": summaries["minimal_absolute_exponent_installation"]["selected_count"] == 1,
            "nonconventional": False,
            "nonproxy_action_density_coupled": False,
            "replays_blocked_lane": True,
            "strict_source_exported": False,
            "current_artifact_available": True,
            "witness": "selects k=0 by size convention, replaying the dimensionless section",
        },
        {
            "candidate": "dimensionless_k_zero_installation",
            "selected_k": summaries["dimensionless_k_zero_installation"]["selected_k"],
            "unique_k": summaries["dimensionless_k_zero_installation"]["selected_count"] == 1,
            "nonconventional": False,
            "nonproxy_action_density_coupled": False,
            "replays_blocked_lane": True,
            "strict_source_exported": False,
            "current_artifact_available": True,
            "witness": "explicitly replays the P2967/P2968 k=0 obstruction",
        },
        {
            "candidate": "minimal_positive_exponent_installation",
            "selected_k": summaries["minimal_positive_exponent_installation"]["selected_k"],
            "unique_k": summaries["minimal_positive_exponent_installation"]["selected_count"] == 1,
            "nonconventional": False,
            "nonproxy_action_density_coupled": False,
            "replays_blocked_lane": False,
            "strict_source_exported": False,
            "current_artifact_available": True,
            "witness": "chooses the smallest positive unit power by ordering convention, not by source",
        },
        {
            "candidate": "maximal_negative_exponent_installation",
            "selected_k": summaries["maximal_negative_exponent_installation"]["selected_k"],
            "unique_k": summaries["maximal_negative_exponent_installation"]["selected_count"] == 1,
            "nonconventional": False,
            "nonproxy_action_density_coupled": False,
            "replays_blocked_lane": False,
            "strict_source_exported": False,
            "current_artifact_available": True,
            "witness": "chooses the closest negative unit power by ordering convention, not by source",
        },
        {
            "candidate": "integer_exponent_subfamily_installation",
            "selected_k": summaries["integer_exponent_subfamily_installation"]["selected_k"],
            "unique_k": summaries["integer_exponent_subfamily_installation"]["selected_count"] == 1,
            "nonconventional": False,
            "nonproxy_action_density_coupled": False,
            "replays_blocked_lane": True,
            "strict_source_exported": False,
            "current_artifact_available": True,
            "witness": "still selects k=0 first under the audited subfamily score",
        },
        {
            "candidate": "nonproxy_euler_stationarity_placeholder",
            "selected_k": summaries["nonproxy_euler_stationarity_placeholder"]["selected_k"],
            "unique_k": summaries["nonproxy_euler_stationarity_placeholder"]["selected_count"] == 1,
            "nonconventional": True,
            "nonproxy_action_density_coupled": False,
            "replays_blocked_lane": False,
            "strict_source_exported": False,
            "current_artifact_available": True,
            "witness": "formal Euler insertion is k-blind without an independent unit/source coupling",
        },
        {
            "candidate": "completed_strict_k_installation_law_schema",
            "selected_k": [],
            "unique_k": True,
            "nonconventional": True,
            "nonproxy_action_density_coupled": True,
            "replays_blocked_lane": False,
            "strict_source_exported": True,
            "current_artifact_available": False,
            "witness": "schema only; no current theorem internally fixes k and installs c_k into P2964",
        },
    ]
    for row in rows:
        row["accepted_current_strict_k_installation"] = row["current_artifact_available"] and row["unique_k"] and row["nonconventional"] and row["nonproxy_action_density_coupled"] and not row["replays_blocked_lane"] and row["strict_source_exported"]
    return rows


def obligation_rows(rows: list[dict[str, Any]], exponents: list[Fraction]) -> list[dict[str, Any]]:
    return [
        {"obligation": "P2968_exponent_set_reused", "satisfied": len(exponents) == 24, "evidence": f"{len(exponents)} distinct k values audited"},
        {"obligation": "unique_k_candidate_exists", "satisfied": any(r["unique_k"] and r["current_artifact_available"] for r in rows), "evidence": "several ordering predicates choose one k"},
        {"obligation": "nonconventional_current_candidate_exists", "satisfied": any(r["nonconventional"] and r["current_artifact_available"] for r in rows), "evidence": "Euler placeholder is nonconventional in form but k-blind"},
        {"obligation": "nonproxy_action_density_coupled", "satisfied": any(r["nonproxy_action_density_coupled"] and r["current_artifact_available"] for r in rows), "evidence": "no current row couples the selected k into P2964 action density"},
        {"obligation": "no_blocked_lane_replay", "satisfied": any((not r["replays_blocked_lane"]) and r["current_artifact_available"] for r in rows), "evidence": "positive/negative ordering avoids explicit k=0 replay but remains conventional"},
        {"obligation": "accepted_current_strict_k_installation", "satisfied": any(r["accepted_current_strict_k_installation"] for r in rows), "evidence": "completed strict k-installation law is unavailable"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["P2968_exponents", "unique_k", "nonconventional", "nonproxy_density_coupling", "strict_source", "no_unit_convention", "no_k0_replay"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_strict_k_installation_law": m == full} for m in range(1 << len(names))]


def build_payload(p2968_path: Any, p2969_path: Any) -> dict[str, Any]:
    exponents = p2968_exponents()
    rows = candidate_installation_rows(exponents)
    obligations = obligation_rows(rows, exponents)
    matrix = acceptance_matrix()
    return {
        "status": "P2970_NONCONVENTIONAL_NONPROXY_K_INSTALLATION_LAW_OBSTRUCTION_NO_STRICT_EXPORT",
        "input_hashes": {
            "P2968": hashlib.sha256(p2968_path.read_bytes()).hexdigest() if p2968_path.exists() else None,
            "P2969": hashlib.sha256(p2969_path.read_bytes()).hexdigest() if p2969_path.exists() else None,
        },
        "constructed_theoretical_objects": {"candidate_object": "NonconventionalNonproxyKInstallationLaw_Obstruction", "exponent_set": [f"{k.numerator}/{k.denominator}" for k in exponents], "installation_law_rows": rows, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "installation_certificate": {"distinct_exponent_count": len(exponents), "candidate_count": len(rows), "accepted_current_strict_installations": [r["candidate"] for r in rows if r["accepted_current_strict_k_installation"]], "acceptance_matrix_rows": len(matrix), "accepted_rows": sum(1 for r in matrix if r["accepts_strict_k_installation_law"])},
        "decision": {
            "positive_progress": "P2970 converts the remaining nonconventional installation-law route into a finite k-selection and coupling-obligation matrix over the 24 P2968 exponents.",
            "breakthrough": "No strict k-installation law is exported: unique k choices are conventional or k=0 replay, while the nonconventional Euler placeholder is k-blind and lacks P2964 nonproxy coupling.",
            "negative_export_flags": {k: False for k in ["strict_k_installation_law_exported", "strict_unit_U_source_exported", "strict_coefficient_source_law_exported", "strict_unit_bearing_nonproxy_coupling_exported", "strict_ratio_package_source_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay minimal-|k|, k=0, positive/negative ordering conventions, integer-subfamily scores, symbolic U, U=1, beta normalization, entropy/reference-cell units, Gamma_9_5 import, or scalar Euler placeholders.  The next proof-grade move must introduce a genuinely new typed structural object outside the current ratio-package lane, or preserve the P2929-P2970 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["installation_certificate"]
    lines = ["# P2970/S1920 nonconventional nonproxy k-installation law obstruction", "", f"Status: `{payload['status']}`", "", "## Installation certificate", f"- distinct exponent count: `{cert['distinct_exponent_count']}`", f"- candidate count: `{cert['candidate_count']}`", f"- accepted current strict installations: `{cert['accepted_current_strict_installations']}`", f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2968)
    read_json(P2969)
    payload = build_payload(P2968, P2969)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2970/S1920 nonconventional nonproxy k-installation law obstruction", "## P2970/S1920 nonconventional nonproxy k-installation law obstruction\n\n`P2970/S1920` attacks the P2969-admissible route: a nonconventional nonproxy installation law internally fixing `k` in `c_k=(9/5)U^k`.  The finite audit reuses the 24 P2968 exponents and tests minimal-|k|, dimensionless `k=0`, smallest positive, closest negative, integer-subfamily, and Euler-placeholder installation candidates.  Unique selections are conventional or replay `k=0`; the nonconventional Euler placeholder remains k-blind and uncoupled from P2964.  Therefore no strict k-installation law, unit `U` source, coefficient source law, unit-bearing nonproxy coupling, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2970/S1920 k-installation `L_total` guard", "## P2970/S1920 k-installation `L_total` guard\n\n`P2970/S1920` shows that finite k-selection rules for `c_k=(9/5)U^k` either choose `k` by convention/replay or remain k-blind without a P2964 coupling theorem.  Since no nonconventional nonproxy installation law fixes `k`, no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current nonconventional nonproxy k-installation obstruction guardrail (P2970/S1920, 2026-06-20)", "## Current nonconventional nonproxy k-installation obstruction guardrail (P2970/S1920, 2026-06-20)\n\n- P2970 audits the last ratio-package route left by P2969: a nonconventional nonproxy installation law internally fixing `k` in `c_k=(9/5)U^k`.\n- Minimal-|k|, `k=0`, positive/negative ordering, and integer-subfamily predicates either replay blocked lanes or select by convention; the Euler-placeholder row is nonconventional in form but k-blind and uncoupled from P2964.\n- Do not promote k-selection ordering predicates, symbolic `U`, `U=1`, beta normalization, entropy/reference-cell units, `Gamma_9_5`, exponent-blind `9/5`, k=0 replay, or scalar Euler placeholders to strict ratio-package source, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must introduce a genuinely new typed structural object outside the current ratio-package lane, or preserve the P2929-P2970 no-strict-export boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
