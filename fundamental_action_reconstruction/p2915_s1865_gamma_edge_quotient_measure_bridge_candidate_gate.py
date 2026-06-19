#!/usr/bin/env python3
"""P2915/S1865: Gamma edge-quotient measure bridge candidate gate.

P2914 proved an exact normalization mismatch for the finite Lambda/Gamma measure
model: site normalization gives 12*m=1 while directed-edge normalization gives
144*m=1.  P2915 constructs the honest missing object suggested by that no-go: a
quotient/renormalization bridge that divides the 144 directed edge carriers into
12 equal fibers before integration.

The finite calculation shows that several quotient maps resolve the arithmetic
mismatch equally well.  Source-site, target-site, and displacement quotients all
have 12 fibers of size 12 and convert site-normalized m=1/12 into a total
quotient-edge measure of 1.  This is a real bridge-candidate family, but not a
strict theorem: current artifacts do not choose one quotient, prove that it is
the continuum/nonproxy integration measure, or couple it to field-variable
provenance and Gamma_9_5 sourcehood.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from typing import Any, Callable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2914 = GEN / "p2914_s1864_gamma_continuum_measure_normalization_obstruction.json"
OUT = GEN / "p2915_s1865_gamma_edge_quotient_measure_bridge_candidate_gate.json"
MD = GEN / "p2915_s1865_gamma_edge_quotient_measure_bridge_candidate_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12
EDGE_COUNT = N * N


def f(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def edges() -> list[tuple[int, int]]:
    return [(i, j) for i in range(N) for j in range(N)]


def quotient_maps() -> dict[str, Callable[[tuple[int, int]], int]]:
    return {
        "source_site_quotient_pi_src(i,j)=i": lambda edge: edge[0],
        "target_site_quotient_pi_tgt(i,j)=j": lambda edge: edge[1],
        "displacement_quotient_pi_diff(i,j)=j-i_mod12": lambda edge: (edge[1] - edge[0]) % N,
    }


def fiber_summary(name: str, fn: Callable[[tuple[int, int]], int], site_weight: Fraction) -> dict[str, Any]:
    fibers: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for edge in edges():
        fibers[fn(edge)].append(edge)
    fiber_sizes = sorted(len(items) for items in fibers.values())
    quotient_edge_weight = site_weight / N
    quotient_total = EDGE_COUNT * quotient_edge_weight
    return {
        "quotient_name": name,
        "fiber_count": len(fibers),
        "fiber_sizes": fiber_sizes,
        "all_fibers_size_12": all(size == N for size in fiber_sizes),
        "site_normalized_m": f(site_weight),
        "per_edge_quotient_weight": f(quotient_edge_weight),
        "total_quotient_edge_measure": f(quotient_total),
        "normalization_mismatch_resolved_arithmetically": quotient_total == 1,
        "sample_fibers": {
            str(label): [[i, j] for i, j in fibers[label][:4]]
            for label in sorted(fibers)[:3]
        },
        "strictly_selected": False,
    }


def build_payload(p2914: dict[str, Any]) -> dict[str, Any]:
    site_weight = Fraction(1, N)
    quotient_rows = [fiber_summary(name, fn, site_weight) for name, fn in quotient_maps().items()]
    arithmetic_passes = [row for row in quotient_rows if row["normalization_mismatch_resolved_arithmetically"]]
    return {
        "status": "P2915_GAMMA_EDGE_QUOTIENT_MEASURE_BRIDGE_CANDIDATE_GATE_NO_EXPORT",
        "input_hashes": {"P2914": hashlib.sha256(P2914.read_bytes()).hexdigest() if P2914.exists() else None},
        "constructed_theoretical_objects": {
            "missing_theorem_name": "Strict_Edge_Quotient_Measure_Bridge_Theorem",
            "quotient_bridge_schema": [
                "choose a strict quotient map pi: Z12 x Z12 -> 12 integration fibers",
                "renormalize each directed-edge contribution by the fiber size 12",
                "prove the quotient is the nonproxy continuum integration measure",
                "couple the quotient measure to P2911 Lambda and P2912 Jacobian field variables",
                "couple the resulting integral to a strict nonzero Gamma_9_5 source theorem",
            ],
            "quotient_candidate_rows": quotient_rows,
            "acceptance_obligation_rows": [
                {"obligation": "finite quotient candidates constructed", "passed": True, "exported_strictly": False},
                {"obligation": "12 equal fibers verified for each candidate", "passed": all(row["all_fibers_size_12"] for row in quotient_rows), "exported_strictly": False},
                {"obligation": "site-normalized measure renormalizes to quotient-edge total one", "passed": len(arithmetic_passes) == len(quotient_rows), "exported_strictly": False},
                {"obligation": "unique strict quotient selected", "passed": False, "exported_strictly": False},
                {"obligation": "strict field-variable/continuum integration theorem", "passed": False, "exported_strictly": False},
                {"obligation": "strict Gamma_9_5 source theorem", "passed": False, "exported_strictly": False},
            ],
        },
        "acceptance_matrix": {
            "p2914_rechecked_measure_obstruction": p2914.get("acceptance_matrix", {}).get("common_site_and_edge_normalization_solution_exists") is False,
            "site_count": N,
            "directed_edge_count": EDGE_COUNT,
            "quotient_candidate_count": len(quotient_rows),
            "arithmetic_normalization_resolving_candidate_count": len(arithmetic_passes),
            "strictly_selected_quotient_count": sum(1 for row in quotient_rows if row["strictly_selected"]),
            "site_normalized_m": f(site_weight),
            "per_edge_quotient_weight": f(site_weight / N),
            "quotient_edge_total": f(EDGE_COUNT * site_weight / N),
            "finite_quotient_bridge_candidate_family_constructed": True,
            "unique_strict_quotient_exported": False,
            "strict_continuum_integration_theorem_exported": False,
            "strict_gamma_9_5_source_exported": False,
            "accepted_as_nonproxy_ltotal_measure_bridge": False,
        },
        "decision": {
            "positive_witnesses": {
                "quotient_bridge_schema_constructed": True,
                "three_equal_fiber_quotients_constructed": True,
                "normalization_mismatch_arithmetically_resolved_by_each_candidate": len(arithmetic_passes) == len(quotient_rows),
            },
            "negative_export_flags": {
                "unique_strict_quotient_exported": False,
                "strict_continuum_integration_theorem_exported": False,
                "strict_field_variable_provenance_exported": False,
                "strict_gamma_9_5_source_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2915 constructs a real quotient-measure bridge candidate family for the P2914 mismatch.  Source-site, target-site, and displacement quotients each split 144 directed edges into 12 fibers of size 12 and renormalize site-normalized m=1/12 to total quotient-edge measure 1.  Because all three pass the finite arithmetic gate, the quotient is not uniquely strict-selected; no continuum integration theorem, field-variable provenance theorem, Gamma_9_5 source theorem, or nonproxy L_total export follows.",
            "next_honest_step": "The next proof-grade move must prove exactly one selector/provenance theorem for the quotient bridge: either a strict reason to choose source-site, target-site, displacement, or another specified quotient as the continuum integration measure, together with its field-variable coupling; or pivot back to a strict nonzero Gamma_9_5 action-unit source theorem.  Without that strict quotient/source theorem, preserve no-new-live-frontier for this lane.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2915/S1865 Gamma edge-quotient measure bridge candidate gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Quotient bridge gate",
        f"- quotient candidates: `{acc['quotient_candidate_count']}`",
        f"- arithmetic normalization-resolving candidates: `{acc['arithmetic_normalization_resolving_candidate_count']}`",
        f"- strictly selected quotients: `{acc['strictly_selected_quotient_count']}`",
        f"- site-normalized m: `{acc['site_normalized_m']}`",
        f"- per-edge quotient weight: `{acc['per_edge_quotient_weight']}`",
        f"- quotient edge total: `{acc['quotient_edge_total']}`",
        f"- accepted as nonproxy L_total measure bridge: `{acc['accepted_as_nonproxy_ltotal_measure_bridge']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2914))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2915/S1865 Gamma edge-quotient measure bridge candidate gate", "## P2915/S1865 Gamma edge-quotient measure bridge candidate gate\n\n`P2915/S1865` constructs a finite quotient-measure bridge candidate family for the P2914 `12` vs `144` mismatch.  Source-site, target-site, and displacement quotients each partition the `144` directed edges into `12` fibers of size `12`; renormalizing by the fiber size sends site-normalized `m=1/12` to quotient-edge total `1`.  This resolves the arithmetic mismatch only as readiness: all three candidates pass, so no unique strict quotient, continuum integration theorem, field-variable provenance theorem, strict `Gamma_9_5` source, nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2915/S1865 quotient-measure bridge `L_total` guard", "## P2915/S1865 quotient-measure bridge `L_total` guard\n\n`P2915/S1865` shows that quotienting the `144` directed edges by a `12`-element fiber can arithmetically repair the P2914 normalization mismatch, but source-site, target-site, and displacement quotients all work.  Without a strict theorem selecting one quotient as the continuum integration measure and coupling it to field variables and `Gamma_9_5`, no nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure follows.\n")
    append_once(AGENTS, "Current Gamma edge-quotient measure bridge candidate guardrail (P2915/S1865, 2026-06-19)", "## Current Gamma edge-quotient measure bridge candidate guardrail (P2915/S1865, 2026-06-19)\n\n- P2915 constructs three quotient-measure bridge candidates for the P2914 normalization mismatch: source-site, target-site, and displacement quotients, each with `12` fibers of size `12`.\n- Each candidate arithmetically repairs the site/edge total mismatch, but none is uniquely strict-selected or proven to be the continuum/nonproxy integration measure.\n- Do not promote the quotient candidates, P2911 pullback readiness, P2912 Jacobian readiness, P2913 Gamma source candidates, or P2914 measure obstruction to nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without a strict quotient/source theorem and field-variable coupling.\n- A next admissible move must prove exactly one quotient selector/provenance theorem or pivot to a strict nonzero `Gamma_9_5` action-unit source theorem; otherwise preserve no-new-live-frontier for the Gamma/Lambda lane.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
