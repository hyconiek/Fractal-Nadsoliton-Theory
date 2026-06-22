#!/usr/bin/env python3
"""P3011/S1961: ToE strict-kernel/selector role-separation matrix.

After P3010 closed the cube-map-functional-graph lane as bounded no-go, this
constructs a new typed explanatory object answering a current conceptual
question: if the theory is a ToE candidate, what does its strict kernel describe,
and what is the selector?

The audit deliberately does not claim ToE closure.  It separates roles:
K_strict_gate describes an operational gate/compression correlation profile of
nadsoliton structure over separation d, while a selector is a symmetry-breaking
choice rule/source that chooses a directed representative or phase/orientation
branch from a quotient/torsor that the kernel values alone do not orient.
"""
from __future__ import annotations

import hashlib, json, math
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3010_s1960_cube_map_action_installation_obstruction import OUT as P3010

OUT = GEN / "p3011_s1961_toe_strict_kernel_selector_role_separation_matrix.json"
MD = GEN / "p3011_s1961_toe_strict_kernel_selector_role_separation_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8
MODULUS = 12
UNITS = [1, 5, 7, 11]
DIRECTED_UNITS = [1, 11]


def strict_kernel(d: int) -> float:
    return math.cos(OMEGA * d + PHI) / (1 + BETA * (d ** ETA))


def kernel_profile() -> dict[str, Any]:
    vals = [strict_kernel(d) for d in range(MODULUS + 1)]
    return {
        "formula": "K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "parameters": {"omega": OMEGA, "phi": PHI, "beta": BETA, "eta": ETA},
        "domain_sample": list(range(MODULUS + 1)),
        "values_rounded_12dp": [round(v, 12) for v in vals],
        "tail_ratio_abs_K12_over_K1": round(abs(vals[12]) / abs(vals[1]), 12),
        "nonlinear_denominator_values_rounded_12dp": [round(1 + d ** ETA, 12) for d in range(MODULUS + 1)],
        "describes": [
            "operational gate-selected correlation/compression profile over separation d",
            "cosine phase/frequency channel plus nonlinear heavy-tail damping/compression",
            "strict working successor of the legacy intermediate kernel, not silent legacy role transfer",
        ],
    }


def selector_torsor_witness() -> dict[str, Any]:
    action_rows = []
    for u, s in product(UNITS, DIRECTED_UNITS):
        image = (u * s) % MODULUS
        action_rows.append({"unit": u, "source_directed_unit": s, "image": image})
    preserving = sorted(u for u in UNITS if (u * 1) % MODULUS == 1 and (u * 11) % MODULUS == 11)
    reversing = sorted(u for u in UNITS if (u * 1) % MODULUS == 11 and (u * 11) % MODULUS == 1)
    invariant_choices = [s for s in DIRECTED_UNITS if all((u * s) % MODULUS == s for u in UNITS)]
    projective_orbit = sorted({1, 11})
    return {
        "selector_problem": "choose a directed representative/phase-orientation branch from the projective pair {+1,-1} in Z/12Z",
        "directed_pair_mod12": DIRECTED_UNITS,
        "aut_z12_units": UNITS,
        "action_rows": action_rows,
        "orientation_preserving_units_on_pair": preserving,
        "orientation_reversing_units_on_pair": reversing,
        "aut_invariant_directed_choices": invariant_choices,
        "projective_choice_available": projective_orbit,
        "directed_selector_available_without_extra_source": bool(invariant_choices),
        "interpretation": "The quotient/projective pair exists, but Aut(Z12) contains the pair-reversing unit 11 (and other units outside the pair), so no nonpremise directed representative is selected by carrier symmetry alone.",
    }


def role_separation_rows(k: dict[str, Any], s: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "row": "ontology",
            "answers": "what exists first",
            "current_content": "nadsoliton itself as primordial information in a solitonic state",
            "satisfied": True,
            "nonpromotion_boundary": "not a separate lower information layer",
        },
        {
            "row": "strict_kernel",
            "answers": "how strict nadsoliton correlations/compression are encoded as a gate-selected working profile",
            "current_content": k["formula"],
            "satisfied": True,
            "nonpromotion_boundary": "does not by itself install SM/GR coefficients, L_total terms, role transfer, or a selector source",
        },
        {
            "row": "selector",
            "answers": "which representative/orientation/phase branch is chosen when quotient/torsor data leave a sign ambiguity",
            "current_content": s["selector_problem"],
            "satisfied": s["directed_selector_available_without_extra_source"],
            "nonpromotion_boundary": "projective availability is not directed strict-core selector closure",
        },
        {
            "row": "toe_closure",
            "answers": "whether the kernel plus selector plus unit-bearing Lagrangian/EOM generate complete physical observables",
            "current_content": "requires strict kernel-to-coefficient map, selector source, unit-bearing L_total/EOM/Hamiltonian, observable generator, and bridge/role discipline",
            "satisfied": False,
            "nonpromotion_boundary": "the current matrix is explanatory/proof-obligation bookkeeping, not ToE closure",
        },
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["kernel_profile", "kernel_to_coefficients", "directed_selector_source", "unit_bearing_ltotal", "eom_hamiltonian", "observable_generator", "bridge_role_audit"]
    rows = []
    for bits in product([False, True], repeat=len(names)):
        present = dict(zip(names, bits))
        rows.append({"present": present, "accepts_toe_strict_closure": all(bits)})
    return rows


def build_payload(p3010_path: Any) -> dict[str, Any]:
    k = kernel_profile()
    s = selector_torsor_witness()
    rows = role_separation_rows(k, s)
    matrix = acceptance_matrix()
    return {
        "status": "P3011_TOE_STRICT_KERNEL_SELECTOR_ROLE_SEPARATION_MATRIX_NO_CLOSURE",
        "input_hashes": {"P3010": hashlib.sha256(p3010_path.read_bytes()).hexdigest() if p3010_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "ToEStrictKernelSelectorRoleSeparation_ObligationMatrix",
            "strict_kernel_profile": k,
            "selector_torsor_witness": s,
            "role_separation_rows": rows,
            "finite_acceptance_matrix": matrix,
        },
        "answer_certificate": {
            "strict_kernel_describes": k["describes"],
            "selector_is": s["selector_problem"],
            "orientation_reversing_units_on_pair": s["orientation_reversing_units_on_pair"],
            "aut_invariant_directed_choice_count": len(s["aut_invariant_directed_choices"]),
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_toe_strict_closure"]),
        },
        "decision": {
            "plain_answer": "If the program is read as a ToE candidate, the strict kernel describes the nadsoliton's operational gate-selected correlation/compression law; the selector is not that kernel, but the missing or explicitly supplied symmetry-breaking rule/source choosing a directed phase/orientation representative from quotient/torsor ambiguity.",
            "breakthrough": "Bounded explanatory progress: the matrix gives a finite role separation and an Aut(Z12) witness that projective ± data do not select a directed sign without an extra source.",
            "negative_export_flags": {k: False for k in ["toe_closure_exported", "strict_selector_closure_exported", "qw2191_discharged", "kernel_to_coefficient_map_exported", "unit_bearing_ltotal_exported", "eom_hamiltonian_exported", "observable_generator_exported", "bridge_closure_exported", "role_transfer_exported"]},
            "next_honest_step": "Attack exactly one missing ToE-closure atom exposed by P3011, preferably a strict kernel-to-coefficient provenance map with units/sign/variational normalization, or a genuinely new directed selector source; do not replay cube-map, selector, bridge, role-transfer, or L_total closure lanes without the named new object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["answer_certificate"]
    lines = [
        "# P3011/S1961 ToE strict-kernel/selector role-separation matrix", "",
        f"Status: `{payload['status']}`", "",
        "## Direct answer",
        payload["decision"]["plain_answer"], "",
        "## Finite certificate",
        f"- strict kernel describes: `{cert['strict_kernel_describes']}`",
        f"- selector is: `{cert['selector_is']}`",
        f"- orientation-reversing units on `{{+1,-1}}`: `{cert['orientation_reversing_units_on_pair']}`",
        f"- Aut-invariant directed choices: `{cert['aut_invariant_directed_choice_count']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3010)
    payload = build_payload(P3010)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3011/S1961 ToE strict-kernel/selector role-separation matrix", "## P3011/S1961 ToE strict-kernel/selector role-separation matrix\n\n`P3011/S1961` answers the current conceptual ToE question without promoting closure.  In the honest current reading, `K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)` describes the strict working gate/correlation/compression profile of the nadsoliton over separation `d`: a cosine phase/frequency channel plus nonlinear heavy-tail damping/compression.  The selector is a different kind of object: a symmetry-breaking rule/source that chooses a directed phase/orientation representative from quotient/torsor ambiguity.  The finite `Aut(Z12)` witness on the projective pair `{+1,-1}` has orientation-reversing unit `11` and zero Aut-invariant directed choices, so projective carrier data alone do not export selector closure.  This matrix exports no `QW-2191` discharge, kernel-to-coefficient map, unit-bearing `L_total`, EOM/Hamiltonian, bridge closure, role transfer, observable generator, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3011/S1961 strict-kernel/selector role-separation `L_total` guard", "## P3011/S1961 strict-kernel/selector role-separation `L_total` guard\n\n`P3011/S1961` is an explanatory obligation matrix, not a new variational term.  It separates the strict kernel profile from the selector role: kernel values encode gate/correlation/compression over `d`, while a selector must source a directed phase/orientation representative.  The finite `Aut(Z12)` pair witness has no invariant directed choice, so no selector source, unit-bearing density, EOM/Hamiltonian term, bridge closure, role transfer, or ToE closure may be installed in `L_total` from this role separation alone.\n")
    append_once(AGENTS, "Current ToE strict-kernel/selector role-separation guardrail (P3011/S1961, 2026-06-22)", "## Current ToE strict-kernel/selector role-separation guardrail (P3011/S1961, 2026-06-22)\n\n- P3011 answers the conceptual ToE question by separating roles: `K_strict_gate` describes the strict working gate/correlation/compression profile of the nadsoliton over separation `d`; the selector is a separate symmetry-breaking rule/source choosing a directed phase/orientation representative from quotient/torsor ambiguity.\n- The finite `Aut(Z12)` witness on `{+1,-1}` has orientation-reversing unit `11` and zero Aut-invariant directed choices; projective carrier availability is not directed selector closure.\n- Do not promote this role-separation matrix to `QW-2191` discharge, strict selector closure, kernel-to-coefficient provenance, unit-bearing `L_total`, EOM/Hamiltonian, bridge closure, role transfer, observable generation, or ToE closure.\n- The next honest proof-grade move should attack exactly one missing ToE-closure atom, preferably strict kernel-to-coefficient provenance with units/sign/variational normalization or a genuinely new directed selector source, without replaying cube-map, selector, bridge, role-transfer, or `L_total` lanes.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
