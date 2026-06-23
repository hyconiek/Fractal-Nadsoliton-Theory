#!/usr/bin/env python3
"""P3035/S1985: Z12 directional branch-selector source obstruction.

Attack exactly one P3028 foundation atom: strict selector/branch source.  This
constructs a finite orientation torsor (+direction/-direction) for sampled
K_strict_gate on Z12 and tests computable branch-selection candidates.  The
candidates are real finite receivers, but they either give identical scores to
both directions or choose only a chart/label convention, so no non-premise strict
selector/branch source is exported.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice import OUT as P3028
from p3034_s1984_z12_finite_difference_action_eom_unit_obstruction import OUT as P3034, formal_action

OUT = GEN / "p3035_s1985_z12_directional_branch_selector_source_obstruction.json"
MD = GEN / "p3035_s1985_z12_directional_branch_selector_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS = [1, 5, 7, 11]


def kernel_vector() -> list[float]:
    return [k_strict(label) for label in range(1, N + 1)]


def directional_gradient(values: list[float], direction: int) -> list[float]:
    return [values[(i + direction) % N] - values[i] for i in range(N)]


def directional_action(values: list[float], direction: int) -> float:
    grad = directional_gradient(values, direction)
    return 0.5 * sum(g * g for g in grad)


def signed_flux(values: list[float], direction: int) -> float:
    # Odd-looking nearest-neighbour flux receiver.  On a cycle with scalar K
    # samples it telescopes to an orientation-blind value, so it is not a source.
    return sum(values[i] * (values[(i + direction) % N] - values[(i - direction) % N]) for i in range(N))


def relabel(values: list[float], unit: int) -> list[float]:
    # U(12) relabel on labels 1..12, returned in the same label order.
    by_label = {label: values[label - 1] for label in range(1, N + 1)}
    return [by_label[((unit * label - 1) % N) + 1] for label in range(1, N + 1)]


def max_label_anchor(values: list[float]) -> int:
    return max(range(1, N + 1), key=lambda label: (abs(values[label - 1]), -label))


def build_matrix() -> dict[str, Any]:
    values = kernel_vector()
    plus = 1
    minus = -1
    plus_action = directional_action(values, plus)
    minus_action = directional_action(values, minus)
    plus_flux = signed_flux(values, plus)
    minus_flux = signed_flux(values, minus)
    candidate_rows = [
        {
            "candidate": "nearest_neighbor_dirichlet_direction_selector",
            "plus_score": round(plus_action, 12),
            "minus_score": round(minus_action, 12),
            "orientation_separating": abs(plus_action - minus_action) > 1e-12,
            "accepted_as_strict_branch_source": False,
            "failure": "forward and backward cyclic Dirichlet energies are exactly equal",
        },
        {
            "candidate": "formal_action_receiver_direction_selector",
            "plus_score": round(formal_action(values), 12),
            "minus_score": round(formal_action(list(reversed(values))), 12),
            "orientation_separating": abs(formal_action(values) - formal_action(list(reversed(values)))) > 1e-12,
            "accepted_as_strict_branch_source": False,
            "failure": "P3034 formal action is reflection/orientation blind and already lacks unit-bearing provenance",
        },
        {
            "candidate": "nearest_neighbor_signed_flux_selector",
            "plus_score": round(plus_flux, 12),
            "minus_score": round(minus_flux, 12),
            "orientation_separating": abs(plus_flux - minus_flux) > 1e-12,
            "accepted_as_strict_branch_source": False,
            "failure": "scalar cyclic flux telescopes to zero/equal opposite-direction scores",
        },
        {
            "candidate": "largest_kernel_peak_label_anchor",
            "anchor_label": max_label_anchor(values),
            "unit_relabel_anchor_labels": {str(u): max_label_anchor(relabel(values, u)) for u in UNITS},
            "orientation_separating": False,
            "accepted_as_strict_branch_source": False,
            "failure": "peak label is a chart anchor, not an Aut(Z12)-compatible directed branch source",
        },
    ]
    obligations = [
        {"obligation": "attacks_single_P3028_foundation_atom", "satisfied": True, "detail": "only strict selector/branch source is tested"},
        {"obligation": "explicit_orientation_torsor", "satisfied": True, "detail": "the two branches are +direction and -direction on cyclic Z12"},
        {"obligation": "finite_candidate_receivers_computable", "satisfied": True, "detail": "Dirichlet, action, flux, and peak-anchor receivers are computed"},
        {"obligation": "orientation_score_separation", "satisfied": any(row["orientation_separating"] for row in candidate_rows), "detail": "no tested receiver separates +direction from -direction without convention"},
        {"obligation": "aut_z12_compatible_nonpremise_source", "satisfied": False, "detail": "inversion units exchange directions and no strict source breaks the torsor"},
        {"obligation": "chart_independent_branch_localizer", "satisfied": False, "detail": "label anchors move under U(12) relabeling and are not physical sectors"},
        {"obligation": "coupling_to_classical_readout_rows", "satisfied": False, "detail": "no theorem couples the branch sign to spacetime/time/matter/energy/observer rows"},
    ]
    return {
        "object": "Z12DirectionalBranchSelectorSource_ObstructionMatrix",
        "tested_foundation_atom": "strict_selector_or_branch_source",
        "orientation_torsor": {"branches": ["+direction", "-direction"], "inversion_units": [7, 11], "aut_units": UNITS},
        "candidate_rows": candidate_rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "candidate_rows": len(candidate_rows),
            "orientation_separating_rows": sum(1 for row in candidate_rows if row["orientation_separating"]),
            "accepted_branch_source_rows": sum(1 for row in candidate_rows if row["accepted_as_strict_branch_source"]),
            "inversion_units_exchange_branches": 2,
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_selector_branch_source_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_strict_branch_source"] for row in candidate_rows),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3028)
    read_json(P3034)
    matrix = build_matrix()
    return {
        "status": "P3035_Z12_DIRECTIONAL_BRANCH_SELECTOR_SOURCE_OBSTRUCTION_NO_SELECTOR_EXPORT",
        "input_hashes": {
            "P3028": hashlib.sha256(P3028.read_bytes()).hexdigest() if P3028.exists() else None,
            "P3034": hashlib.sha256(P3034.read_bytes()).hexdigest() if P3034.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "A finite +direction/-direction orientation torsor for sampled K_strict_gate on Z12 is explicit, and four selector candidates are computable.  Dirichlet/action/flux receivers are orientation-blind, while the largest-peak anchor is chart-label dependent under U(12) relabeling.  Thus no non-premise Aut-compatible strict selector/branch source is exported.",
            "negative_export_flags": {k: False for k in ["strict_selector_branch_source_exported", "qw2191_discharged", "classical_transition_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay orientation-blind quadratic/action/flux receivers or chart-label peak anchors as strict selector sources.  A next move must supply a genuinely new non-premise inversion-odd/chiral branch source with a fixed sign and an explicit coupling to one P3028 readout row, or pivot to the remaining external physical unit-source atom.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3035/S1985 Z12 directional branch-selector source obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- candidate rows: `{c['candidate_rows']}`",
        f"- orientation-separating rows: `{c['orientation_separating_rows']}`",
        f"- accepted branch-source rows: `{c['accepted_branch_source_rows']}`",
        f"- inversion units exchanging branches: `{c['inversion_units_exchange_branches']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict selector/branch source exported: `{c['strict_selector_branch_source_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3035/S1985 Z12 directional branch-selector source obstruction", "## P3035/S1985 Z12 directional branch-selector source obstruction\n\n`P3035/S1985` attacks exactly one P3028 foundation atom: strict selector/branch source.  It constructs the finite `+direction/-direction` orientation torsor on cyclic `Z12` for sampled `K_strict_gate` and tests four computable selector candidates: nearest-neighbor Dirichlet direction, P3034 formal action under reflection, nearest-neighbor signed flux, and largest-kernel-peak label anchor.  The first three are orientation-blind, and the last is a chart-label anchor rather than an `Aut(Z12)`-compatible source.  Hence `0/4` candidates export a strict branch source; no `QW-2191` discharge, classical transition, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3035/S1985 branch-selector source `L_total` guard", "## P3035/S1985 branch-selector source `L_total` guard\n\n`P3035/S1985` adds no physical `L_total` term.  The directional/action/flux receivers are branch-blind, and the peak-label anchor is chart-conventional, so no strict variational branch source or readout-sector coupling is installed.\n")
    append_once(AGENTS, "Current Z12 directional branch-selector source guardrail (P3035/S1985, 2026-06-23)", "## Current Z12 directional branch-selector source guardrail (P3035/S1985, 2026-06-23)\n\n- P3035 attacks exactly one P3028 foundation atom: strict selector/branch source.\n- The finite `+direction/-direction` orientation torsor is explicit, but Dirichlet/action/flux receivers are orientation-blind and peak-label anchors are chart-dependent under `U(12)` relabeling.\n- Do not promote orientation-blind quadratic/action/flux receivers, largest-peak label anchors, or internal branch conventions to strict selector closure, `QW-2191` discharge, classical transition, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must supply a genuinely new non-premise inversion-odd/chiral branch source with fixed sign and explicit readout coupling, or pivot to the remaining external physical unit-source atom.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
