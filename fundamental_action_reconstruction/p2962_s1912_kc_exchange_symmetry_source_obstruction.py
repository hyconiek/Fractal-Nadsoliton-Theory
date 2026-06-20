#!/usr/bin/env python3
"""P2962/S1912: K/C exchange-symmetry source obstruction.

P2961 proved a clean quotient-fixed candidate, but its missing premise was that
the nadsoliton exports the K/C exchange symmetry.  P2962 tests that premise on
the actual P2938 ingredients, rather than on the abstract coefficient pair
(a,b).  The audit asks whether the kernel-excess vector K and the
character-negativity vector C are presently interchangeable as typed artifacts.

The finite answer is negative: K and C have different support cardinalities,
different nonzero coordinate multisets, different total weights, and different
provenance labels.  Therefore P2961's exchange is a useful developmental
symmetry of the coefficient quotient, but it is not yet a strict artifact
symmetry of the P2938 ingredients.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go import K, C, TARGET
from p2961_s1911_exchange_balanced_scale_quotient_source_candidate import OUT as P2961

OUT = GEN / "p2962_s1912_kc_exchange_symmetry_source_obstruction.json"
MD = GEN / "p2962_s1912_kc_exchange_symmetry_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LABELS = ["p1", "p2", "p3", "p4", "p5"]


def support(v: list[int]) -> list[int]:
    return [i for i, x in enumerate(v) if x != 0]


def permutation_maps_source_to_target(source: list[int], target: list[int]) -> list[list[int]]:
    maps = []
    for perm in itertools.permutations(range(len(source))):
        image = [0] * len(source)
        for i, j in enumerate(perm):
            image[j] = source[i]
        if image == target:
            maps.append(list(perm))
    return maps


def invariant_rows() -> list[dict[str, Any]]:
    k_supp, c_supp = support(K), support(C)
    return [
        {"invariant": "support_cardinality", "K": len(k_supp), "C": len(c_supp), "matches": len(k_supp) == len(c_supp)},
        {"invariant": "nonzero_coordinate_multiset", "K": sorted(K[i] for i in k_supp), "C": sorted(C[i] for i in c_supp), "matches": sorted(K[i] for i in k_supp) == sorted(C[i] for i in c_supp)},
        {"invariant": "total_weight", "K": sum(K), "C": sum(C), "matches": sum(K) == sum(C)},
        {"invariant": "typed_provenance_label", "K": "kernel_excess", "C": "character_negativity", "matches": False},
        {"invariant": "coordinate_permutation_to_counterpart", "K_to_C_maps": len(permutation_maps_source_to_target(K, C)), "C_to_K_maps": len(permutation_maps_source_to_target(C, K)), "matches": bool(permutation_maps_source_to_target(K, C) and permutation_maps_source_to_target(C, K))},
    ]


def obstruction_rows(invariants: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"obligation": "P2961_coefficient_quotient_exchange_candidate_available", "satisfied": True, "evidence": "P2961 exports the primitive coefficient-quotient fixed-point candidate"},
        {"obligation": "K_and_C_support_cardinalities_match", "satisfied": next(r for r in invariants if r["invariant"] == "support_cardinality")["matches"], "evidence": "K has support size 2 while C has support size 3"},
        {"obligation": "K_and_C_nonzero_multisets_match", "satisfied": next(r for r in invariants if r["invariant"] == "nonzero_coordinate_multiset")["matches"], "evidence": "K nonzero multiset [1,2] differs from C [2,2,2]"},
        {"obligation": "K_and_C_typed_provenance_equivalence_exported", "satisfied": False, "evidence": "current artifacts label K as kernel-excess and C as character-negativity, not as two typed copies of one object"},
        {"obligation": "strict_nadsoliton_artifact_symmetry_source_exported", "satisfied": False, "evidence": "no current theorem exports an internal nadsoliton automorphism swapping the two typed ingredients"},
        {"obligation": "unit_bearing_nonproxy_coupling_exported", "satisfied": False, "evidence": "no action-density coupling receiving the exchange-selected vector is constructed"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["coefficient_quotient_candidate", "support_cardinalities_match", "nonzero_multisets_match", "typed_provenance_equivalence", "strict_artifact_symmetry_source", "unit_bearing_nonproxy_coupling"]
    rows = []
    for mask in range(1 << len(names)):
        present = {name: bool(mask & (1 << i)) for i, name in enumerate(names)}
        rows.append({"mask": mask, "present": present, "accepts_strict_KC_exchange_source": all(present.values())})
    return rows


def build_payload(p2961: dict[str, Any]) -> dict[str, Any]:
    invariants = invariant_rows()
    obligations = obstruction_rows(invariants)
    matrix = acceptance_matrix()
    return {
        "status": "P2962_KC_EXCHANGE_SYMMETRY_SOURCE_OBSTRUCTION_NO_STRICT_EXPORT",
        "input_hashes": {"P2961": hashlib.sha256(P2961.read_bytes()).hexdigest() if P2961.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "KCExchangeSymmetry_StrictSourceObstruction",
            "artifact_vectors": {"K_kernel_excess": K, "C_character_negativity": C, "K_plus_C_target": TARGET, "coordinate_labels": LABELS},
            "artifact_invariant_rows": invariants,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "obstruction_certificate": {
            "coefficient_quotient_candidate_available": obligations[0]["satisfied"],
            "support_cardinality_match": obligations[1]["satisfied"],
            "nonzero_multiset_match": obligations[2]["satisfied"],
            "typed_provenance_equivalence_exported": obligations[3]["satisfied"],
            "strict_artifact_symmetry_source_exported": obligations[4]["satisfied"],
            "unit_bearing_nonproxy_coupling_exported": obligations[5]["satisfied"],
            "strict_KC_exchange_source_exported": all(row["satisfied"] for row in obligations),
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_strict_KC_exchange_source"]),
        },
        "decision": {
            "positive_progress": "P2962 sharpens P2961: the quotient exchange is mathematically clean, but the actual P2938 K and C artifacts are not currently typed-isomorphic.",
            "breakthrough": "No strict exchange source yet.  The obstruction is concrete: support sizes 2 vs 3, nonzero multisets [1,2] vs [2,2,2], and distinct provenance labels block an artifact-level K/C symmetry on current data.",
            "negative_export_flags": {k: False for k in ["strict_KC_exchange_source_exported", "strict_p2938_provenance_exported", "strict_ratio_package_source_exported", "damping_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay coefficient-quotient exchange, target_sum=9 cuts, primitive equal-summand convention, K+C decompositions, or beta-scale normalization.  A next proof-grade move must introduce a new typed mediator/functor that makes K and C comparable without erasing their support/provenance mismatch, or construct a unit-bearing nonproxy coupling that does not require artifact-level K/C exchange symmetry; otherwise pivot outside the ratio-package lane while preserving the P2929-P2962 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["obstruction_certificate"]
    lines = [
        "# P2962/S1912 K/C exchange-symmetry source obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Obstruction certificate",
        f"- coefficient quotient candidate available: `{cert['coefficient_quotient_candidate_available']}`",
        f"- support cardinality match: `{cert['support_cardinality_match']}`",
        f"- nonzero multiset match: `{cert['nonzero_multiset_match']}`",
        f"- typed provenance equivalence exported: `{cert['typed_provenance_equivalence_exported']}`",
        f"- strict artifact symmetry source exported: `{cert['strict_artifact_symmetry_source_exported']}`",
        f"- unit-bearing nonproxy coupling exported: `{cert['unit_bearing_nonproxy_coupling_exported']}`",
        f"- strict K/C exchange source exported: `{cert['strict_KC_exchange_source_exported']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
        "",
        "## Lay summary",
        payload["decision"]["positive_progress"],
        payload["decision"]["breakthrough"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2961))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2962/S1912 K/C exchange-symmetry source obstruction", "## P2962/S1912 K/C exchange-symmetry source obstruction\n\n`P2962/S1912` audits the missing strict source premise behind P2961.  The P2961 exchange is clean on coefficient quotient classes, but the actual P2938 ingredients are not currently interchangeable typed artifacts: `K=[1,2,0,0,0]` has support size `2` and nonzero multiset `[1,2]`, while `C=[0,0,2,2,2]` has support size `3` and nonzero multiset `[2,2,2]`; their labels are also kernel-excess versus character-negativity.  No coordinate permutation maps K to C.  Therefore no strict nadsoliton K/C exchange-symmetry source, strict P2938 provenance, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2962/S1912 K/C exchange-symmetry `L_total` guard", "## P2962/S1912 K/C exchange-symmetry `L_total` guard\n\n`P2962/S1912` shows that P2961's coefficient-level exchange quotient is not yet an artifact-level K/C symmetry of the actual P2938 ingredients.  Because support shape, nonzero coordinate multiset, and typed provenance labels differ, no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE from this exchange route without a new typed mediator/functor or an independent unit-bearing nonproxy coupling.\n")
    append_once(AGENTS, "Current K/C exchange-symmetry source obstruction guardrail (P2962/S1912, 2026-06-20)", "## Current K/C exchange-symmetry source obstruction guardrail (P2962/S1912, 2026-06-20)\n\n- P2962 audits the missing strict source behind P2961 and finds that the coefficient-level exchange quotient is not yet an artifact-level K/C symmetry.\n- The obstruction is finite and concrete: K has support size `2` and nonzero multiset `[1,2]`, C has support size `3` and nonzero multiset `[2,2,2]`, the provenance labels differ, and no coordinate permutation maps K to C.\n- Do not promote P2961/P2962 to strict K/C exchange source, strict P2938 provenance, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must introduce a new typed mediator/functor making K and C comparable without erasing their support/provenance mismatch, construct a unit-bearing nonproxy coupling that does not require artifact-level K/C exchange symmetry, or pivot outside the ratio-package lane while preserving the P2929-P2962 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
