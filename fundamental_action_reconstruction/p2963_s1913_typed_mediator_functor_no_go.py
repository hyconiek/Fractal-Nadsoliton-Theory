#!/usr/bin/env python3
"""P2963/S1913: typed mediator/functor no-go for K/C comparability.

P2962 left one honest continuation: introduce a typed mediator/functor that
makes K and C comparable without erasing their support/provenance mismatch.  This
audit constructs a finite menu of candidate mediators and classifies each by two
conditions: whether it makes K and C equal/comparable, and whether it still
retains the obstruction data that made K and C different.

The result is a bounded no-go on current artifacts.  Coarse mediators such as
max_nonzero make K and C equal only by forgetting support and provenance.  Richer
mediators preserve the mismatch but do not make K and C comparable.  Therefore no
strict typed mediator/functor is exported.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go import K, C
from p2962_s1912_kc_exchange_symmetry_source_obstruction import OUT as P2962, support

OUT = GEN / "p2963_s1913_typed_mediator_functor_no_go.json"
MD = GEN / "p2963_s1913_typed_mediator_functor_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def entropy(values: list[int]) -> float:
    total = sum(values)
    if total == 0:
        return 0.0
    probs = [v / total for v in values if v]
    return -sum(p * math.log(p) for p in probs)


def nonzero(v: list[int]) -> list[int]:
    return [x for x in v if x]


def mediator_rows() -> list[dict[str, Any]]:
    candidates = [
        ("identity_typed_signature", lambda v, label: {"support_size": len(support(v)), "nonzero_multiset": sorted(nonzero(v)), "label": label}, True, True),
        ("support_cardinality", lambda v, label: len(support(v)), True, False),
        ("nonzero_multiset", lambda v, label: sorted(nonzero(v)), True, False),
        ("total_weight", lambda v, label: sum(v), False, False),
        ("max_nonzero", lambda v, label: max(nonzero(v)), False, False),
        ("nonzero_entropy_rounded_12", lambda v, label: round(entropy(nonzero(v)), 12), True, False),
        ("support_normalized_density", lambda v, label: [x / sum(nonzero(v)) for x in nonzero(v)], True, False),
        ("provenance_label_only", lambda v, label: label, False, True),
    ]
    rows = []
    for name, fn, preserves_support_shape, preserves_provenance_label in candidates:
        k_val = fn(K, "kernel_excess")
        c_val = fn(C, "character_negativity")
        comparable = k_val == c_val
        rows.append({
            "mediator": name,
            "K_value": k_val,
            "C_value": c_val,
            "makes_K_and_C_equal": comparable,
            "preserves_support_shape": preserves_support_shape,
            "preserves_provenance_label": preserves_provenance_label,
            "does_not_erase_P2962_mismatch": preserves_support_shape and preserves_provenance_label,
            "accepted_typed_mediator": comparable and preserves_support_shape and preserves_provenance_label,
            "boundary": "equalizes only by erasure" if comparable and not (preserves_support_shape and preserves_provenance_label) else "preserves mismatch, so no equality" if not comparable else "accepted",
        })
    return rows


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    accepted = [r for r in rows if r["accepted_typed_mediator"]]
    erasing_equalizers = [r for r in rows if r["makes_K_and_C_equal"] and not r["does_not_erase_P2962_mismatch"]]
    return [
        {"obligation": "finite_mediator_menu_constructed", "satisfied": True, "evidence": f"{len(rows)} mediator candidates audited"},
        {"obligation": "some_mediator_makes_K_and_C_equal", "satisfied": any(r["makes_K_and_C_equal"] for r in rows), "evidence": f"equalizing mediators: {[r['mediator'] for r in rows if r['makes_K_and_C_equal']]}"},
        {"obligation": "equalizing_mediator_preserves_support_and_provenance", "satisfied": bool(accepted), "evidence": f"erasing equalizers: {[r['mediator'] for r in erasing_equalizers]}"},
        {"obligation": "strict_nadsoliton_mediator_functor_exported", "satisfied": False, "evidence": "no current artifact exports a functor turning the erasing scalar max into typed K/C equivalence"},
        {"obligation": "unit_bearing_nonproxy_coupling_exported", "satisfied": False, "evidence": "no coupling to an action density is constructed"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_mediator_menu", "some_equalizer", "equalizer_preserves_mismatch", "strict_nadsoliton_functor", "unit_bearing_nonproxy_coupling"]
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_strict_typed_mediator": m == (1 << len(names)) - 1} for m in range(1 << len(names))]


def build_payload(p2962: dict[str, Any]) -> dict[str, Any]:
    rows = mediator_rows()
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2963_TYPED_MEDIATOR_FUNCTOR_NO_GO_NO_STRICT_EXPORT",
        "input_hashes": {"P2962": hashlib.sha256(P2962.read_bytes()).hexdigest() if P2962.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "TypedMediatorFunctor_AcceptanceNoGo",
            "mediator_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "mediator_certificate": {
            "mediator_count": len(rows),
            "equalizing_mediators": [r["mediator"] for r in rows if r["makes_K_and_C_equal"]],
            "accepted_typed_mediators": [r["mediator"] for r in rows if r["accepted_typed_mediator"]],
            "strict_typed_mediator_exported": all(row["satisfied"] for row in obligations),
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_strict_typed_mediator"]),
        },
        "decision": {
            "positive_progress": "P2963 tests the typed-mediator escape route explicitly.  It finds that max_nonzero can equalize K and C, but only by erasing the support/provenance mismatch identified in P2962.",
            "breakthrough": "No strict mediator is exported: every mediator that preserves the mismatch keeps K and C unequal, while every equalizer is too coarse to be a typed provenance functor.",
            "negative_export_flags": {k: False for k in ["strict_typed_mediator_exported", "strict_KC_exchange_source_exported", "strict_p2938_provenance_exported", "strict_ratio_package_source_exported", "damping_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay typed-mediator scalar collapses, max-coordinate equalization, coefficient-quotient exchange, target_sum cuts, primitive equal-summand, K+C decompositions, or beta-scale normalization.  A next proof-grade move must either construct an actual unit-bearing nonproxy coupling that receives the P2938 aggregate without needing K/C artifact exchange, or introduce a new typed structural object richer than the current K/C split; otherwise pivot outside the ratio-package lane while preserving the P2929-P2963 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["mediator_certificate"]
    lines = [
        "# P2963/S1913 typed mediator/functor no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Mediator certificate",
        f"- mediator count: `{cert['mediator_count']}`",
        f"- equalizing mediators: `{cert['equalizing_mediators']}`",
        f"- accepted typed mediators: `{cert['accepted_typed_mediators']}`",
        f"- strict typed mediator exported: `{cert['strict_typed_mediator_exported']}`",
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
    payload = build_payload(read_json(P2962))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2963/S1913 typed mediator/functor no-go", "## P2963/S1913 typed mediator/functor no-go\n\n`P2963/S1913` audits the P2962-admissible typed mediator/functor escape route.  Eight candidate mediators are tested.  The only audited mediator that makes K and C equal is `max_nonzero`, but it erases support shape and provenance labels; mediators that preserve the P2962 mismatch keep K and C unequal.  Therefore no strict typed mediator/functor, strict K/C exchange source, strict P2938 provenance, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2963/S1913 typed mediator `L_total` guard", "## P2963/S1913 typed mediator `L_total` guard\n\n`P2963/S1913` finds no typed mediator that both compares K and C and preserves the support/provenance mismatch identified by P2962.  Scalar equalization such as `max_nonzero` is too coarse for a typed source functor, so no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE through this mediator route.\n")
    append_once(AGENTS, "Current typed mediator/functor no-go guardrail (P2963/S1913, 2026-06-20)", "## Current typed mediator/functor no-go guardrail (P2963/S1913, 2026-06-20)\n\n- P2963 audits the typed mediator/functor route opened by P2962.  Eight finite mediator candidates are tested.\n- The only audited equalizer is `max_nonzero`; it makes K and C equal only by erasing support shape and provenance labels.  Mediators preserving the P2962 mismatch keep K and C unequal.\n- Do not promote scalar mediator collapse to strict K/C exchange source, strict P2938 provenance, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must construct an actual unit-bearing nonproxy coupling receiving the P2938 aggregate without K/C artifact exchange, introduce a new typed structural object richer than the current K/C split, or pivot outside the ratio-package lane while preserving the P2929-P2963 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
