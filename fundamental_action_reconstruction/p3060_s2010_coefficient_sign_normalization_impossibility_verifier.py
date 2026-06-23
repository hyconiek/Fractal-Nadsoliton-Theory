#!/usr/bin/env python3
"""P3060/S2010: coefficient-sign normalization impossibility verifier.

P3059 left one precise missing object: a strict non-premise global
coefficient-sign normalization/source rule for the polarity-odd synthesis
module.  P3060 does not add another clue.  It constructs a larger named class
of candidate normalizers and proves, by finite enumeration, that every
currently admissible equivariant/magnitude/support rule remains sign-even and
therefore cannot choose between a coefficient vector and its negative.
"""
from __future__ import annotations

import hashlib, itertools, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3059_s2009_polarity_odd_source_law_synthesis_verifier import OUT as P3059, enumerate_synthesis_module

OUT = GEN / "p3060_s2010_coefficient_sign_normalization_impossibility_verifier.json"
MD = GEN / "p3060_s2010_coefficient_sign_normalization_impossibility_verifier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "global_coefficient_sign_normalization": r"global coefficient-sign|coefficient-sign normalization|global coefficient sign|sign-normalization",
    "sign_paired_source_law_obstruction": r"sign-pair|sign-paired source law|paired polarity|coefficient.*negative",
    "nonpremise_normalizer_impossibility": r"non-premise.*normalization|normalization.*impossible|impossibility.*coefficient|larger named coefficient/source class",
}

WEIGHTS = [-2, -1, 0, 1, 2]
EVEN_FEATURES = ("l1_norm", "l2_norm", "support_size", "max_abs", "abs_sum_mod2")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def features(coeffs: tuple[int, ...]) -> dict[str, int]:
    abs_values = [abs(c) for c in coeffs]
    return {
        "l1_norm": sum(abs_values),
        "l2_norm": sum(c * c for c in coeffs),
        "support_size": sum(1 for c in coeffs if c),
        "max_abs": max(abs_values),
        "abs_sum_mod2": sum(abs_values) % 2,
    }


def linear_score(coeffs: tuple[int, ...], weights: tuple[int, ...]) -> int:
    feat = features(coeffs)
    return sum(weights[i] * feat[name] for i, name in enumerate(EVEN_FEATURES))


def enumerate_normalizers() -> list[dict[str, Any]]:
    rows = []
    nonzero_coeffs = [tuple(row["coefficient_tuple"]) for row in enumerate_synthesis_module() if row["nonzero_signed_value"]]
    for weights in itertools.product(WEIGHTS, repeat=len(EVEN_FEATURES)):
        if all(w == 0 for w in weights):
            continue
        separating_pairs = 0
        unique_positive_choices = 0
        for coeffs in nonzero_coeffs:
            neg = tuple(-c for c in coeffs)
            s_pos = linear_score(coeffs, weights)
            s_neg = linear_score(neg, weights)
            if s_pos != s_neg:
                separating_pairs += 1
            if s_pos > s_neg:
                unique_positive_choices += 1
        rows.append({
            "weights": dict(zip(EVEN_FEATURES, weights)),
            "score_family": "integer linear score on sign-even magnitude/support invariants",
            "tested_sign_pairs": len(nonzero_coeffs),
            "separating_sign_pairs": separating_pairs,
            "unique_positive_choices": unique_positive_choices,
            "accepted_as_nonpremise_sign_normalizer": False,
            "blocker": "all allowed features are invariant under c -> -c, so every score is identical on each sign pair",
        })
    return rows


def build_payload() -> dict[str, Any]:
    p3059 = read_json(P3059)
    grep_rows = content_grep()
    normalizers = enumerate_normalizers()
    total = len(normalizers)
    separating = sum(1 for row in normalizers if row["separating_sign_pairs"])
    obligations = [
        {"obligation": "content_first_grep_before_impossibility_test", "satisfied": True, "detail": "searched by coefficient-sign normalization and sign-pair obstruction content"},
        {"obligation": "construct_larger_named_candidate_class", "satisfied": True, "detail": "SignEvenMagnitudeSupportNormalizerClass over five invariants and weights [-2,2]^5"},
        {"obligation": "exhaust_candidate_normalizers", "satisfied": True, "detail": "all 3124 nonzero linear normalizer score rules are enumerated"},
        {"obligation": "find_nonpremise_global_sign_rule", "satisfied": False, "detail": "zero normalizers separate any c/-c sign pair"},
        {"obligation": "export_selector_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE closure follows"},
    ]
    return {
        "status": "P3060_COEFFICIENT_SIGN_NORMALIZATION_CLASS_IMPOSSIBILITY_NO_EXPORT",
        "input_hashes": {"P3059": hashlib.sha256(P3059.read_bytes()).hexdigest() if P3059.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "candidate_class": {
                "object": "SignEvenMagnitudeSupportNormalizerClass",
                "target_missing_object": "strict_nonpremise_global_coefficient_sign_normalization_source_rule",
                "features": list(EVEN_FEATURES),
                "weight_domain": "integer weights in [-2,2]^5 excluding the zero score",
                "acceptance_rule": "accept only a score rule that separates at least one c/-c pair and is exported by a strict non-premise source law",
            },
            "normalizer_row_summary": {
                "row_count": total,
                "all_rows_have_zero_separating_pairs": all(row["separating_sign_pairs"] == 0 for row in normalizers),
                "sample_rows": normalizers[:5],
                "row_digest": hashlib.sha256(json.dumps(normalizers, sort_keys=True).encode("utf-8")).hexdigest(),
            },
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "sign_even_features": len(EVEN_FEATURES),
            "weight_values_per_feature": len(WEIGHTS),
            "candidate_normalizers": total,
            "normalizers_separating_any_sign_pair": separating,
            "accepted_nonpremise_sign_normalizers": 0,
            "p3059_status_seen": p3059.get("status"),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3060 proves impossibility only inside the named SignEvenMagnitudeSupportNormalizerClass: all 3124 nonzero integer linear score rules built from current sign-even magnitude/support invariants are invariant under c -> -c.  Therefore 0 rules separate any P3059 sign pair and 0 export a strict non-premise global coefficient-sign normalization/source rule.  This closes the requested larger-class replay without exporting G_selector, QW-2191 discharge, selector closure, L_total, bridge/role transfer, or ToE closure.",
            "negative_export_flags": {k: False for k in ["global_coefficient_sign_normalization_exported", "strict_polarity_odd_source_law_boundary_condition_exported", "actual_g_selector_formula_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay sign-even magnitude/support normalizers.  The next proof-grade move must introduce one genuinely sign-odd, strict-sourced coefficient-sign normalizer with an exported non-premise signed value, or pivot to a different named P3057 atom and run the same content-first no-replay discipline while preserving no selector/readout/L_total/bridge/ToE export.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3060/S2010 coefficient-sign normalization impossibility verifier", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- sign-even features: `{c['sign_even_features']}`", f"- weight values per feature: `{c['weight_values_per_feature']}`", f"- candidate normalizers: `{c['candidate_normalizers']}`", f"- normalizers separating any sign pair: `{c['normalizers_separating_any_sign_pair']}`", f"- accepted nonpremise sign normalizers: `{c['accepted_nonpremise_sign_normalizers']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3060/S2010 coefficient-sign normalization impossibility verifier", "## P3060/S2010 coefficient-sign normalization impossibility verifier\n\n`P3060/S2010` responds to the P3059 recommendation by proving impossibility in one larger named class: `SignEvenMagnitudeSupportNormalizerClass`.  It enumerates all `3124` nonzero integer linear score rules over five sign-even magnitude/support invariants (`l1_norm`, `l2_norm`, `support_size`, `max_abs`, `abs_sum_mod2`) with weights in `[-2,2]^5`.  Because every feature is invariant under `c -> -c`, `0` normalizers separate any P3059 sign pair and `0` export a strict non-premise global coefficient-sign normalization/source rule.  No actual `G_selector`, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3060/S2010 coefficient-sign normalization impossibility `L_total` guard", "## P3060/S2010 coefficient-sign normalization impossibility `L_total` guard\n\n`P3060/S2010` adds no physical `L_total` term.  It is a finite obstruction inside a sign-even normalizer class; it exports no unit-bearing signed source/action/EOM carrier and no nonproxy variational chain rule.\n")
    append_once(AGENTS, "Current coefficient-sign normalization impossibility guardrail (P3060/S2010, 2026-06-23)", "## Current coefficient-sign normalization impossibility guardrail (P3060/S2010, 2026-06-23)\n\n- P3060 proves a bounded impossibility result only inside `SignEvenMagnitudeSupportNormalizerClass`, the larger named coefficient/source class built from current sign-even magnitude/support invariants.\n- The finite enumeration checks all `3124` nonzero integer linear normalizer scores over five invariants with weights in `[-2,2]^5`; `0` scores separate any `c/-c` sign pair and `0` export a non-premise global coefficient-sign source rule.\n- Do not promote sign-even magnitude/support scoring, coefficient conventions, or this bounded impossibility result to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must introduce a genuinely sign-odd strict-sourced coefficient-sign normalizer with an exported non-premise signed value, or pivot to a different named P3057 atom while preserving the P3048-P3060 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
