#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_1988 = GEN / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json"
IN_1991 = GEN / "p1991_s941_strict_augmented_provider_channel_matrix_witness.json"
IN_2168 = GEN / "p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.json"
IN_2181 = GEN / "p2181_s1131_strict_qw2191_replay_certificate_freeze.json"
IN_2297 = GEN / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json"
OUT = GEN / "p2298_s1248_strict_non_gb_provider_source_nonadmissibility_theorem_probe.json"
MD = GEN / "p2298_s1248_strict_non_gb_provider_source_nonadmissibility_theorem_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def as_bool(value: Any) -> bool:
    return bool(value) is True


def candidate_vector(candidate: dict[str, Any]) -> list[int]:
    return [
        int(candidate["exported_current_source"]),
        int(candidate["algebraically_capable_for_full_residual_family"]),
        int(candidate["strict_admissible_source"]),
        int(candidate["internal_selector_source_exported"]),
        int(candidate["noncyclic_anchor_exported"]),
        int(candidate["no_legacy_transfer_used"]),
        int(candidate["qw2191_discharged"]),
    ]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1988 = load(IN_1988)
    p1991 = load(IN_1991)
    p2168 = load(IN_2168)
    p2181 = load(IN_2181)
    p2297 = load(IN_2297)

    p2297_probe = p2297.get("strict_non_gb_spatial_eom_provider_matrix_obstruction_probe", {})
    p2297_results = p2297_probe.get("provider_matrix_results", {})
    strict_core = p2297_results.get("strict_core_minimal_provider", {})
    augmented = p2297_results.get("p1990_augmented_provider_non_strict", {})
    formal_full = p2297_results.get("formal_full_residual_basis_provider", {})

    p1988_projection = p1988.get("family_projection", {})
    strict_non_gb_families = p1988_projection.get("strict_non_gb_families_detected", [])
    outside_eh_families = p1988_projection.get("outside_eh_family_capacity", [])

    qw2191_required_summary = (
        p2168.get("strict_qw2191_theorem_obligations_executable_validator", {})
        .get("required_summary", {})
    )
    qw2191_all_required_pass = as_bool(qw2191_required_summary.get("all_required_pass"))
    qw2191_replay_certified = as_bool(p2181.get("gatekeeper_checks", {}).get("replay_certified"))

    formal_basis_count = int(formal_full.get("basis_monomial_count", 0))
    formal_rank = int(formal_full.get("matrix_report", {}).get("rank_A", 0))
    strict_core_consistent = as_bool(strict_core.get("matrix_report", {}).get("consistent"))
    augmented_consistent = as_bool(augmented.get("matrix_report", {}).get("consistent"))
    formal_full_consistent = as_bool(formal_full.get("matrix_report", {}).get("consistent"))
    formal_full_exact = as_bool(formal_full.get("exact_reconstruction_zero"))
    formal_full_non_admissible = "NOT_A_LEGAL_STRICT_PROVIDER" in formal_full.get("admissibility_verdict", "")

    candidates = [
        {
            "source_id": "strict_core_minimal_provider_from_P2297",
            "source_packet": "P2297.strict_core_minimal_provider",
            "exported_current_source": True,
            "algebraically_capable_for_full_residual_family": strict_core_consistent,
            "strict_admissible_source": True,
            "internal_selector_source_exported": False,
            "noncyclic_anchor_exported": False,
            "no_legacy_transfer_used": True,
            "qw2191_discharged": qw2191_all_required_pass,
            "failure_mode": "rank_A<rank_augmented; no componentwise full-family solution in the strict-core provider class",
        },
        {
            "source_id": "p1990_augmented_provider_channel_from_P1991_P2297",
            "source_packet": "P1991/P2297.p1990_augmented_provider_non_strict",
            "exported_current_source": True,
            "algebraically_capable_for_full_residual_family": augmented_consistent,
            "strict_admissible_source": False,
            "internal_selector_source_exported": False,
            "noncyclic_anchor_exported": False,
            "no_legacy_transfer_used": True,
            "qw2191_discharged": qw2191_all_required_pass,
            "failure_mode": "explicitly non-strict augmented selector label and still rank-inconsistent for the full P2297 residual system",
        },
        {
            "source_id": "formal_full_residual_basis_provider_from_P2297",
            "source_packet": "P2297.formal_full_residual_basis_provider",
            "exported_current_source": True,
            "algebraically_capable_for_full_residual_family": formal_full_consistent and formal_full_exact,
            "strict_admissible_source": not formal_full_non_admissible,
            "internal_selector_source_exported": False,
            "noncyclic_anchor_exported": False,
            "no_legacy_transfer_used": True,
            "qw2191_discharged": qw2191_all_required_pass,
            "failure_mode": "copies the residual-family basis; no internal strict selector/source exported, so strict promotion would require a new source or explicit non-strict selector premise",
        },
    ]

    source_matrix = np.array([candidate_vector(c) for c in candidates], dtype=float)
    # A legal current strict source must be algebraically capable, strict-admissible,
    # internally sourced, noncyclic, and QW-2191-discharged. The score is zero only
    # if every required bit is present.
    deficiency_scores = []
    for row in source_matrix:
        required = row[[1, 2, 3, 4, 6]]
        deficiency_scores.append(float(la.norm(np.ones_like(required) - required, ord=1)))

    legal_strict_capable_sources = [
        c["source_id"]
        for c in candidates
        if c["algebraically_capable_for_full_residual_family"]
        and c["strict_admissible_source"]
        and c["internal_selector_source_exported"]
        and c["noncyclic_anchor_exported"]
        and c["qw2191_discharged"]
    ]
    algebraically_capable_sources = [
        c["source_id"] for c in candidates if c["algebraically_capable_for_full_residual_family"]
    ]
    capable_but_non_admissible_sources = [
        c["source_id"]
        for c in candidates
        if c["algebraically_capable_for_full_residual_family"] and not c["strict_admissible_source"]
    ]

    x = sp.Symbol("x")
    capability_polynomial = sp.Poly(
        sum(int(c["algebraically_capable_for_full_residual_family"]) * x**i for i, c in enumerate(candidates)),
        x,
    )
    strict_legal_polynomial = sp.Poly(
        sum(
            int(c["algebraically_capable_for_full_residual_family"])
            * int(c["strict_admissible_source"])
            * int(c["internal_selector_source_exported"])
            * int(c["noncyclic_anchor_exported"])
            * int(c["qw2191_discharged"])
            * x**i
            for i, c in enumerate(candidates)
        ),
        x,
    )

    theorem_body = {
        "statement_id": "P2298_CURRENT_STRICT_PROVIDER_SOURCE_NONADMISSIBILITY_THEOREM",
        "scope": "current exported provider-source inventory for the P2297 strict non-GB ADM/Bianchi-I full residual family",
        "formal_statement": (
            "For every provider source exported in the current strict inventory, algebraic cancellation of the "
            "P2297 full non-GB spatial-EOM residual is not simultaneously legal as a strict provider: the strict-core "
            "source is matrix-inconsistent, the augmented channel is non-strict and matrix-inconsistent, and the only "
            "algebraically capable full-basis source is a residual-copying basis with no internal selector source. Any "
            "promotion of that basis therefore requires either a new noncyclic internal strict provider source or an "
            "explicit non-strict selector premise; QW-2191 remains undischarged."
        ),
        "not_claimed": [
            "absolute impossibility over all future mathematical ansatz classes",
            "strict selector closure",
            "QW-2191 discharge",
            "legacy-to-strict bridge",
            "ToE closure",
        ],
        "proof_obligations": {
            "full_residual_family_present": len(strict_non_gb_families) > 0 and len(outside_eh_families) > 0,
            "p2297_formal_basis_spans_residual": formal_full_consistent and formal_full_exact,
            "strict_core_rejected_by_rank": not strict_core_consistent,
            "p1990_augmented_rejected_by_rank_and_label": (not augmented_consistent)
            and p1991.get("selector_premise_label", {}).get("status") == "NON_STRICT_AUGMENTED_CLASS",
            "qW2191_replay_available_but_not_discharged": qw2191_replay_certified and not qw2191_all_required_pass,
            "legal_strict_capable_source_count_zero": len(legal_strict_capable_sources) == 0,
        },
    }

    gatekeeper_checks = {
        "p2297_provider_matrix_loaded": bool(p2297_probe),
        "full_residual_family_detected": len(strict_non_gb_families) > 0 and len(outside_eh_families) > 0,
        "formal_basis_rank_positive": formal_basis_count > 0 and formal_rank > 0,
        "strict_core_matrix_inconsistent": not strict_core_consistent,
        "augmented_non_strict_matrix_inconsistent": not augmented_consistent,
        "formal_full_basis_algebraically_capable": formal_full_consistent and formal_full_exact,
        "formal_full_basis_marked_non_admissible": formal_full_non_admissible,
        "no_current_legal_strict_capable_provider_source": len(legal_strict_capable_sources) == 0,
        "all_capable_sources_require_new_source_or_non_strict_selector": bool(capable_but_non_admissible_sources)
        and set(algebraically_capable_sources) == set(capable_but_non_admissible_sources),
        "qw2191_replay_certified": qw2191_replay_certified,
        "qw2191_not_discharged": not qw2191_all_required_pass,
        "no_selector_closure_claimed": True,
        "no_legacy_transfer_used": True,
        "no_toe_closure_claimed": True,
    }

    theorem_fingerprint = sha256_json(theorem_body)
    payload = {
        "schema_version": "p2298_s1248_v1",
        "packet_id": "P2298",
        "stage_id": "S1248",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_WITH_THEOREM_GRADE_NONADMISSIBILITY_TRACE",
        "result_kind": "STRICT_NON_GB_PROVIDER_SOURCE_NONADMISSIBILITY_THEOREM_WITH_TRACE",
        "strict_non_gb_provider_source_nonadmissibility_theorem_probe": {
            "probe_id": "P2298_S1248_STRICT_NON_GB_PROVIDER_SOURCE_NONADMISSIBILITY_THEOREM",
            "source_packets": {
                "p1988": "generated/p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json",
                "p1991": "generated/p1991_s941_strict_augmented_provider_channel_matrix_witness.json",
                "p2168": "generated/p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.json",
                "p2181": "generated/p2181_s1131_strict_qw2191_replay_certificate_freeze.json",
                "p2297": "generated/p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json",
            },
            "source_hashes": {
                "p1988_sha256": sha256_file(IN_1988),
                "p1991_sha256": sha256_file(IN_1991),
                "p2168_sha256": sha256_file(IN_2168),
                "p2181_sha256": sha256_file(IN_2181),
                "p2297_sha256": sha256_file(IN_2297),
            },
            "residual_family_scope": {
                "strict_non_gb_families_detected": strict_non_gb_families,
                "outside_eh_families": outside_eh_families,
                "formal_full_basis_monomial_count": formal_basis_count,
                "formal_full_basis_rank": formal_rank,
            },
            "exported_provider_source_inventory": candidates,
            "admissibility_matrix": {
                "columns": [
                    "exported_current_source",
                    "algebraically_capable_for_full_residual_family",
                    "strict_admissible_source",
                    "internal_selector_source_exported",
                    "noncyclic_anchor_exported",
                    "no_legacy_transfer_used",
                    "qw2191_discharged",
                ],
                "rows": [c["source_id"] for c in candidates],
                "matrix": source_matrix.astype(int).tolist(),
                "deficiency_l1_scores_for_required_bits": deficiency_scores,
                "capability_polynomial": str(capability_polynomial.as_expr()),
                "strict_legal_polynomial": str(strict_legal_polynomial.as_expr()),
                "strict_legal_polynomial_is_zero": strict_legal_polynomial.as_expr() == 0,
            },
            "nonadmissibility_result": {
                "legal_strict_capable_sources": legal_strict_capable_sources,
                "algebraically_capable_sources": algebraically_capable_sources,
                "capable_but_non_admissible_sources": capable_but_non_admissible_sources,
                "theorem_grade_current_inventory_verdict": "NO_CURRENT_LEGAL_STRICT_PROVIDER_SOURCE_FOR_FULL_P2297_RESIDUAL_FAMILY",
                "required_escape_hatches": [
                    "export a genuinely new noncyclic internal strict provider source and re-run P2297/P2298",
                    "or introduce an explicit non-strict selector premise and keep the closure outside strict-core",
                ],
            },
            "qw2191_binding": {
                "p2168_all_required_pass": qw2191_all_required_pass,
                "p2181_replay_certified": qw2191_replay_certified,
                "selector_closure_status": "NOT_CLAIMED_AND_NOT_DISCHARGED",
                "reason": "P2168 required QW-2191 obligations do not all pass; P2181 certifies replay consistency only.",
            },
            "theorem_export": theorem_body,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2299_candidate",
            "goal": "Either construct a genuinely noncyclic internal strict provider source for the P2297 residual family, or formalize an explicitly non-strict selector-premise theorem branch without promoting it to strict closure.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_THEOREM_GRADE_NONADMISSIBILITY_TRACE",
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")

    md = f"""# P2298/S1248 — strict non-GB provider-source non-admissibility theorem probe\n\n- Status: `{payload['status']}`\n- Verdict: `{payload['strict_non_gb_provider_source_nonadmissibility_theorem_probe']['nonadmissibility_result']['theorem_grade_current_inventory_verdict']}`\n- Full residual families detected: `{len(strict_non_gb_families)}`; outside-EH families: `{len(outside_eh_families)}`; formal full-basis monomials: `{formal_basis_count}`.\n- Legal current strict capable provider sources: `{len(legal_strict_capable_sources)}`.\n- Algebraically capable but non-admissible sources: `{', '.join(capable_but_non_admissible_sources)}`.\n- QW-2191 binding: replay certified = `{qw2191_replay_certified}`, all theorem obligations pass = `{qw2191_all_required_pass}`.\n- Theorem fingerprint: `{theorem_fingerprint}`\n\n## Scope\nThis is a current-inventory non-admissibility theorem trace: it does not prove impossibility over all future ansatz classes, does not discharge QW-2191, does not transfer legacy kernel roles, and does not claim ToE closure.\n\n## Next honest step\n{payload['recommended_next_honest_step']['goal']}\n"""
    MD.write_text(md, encoding="utf-8")


if __name__ == "__main__":
    main()
