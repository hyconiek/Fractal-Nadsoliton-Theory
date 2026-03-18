#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P723 = GENERATED / "p723_current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_audit_probe_summary.json"
IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"

OUT_JSON = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe.json"
OUT_SUMMARY = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P723, IN_P715, IN_F141, IN_F143, IN_F147]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P724",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p723 = load_json(IN_P723)
    p715 = load_json(IN_P715)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    f147 = load_json(IN_F147)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    p715_result = p715.get("result") or {}
    rooted_results = p715.get("rooted_results") or {}

    source_positive_polarity_available = (
        ((f141.get("support_packet") or {}).get("psi_src_barrier_sign_component_witness_v1") == 1)
        and (float(f143.get("scalar_component_witness_value") or 0.0) > 0.0)
        and bool((((f147.get("support_packet") or {}).get("selector_signed_split_realization") or {}).get("positive_signed_selector_response")))
        and bool((((f147.get("support_packet") or {}).get("preobserver_scope_realization") or {}).get("positive_plus_output")))
        and bool((((f147.get("support_packet") or {}).get("preobserver_scope_realization") or {}).get("vanishing_minus_output")))
    )

    positive_roots: list[str] = []
    negative_roots: list[str] = []
    for pair_name, rooted in sorted(rooted_results.items()):
        anchor_scalar = rooted.get("root_anchor_scalar")
        if not isinstance(anchor_scalar, (int, float)):
            continue
        if float(anchor_scalar) > 0.0:
            positive_roots.append(pair_name)
        elif float(anchor_scalar) < 0.0:
            negative_roots.append(pair_name)

    anchor_source_by_root = p715.get("anchor_source_by_root") or p715_result.get("anchor_source_by_root") or {}
    surviving_positive_odd_lane = [pair for pair in positive_roots if anchor_source_by_root.get(pair) == "w_break"]
    surviving_positive_even_fallback_lane = [
        pair for pair in positive_roots if anchor_source_by_root.get(pair) == "w_ref_unnormalized"
    ]

    add_check(
        "p723_chart_seed_bridge_still_not_exported",
        bool(p723.get("t178_target_exported_on_current_repo_state")),
        False,
        "P723 already confirms that the source-topology-to-atlas chart-seed bridge is not yet exported.",
    )
    add_check(
        "source_positive_polarity_available",
        source_positive_polarity_available,
        True,
        "The current source-topology lane exports a positive barrier sign, positive flow, and positive selector plus-channel polarity.",
    )
    add_check(
        "p715_all_roots_supported",
        bool(p715_result.get("all_roots_supported")),
        True,
        "The currently strongest atlas-entry candidate family does support all five roots.",
    )
    add_check(
        "positive_source_polarity_compatible_roots",
        positive_roots,
        ["pair1", "pair2", "pair3", "pair5"],
        "Within the current all-root transported-family candidate, positive source polarity is compatible exactly with the four positive-root entries.",
    )
    add_check(
        "positive_source_polarity_incompatible_roots",
        negative_roots,
        ["pair4"],
        "The same positive source polarity excludes the unique negative atlas-entry branch pair4 on current exports.",
    )
    add_check(
        "positive_polarity_compatible_roots_match_same_projective_orbit_roots",
        positive_roots,
        list(p715_result.get("same_orbit_roots_relative_to_reference") or []),
        "The positive-polarity-compatible roots are exactly the same-orbit roots in the current projective atlas family.",
    )
    add_check(
        "unique_chart_seed_selected",
        len(positive_roots) == 1,
        False,
        "Positive source polarity alone still does not select one unique chart seed.",
    )
    add_check(
        "residual_positive_corridor_splits_into_odd_and_even_fallback_lanes",
        {
            "odd_lane": surviving_positive_odd_lane,
            "even_fallback_lane": surviving_positive_even_fallback_lane,
        },
        {
            "odd_lane": ["pair1", "pair5"],
            "even_fallback_lane": ["pair2", "pair3"],
        },
        "The surviving positive corridor still splits into an odd-anchor lane and an even-fallback lane on current exports.",
    )

    status = (
        "PARTIAL_POSITIVE_SOURCE_POLARITY_REDUCES_ATLAS_ENTRY_CORRIDOR_ONLY"
        if not blocking
        else "P724_REQUIRES_REVIEW_CHANGED_POSITIVE_POLARITY_CORRIDOR_STATE"
    )

    artifact = {
        "stage": "P724",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_only",
        "inputs": {
            "P723": str(IN_P723.relative_to(REPO)),
            "P715": str(IN_P715.relative_to(REPO)),
            "F141": str(IN_F141.relative_to(REPO)),
            "F143": str(IN_F143.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "source_positive_polarity_available": source_positive_polarity_available,
        "atlas_entry_roots_compatible_with_current_positive_source_polarity": positive_roots,
        "atlas_entry_roots_incompatible_with_current_positive_source_polarity": negative_roots,
        "unique_chart_seed_selected": False,
        "residual_positive_corridor_split": {
            "odd_anchor_lane": surviving_positive_odd_lane,
            "even_fallback_lane": surviving_positive_even_fallback_lane,
        },
        "audit_conclusion": {
            "current_positive_source_polarity_rules_out_pair4": negative_roots == ["pair4"],
            "current_positive_source_polarity_reduces_but_does_not_close_t178_gap": True,
            "next_honest_move": (
                "export_or_attack_a_finer_source_to_atlas_seed_rule_distinguishing_within_the_positive_entry_corridor"
            ),
        },
        "hard_limits": [
            "No T178 discharge claim.",
            "No T177 discharge claim.",
            "No unique atlas chart-seed selection claim.",
            "No strict physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P724",
        "status": status,
        "as_of": AS_OF,
        "source_positive_polarity_available": source_positive_polarity_available,
        "atlas_entry_roots_compatible_with_current_positive_source_polarity": positive_roots,
        "atlas_entry_roots_incompatible_with_current_positive_source_polarity": negative_roots,
        "unique_chart_seed_selected": False,
        "next_honest_move": "finer_source_to_atlas_seed_rule_within_positive_corridor",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
