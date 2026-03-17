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

IN_P709 = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json"
IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe.json"
IN_P716 = GENERATED / "p716_current_strict_t176_pair4_negative_cosine_polarity_global_z2_orbit_split_audit_probe_summary.json"

OUT_JSON = GENERATED / "p717_current_strict_t176_pair4_exact_branch_split_release_7_os_gauge_irrelevance_bridge_audit_probe.json"
OUT_SUMMARY = GENERATED / "p717_current_strict_t176_pair4_exact_branch_split_release_7_os_gauge_irrelevance_bridge_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def negate_sign_vector(sign_vector: dict[str, int]) -> dict[str, int]:
    return {pair: int(-value) for pair, value in sign_vector.items()}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P709, IN_P715, IN_P716]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P717",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p709 = load_json(IN_P709)
    p715 = load_json(IN_P715)
    p716 = load_json(IN_P716)

    result_715 = p715.get("result") or {}
    rooted_results = p715.get("rooted_results") or {}
    reference_root = result_715.get("reference_root")
    reference_sign_vector = (rooted_results.get(reference_root) or {}).get("sign_vector_by_pair") if isinstance(reference_root, str) else None
    pair4_sign_vector = (rooted_results.get("pair4") or {}).get("sign_vector_by_pair")

    pair4_is_global_negation = False
    if isinstance(reference_sign_vector, dict) and isinstance(pair4_sign_vector, dict):
        pair4_is_global_negation = pair4_sign_vector == negate_sign_vector(reference_sign_vector)

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

    add_check(
        "p715_exact_branch_split_still_present",
        bool(result_715.get("all_roots_supported"))
        and bool(result_715.get("projective_root_independent_sign_orbit"))
        and not bool(result_715.get("exact_root_independent_sign_vector")),
        True,
        "P715 still exhibits the all-root projective orbit together with failure of one exact directed branch.",
    )
    add_check(
        "p715_pair4_branch_is_global_negation_of_reference",
        pair4_is_global_negation,
        True,
        "The exact pair4 branch is still the global negation of the reference branch on current exports.",
    )
    add_check(
        "p716_localizes_split_to_pair4_negative_cosine_polarity",
        bool(p716.get("current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity")),
        True,
        "P716 still localizes the exact split to the pair4 negative cosine-axis role.",
    )
    add_check(
        "p709_release_7_os_sign_gauge_irrelevance_still_passes",
        p709.get("status"),
        "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED",
        "P709 still proves sign-gauge irrelevance for the concrete Release-7 OS observables.",
    )
    add_check(
        "p709_full_sign_family_audit_still_ok",
        bool(p709.get("sign_ok")),
        True,
        "P709 still passes over the full residual sign-pattern family, so in particular it covers the present exact two-branch split.",
    )

    status = (
        "PASS_PAIR4_EXACT_BRANCH_SPLIT_RELEASE_7_OS_GAUGE_IRRELEVANCE_BRIDGED"
        if not blocking
        else "P717_REQUIRES_REVIEW_CHANGED_PAIR4_BRANCH_OR_OS_GAUGE_STATE"
    )

    artifact = {
        "stage": "P717",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit whether the exact pair4-induced branch split isolated by P715/P716 is already gauge-irrelevant "
            "for the concrete Release-7 OS observables covered by P709"
        ),
        "inputs": {
            "P709_summary": str(IN_P709.relative_to(REPO)),
            "P715": str(IN_P715.relative_to(REPO)),
            "P716_summary": str(IN_P716.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "result": {
            "pair4_exact_branch_split_present": bool(
                bool(result_715.get("all_roots_supported"))
                and bool(result_715.get("projective_root_independent_sign_orbit"))
                and not bool(result_715.get("exact_root_independent_sign_vector"))
            ),
            "pair4_exact_branch_is_global_negation_of_reference": pair4_is_global_negation,
            "pair4_split_localized_to_negative_cosine_polarity": bool(
                p716.get("current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity")
            ),
            "release_7_os_sign_gauge_irrelevance_available": p709.get("status")
            == "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED",
            "pair4_exact_branch_split_gauge_irrelevant_for_release_7_os_observables": len(blocking) == 0,
            "covered_release_7_observables": [
                "P694_mass_proxy",
                "P696_selector_aligned_channel_spectrum_proxy_invariants",
                "F704_basis_invariant_mass_observable",
            ],
            "counts_as_strict_physical_orientation_datum": False,
            "implies_t176_discharge": False,
            "remaining_gap_after_positive_result": (
                "exact branch-fixing in strict core remains open if one wants more than projective/gauge-safe Release-7 OS scope"
            ),
        },
        "hard_limits": [
            "No strict-core directed/sign-sensitive physical orientation datum claim.",
            "No T176 discharge claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No gauge-irrelevance claim for arbitrary future observables.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P717",
        "status": status,
        "as_of": AS_OF,
        "pair4_exact_branch_split_gauge_irrelevant_for_release_7_os_observables": len(blocking) == 0,
        "counts_as_strict_physical_orientation_datum": False,
        "implies_t176_discharge": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
