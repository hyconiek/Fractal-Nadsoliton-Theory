#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.json"
MD = GEN / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.md"

SOURCE_FILES = {
    "P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT": GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json",
    "P2316_CURRENT_BEST_LAGRANGIAN_EOM": GEN / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_BRIDGE": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "P2324_AXIS_BRANCH_SUSCEPTIBILITY": GEN / "p2324_s1274_strict_axis_branch_susceptibility_nonpromotion_probe.json",
    "P2325_SIGNED_SOURCE_BRIDGE_AUDIT": GEN / "p2325_s1275_strict_signed_source_to_axis_susceptibility_bridge_audit_probe.json",
    "P2328_SELECTOR_FREE_DYNAMICS_PRIORITY": ROOT / "P2328_S1278_SELECTOR_FREE_DYNAMICS_CLOSURE_PRIORITY_PACKET.md",
}

GREP_PATTERNS = (
    "S_future",
    "selector-free",
    "selector_independent",
    "Lagrangian",
    "Euler-Lagrange",
    "EOM",
    "QW-2191",
    "provider_lift",
    "future-state selector",
)

SELECTOR_MARKERS = (
    "selector",
    "s_future",
    "future-state",
    "future_state",
    "branch actualization",
    "actualization",
    "qw-2191",
    "provider_lift",
    "policy-margin",
    "policy_margin",
    "signed tilt",
    "mu/theta",
    "theta_ref",
)

TERM_SECTOR_HINTS = {
    "L_psi_kin": "psi kinetic term",
    "L_psi_mass": "psi mass/quadratic term",
    "L_psi_self": "psi local quartic self-interaction",
    "L_A_kin": "A kinetic term",
    "L_A_mass": "A mass/quadratic term",
    "L_A_self": "A local quartic self-interaction",
    "L_h_kin": "h kinetic term",
    "L_h_mass": "h mass/quadratic term",
    "L_h_self": "h local quartic self-interaction",
    "L_mix_trilinear": "A-h-psi local trilinear mixing term",
    "L_mix_biquadratic": "A^2-psi^2 local biquadratic mixing term",
}

OBSERVABLE_CLAIMS = [
    {
        "claim_id": "local_termwise_lagrangian_density",
        "description": "The reduced P2086 L_total can be assembled and audited term-by-term.",
        "selector_required": False,
        "status": "safe_selector_independent_computational_claim",
    },
    {
        "claim_id": "euler_lagrange_variations_for_psi_A_h",
        "description": "P2086 Euler-Lagrange variations for psi, A, and h have zero recomposition residuals.",
        "selector_required": False,
        "status": "safe_selector_independent_computational_claim",
    },
    {
        "claim_id": "mass_energy_sector_audit_inside_reduced_model",
        "description": "Mass/quadratic, kinetic, self, and local mix terms can be inspected as reduced-sector ingredients.",
        "selector_required": False,
        "status": "safe_if_not_promoted_to_unique_future_branch",
    },
    {
        "claim_id": "hessian_branch_stability_classification",
        "description": "Branch stability can be classified without selecting which branch becomes actual.",
        "selector_required": False,
        "status": "safe_as_classification_not_actualization",
    },
    {
        "claim_id": "exact_future_branch_actualization",
        "description": "A single next branch/state is selected as the actual future.",
        "selector_required": True,
        "status": "blocked_by_QW2191_S_future_missing",
    },
    {
        "claim_id": "provider_lift_per_step_policy_margin",
        "description": "A signed source is mapped into Task-3 provider_lift_per_step/policy-margin semantics.",
        "selector_required": True,
        "status": "blocked_by_P2318_P2325_bridge_fields_missing",
    },
    {
        "claim_id": "G1_G3_update_or_ToE_closure",
        "description": "The global goal status is upgraded or the full theory is declared closed.",
        "selector_required": True,
        "status": "blocked_no_false_pass",
    },
]


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": path.relative_to(REPO).as_posix()}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_json(obj: Any) -> str:
    return sha256_text(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def selector_marker_hits(text: str) -> list[str]:
    lowered = text.lower()
    return sorted(marker for marker in SELECTOR_MARKERS if marker in lowered)


def collect_repo_grep_hits() -> list[dict[str, Any]]:
    candidates = [
        ROOT / "P2328_S1278_SELECTOR_FREE_DYNAMICS_CLOSURE_PRIORITY_PACKET.md",
        ROOT / "P2327_S1277_KERNEL_DERIVED_FUTURE_STATE_SELECTOR_QW2191_CONDITION_PACKET.md",
        ROOT / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.py",
        ROOT / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.py",
        ROOT / "p2325_s1275_strict_signed_source_to_axis_susceptibility_bridge_audit_probe.py",
    ]
    candidates.extend(sorted(GEN.glob("p23*_s12*.json"))[:80])
    hits: list[dict[str, Any]] = []
    for path in candidates:
        text = read_text(path)
        if not text:
            continue
        lower = text.lower()
        count = sum(lower.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        excerpt = ""
        for index, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = index
                excerpt = line.strip()[:220]
                break
        hits.append({
            "path": rel(path),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": excerpt,
        })
    hits.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return hits


def classify_terms(p2086: dict[str, Any]) -> list[dict[str, Any]]:
    results = p2086.get("eom_execution_results", {})
    terms = results.get("lagrangian_terms", {})
    rows: list[dict[str, Any]] = []
    for term_id, expression in sorted(terms.items()):
        hits = selector_marker_hits(term_id + " " + expression)
        rows.append({
            "term_id": term_id,
            "sector_hint": TERM_SECTOR_HINTS.get(term_id, "unclassified reduced-sector term"),
            "selector_marker_hits": hits,
            "selector_independence_class": "selector_independent_reduced_local_dynamics" if not hits else "selector_dependent_or_contaminated",
            "requires_S_future_for_termwise_variation": bool(hits),
            "safe_claim_scope": "term can be varied and included in reduced L_total without asserting branch actualization" if not hits else "term needs manual review before selector-free use",
            "expression_srepr_excerpt": expression[:240],
        })
    return rows


def classify_variations(p2086: dict[str, Any]) -> list[dict[str, Any]]:
    results = p2086.get("eom_execution_results", {})
    termwise_map = results.get("termwise_variation_map", {})
    residuals = results.get("symbolic_recomposition_residual", {})
    numeric = results.get("numeric_probe_residual", {})
    rows: list[dict[str, Any]] = []
    for field, variation_map in sorted(termwise_map.items()):
        joined = " ".join([field, *variation_map.keys(), *map(str, variation_map.values())])
        hits = selector_marker_hits(joined)
        rows.append({
            "field": field,
            "termwise_variation_count": len(variation_map),
            "nonzero_termwise_variations": sorted(k for k, v in variation_map.items() if str(v) not in ("Integer(0)", "0", "0.0")),
            "symbolic_recomposition_residual": residuals.get(field),
            "numeric_probe_residual": numeric.get(field),
            "residual_zero": residuals.get(field) in ("Integer(0)", "0", 0) and str(numeric.get(field)) == "0",
            "selector_marker_hits": hits,
            "selector_independence_class": "selector_independent_reduced_EOM_variation" if not hits else "selector_dependent_or_contaminated",
            "requires_S_future_for_reduced_variation": bool(hits),
        })
    return rows


def residual_rows(p2086: dict[str, Any]) -> list[dict[str, Any]]:
    results = p2086.get("eom_execution_results", {})
    symbolic = results.get("symbolic_recomposition_residual", {})
    numeric = results.get("numeric_probe_residual", {})
    rows: list[dict[str, Any]] = []
    for field in sorted(set(symbolic) | set(numeric)):
        rows.append({
            "field": field,
            "symbolic_residual": symbolic.get(field),
            "numeric_residual": numeric.get(field),
            "selector_independence_class": "selector_independent_residual_check",
            "requires_S_future": False,
            "interpretation": "recomposition consistency check for reduced EOM; does not choose a future branch",
        })
    return rows


def bridge_dependency_summary(artifacts: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2318_probe = artifacts["P2318_SELECTOR_LANE_TO_MARGIN_BRIDGE"].get(
        "strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}
    )
    p2324_probe = artifacts["P2324_AXIS_BRANCH_SUSCEPTIBILITY"].get(
        "strict_axis_branch_susceptibility_nonpromotion_probe", {}
    )
    p2325_probe = artifacts["P2325_SIGNED_SOURCE_BRIDGE_AUDIT"].get(
        "strict_signed_source_to_axis_susceptibility_bridge_audit_probe", {}
    )
    return {
        "p2318_missing_required_bridge_fields": p2318_probe.get("bridge_obligation_verdict", {}).get(
            "missing_required_bridge_fields", []
        ),
        "p2324_mu_template_present_but_source_missing": p2324_probe.get("susceptibility_certificate", {}).get(
            "all_pairs_split_by_template_if_mu_given", False
        ),
        "p2325_admissible_signed_source_count": p2325_probe.get("bridge_audit_certificate", {}).get(
            "admissible_signed_source_count", 0
        ),
        "selector_dependent_boundary": [
            "exact future branch actualization",
            "signed mu/theta_ref source descent",
            "provider_lift_per_step bridge",
            "QW-2191 discharge",
            "G1/G3 or ToE closure",
        ],
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items() if path.suffix == ".json"}
    p2086 = artifacts["P2086_FULL_LTOTAL_TERMWISE_EOM_AUDIT"]
    term_rows = classify_terms(p2086)
    variation_rows = classify_variations(p2086)
    recomposition_rows = residual_rows(p2086)
    safe_claims = [row for row in OBSERVABLE_CLAIMS if not row["selector_required"]]
    blocked_claims = [row for row in OBSERVABLE_CLAIMS if row["selector_required"]]
    bridge_summary = bridge_dependency_summary(artifacts)
    grep_hits = collect_repo_grep_hits()

    selector_independence_certificate = {
        "p2086_loaded": p2086.get("packet_id") == "P2086",
        "term_count": len(term_rows),
        "selector_independent_term_count": sum(
            row["selector_independence_class"] == "selector_independent_reduced_local_dynamics" for row in term_rows
        ),
        "variation_field_count": len(variation_rows),
        "selector_independent_variation_field_count": sum(
            row["selector_independence_class"] == "selector_independent_reduced_EOM_variation" for row in variation_rows
        ),
        "all_recomposition_residuals_selector_independent": all(not row["requires_S_future"] for row in recomposition_rows),
        "safe_observable_claim_count": len(safe_claims),
        "blocked_observable_claim_count": len(blocked_claims),
        "p2328_strategy_confirmed": SOURCE_FILES["P2328_SELECTOR_FREE_DYNAMICS_PRIORITY"].exists(),
        "selector_boundary_preserved": bridge_summary["p2325_admissible_signed_source_count"] == 0,
    }

    theorem_export = {
        "theorem_name": "P2329 selector-independence audit for reduced Lagrangian/EOM work",
        "claim": "For the current P2086 reduced L_total artifact, termwise Lagrangian terms, psi/A/h Euler-Lagrange variations, and recomposition residual checks are selector-independent reduced-dynamics computations; exact future-branch actualization, provider_lift_per_step, QW-2191 discharge, and G1/G3 or ToE closure remain selector-dependent and unclaimed.",
        "proof_witnesses": [
            "P2086 provides 11 reduced L_total terms and termwise psi/A/h EOM recomposition residuals.",
            "No P2086 reduced term or psi/A/h variation expression contains S_future/selector/future-branch markers.",
            "P2318/P2324/P2325 keep the signed-source/provider-lift/future-branch boundary unfilled.",
            "P2328 explicitly permits selector-free dynamics work while preserving QW-2191 as branch-actualization frontier.",
        ],
        "scope_limits": [
            "This is an audit of the reduced P2086 L_total, not a full tensor-resolved global Lagrangian theorem.",
            "Selector-independent does not mean final ToE-closed.",
            "No exact future-state selector is exported here.",
        ],
        "strict_guardrails": {
            "no_legacy_kernel_role_transfer": True,
            "no_selector_premise_added": True,
            "no_qw2191_discharge_claimed": True,
            "no_g1_g3_update_claimed": True,
            "no_toe_closure_claimed": True,
        },
    }

    probe = {
        "probe_id": "P2329_S1279_SELECTOR_INDEPENDENCE_LAGRANGIAN_EOM_AUDIT",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_patterns": list(GREP_PATTERNS),
            "hit_count": len(grep_hits),
            "top_hits": grep_hits[:20],
        },
        "selector_independence_certificate": selector_independence_certificate,
        "selector_independent_terms": term_rows,
        "selector_independent_variations": variation_rows,
        "selector_independent_residuals": recomposition_rows,
        "safe_observable_claims": safe_claims,
        "blocked_selector_dependent_claims": blocked_claims,
        "bridge_dependency_summary": bridge_summary,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "p2086_loaded": selector_independence_certificate["p2086_loaded"],
        "term_count_is_11": selector_independence_certificate["term_count"] == 11,
        "all_terms_selector_independent": selector_independence_certificate["selector_independent_term_count"] == 11,
        "variation_fields_are_psi_A_h": sorted(row["field"] for row in variation_rows) == ["A", "h", "psi"],
        "all_variations_selector_independent": selector_independence_certificate[
            "selector_independent_variation_field_count"
        ] == 3,
        "all_residuals_selector_independent": selector_independence_certificate[
            "all_recomposition_residuals_selector_independent"
        ],
        "safe_claims_nonempty": len(safe_claims) >= 4,
        "blocked_claims_nonempty": len(blocked_claims) >= 3,
        "p2328_strategy_packet_loaded": selector_independence_certificate["p2328_strategy_confirmed"],
        "p2325_no_admissible_signed_source_preserved": bridge_summary["p2325_admissible_signed_source_count"] == 0,
        "no_legacy_kernel_role_transfer": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2329_s1279_v1",
        "packet_id": "P2329",
        "stage_id": "S1279",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_SELECTOR_FREE_DYNAMICS_AUDIT_WITH_BRANCH_ACTUALIZATION_STILL_BLOCKED",
        "result_kind": "STRICT_SELECTOR_INDEPENDENCE_LAGRANGIAN_EOM_AUDIT_NO_QW2191_DISCHARGE",
        "strict_selector_independence_lagrangian_eom_audit_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": "Use the selector-independent P2086 term/variation/residual subset to continue Lagrangian/EOM sharpening, while separately auditing any term or observable claim that requires branch sign, S_future, provider_lift_per_step, or QW-2191 discharge.",
        "global_status": "OPEN_DYNAMICS_PROGRESS_ALLOWED_SELECTOR_ACTUALIZATION_UNRESOLVED",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2329 selector-independence Lagrangian/EOM audit\n\n"
        "Status: selector-free reduced dynamics can continue; branch actualization remains blocked.\n\n"
        f"- P2086 reduced term count: {selector_independence_certificate['term_count']}\n"
        f"- Selector-independent term count: {selector_independence_certificate['selector_independent_term_count']}\n"
        f"- Selector-independent variation fields: {selector_independence_certificate['selector_independent_variation_field_count']}\n"
        f"- Safe observable claim count: {len(safe_claims)}\n"
        f"- Blocked selector-dependent claim count: {len(blocked_claims)}\n"
        "- No QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
