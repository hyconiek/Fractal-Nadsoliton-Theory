#!/usr/bin/env python3
"""P3037/S1987: selector-mechanism hint sheaf acceptance matrix.

Assume the right selector mechanism may not fit a familiar human schema, while
existing research still contains hints.  This step does not replay a concrete
selector candidate.  It constructs a finite typed hint sheaf: each hint source is
mapped to feature atoms that an unknown selector mechanism would have to glue.
The computation checks whether current hints already form one accepted glued
mechanism.  They do not; the object is a proof obligation map for the next move.
"""
from __future__ import annotations

import hashlib, json, re
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3035_s1985_z12_directional_branch_selector_source_obstruction import OUT as P3035
from p3036_s1986_external_physical_unit_source_scale_orbit_obstruction import OUT as P3036

OUT = GEN / "p3037_s1987_selector_mechanism_hint_sheaf_acceptance_matrix.json"
MD = GEN / "p3037_s1987_selector_mechanism_hint_sheaf_acceptance_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
V1 = ROOT / "V1_INFORMATIONAL_VISCOSITY_HYPOTHESIS_AUDIT.md"
H42 = ROOT / "H42_C_BASED_RETARDATION_OPERATOR_ON_PAIR1_AUDIT.md"

FEATURES = [
    "damping_memory_feedback_hint",
    "retardation_anisotropy_anchor_hint",
    "inversion_odd_or_chiral_sign_hint",
    "nonpremise_source_localizer_needed",
    "coupling_polarity_needed",
    "readout_unit_coupling_needed",
]

HINT_SOURCES = [
    {
        "source": "V1 informational viscosity hypothesis",
        "path": V1,
        "patterns": {
            "damping_memory_feedback_hint": r"viscosity|damping|memory|feedback|lepko",
            "nonpremise_source_localizer_needed": r"no explicit operator|no current file exports|not an existing solved lane|not.*already_supported_selector_mechanism",
        },
    },
    {
        "source": "H42 c-based retardation operator",
        "path": H42,
        "patterns": {
            "retardation_anisotropy_anchor_hint": r"retardation|anisotropic|psi0|spectral/response split|Delta_ret",
            "nonpremise_source_localizer_needed": r"without an imported anchor|does not by itself generate selector orientation|no claim.*QW-2191",
        },
    },
    {
        "source": "P3035 directional branch-selector obstruction",
        "path": P3035,
        "patterns": {
            "inversion_odd_or_chiral_sign_hint": r"inversion|chiral|orientation torsor|\+direction|-direction",
            "coupling_polarity_needed": r"readout coupling|coupling|fixed sign|branch source",
        },
    },
    {
        "source": "P3036 external unit-source obstruction",
        "path": P3036,
        "patterns": {
            "readout_unit_coupling_needed": r"physical length/action/clock unit|unit-bearing action|readout|unit theorem|scale torsor",
            "nonpremise_source_localizer_needed": r"conditional gauges|no.*unit export|dimensionless representative",
        },
    },
]


def text_for(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def hit(pattern: str, text: str) -> bool:
    return re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE) is not None


def build_hint_rows() -> list[dict[str, Any]]:
    rows = []
    for item in HINT_SOURCES:
        text = text_for(item["path"])
        features = {feature: False for feature in FEATURES}
        evidence = {}
        for feature, pattern in item["patterns"].items():
            matched = hit(pattern, text)
            features[feature] = matched
            evidence[feature] = pattern if matched else None
        rows.append({
            "source": item["source"],
            "path": str(item["path"].relative_to(REPO)),
            "feature_hits": features,
            "hit_count": sum(features.values()),
            "evidence_patterns": evidence,
            "exports_complete_selector_mechanism": False,
        })
    return rows


def glue_profiles(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    profiles = []
    n = len(rows)
    for mask in range(1, 1 << n):
        included = [rows[i] for i in range(n) if mask & (1 << i)]
        union = {feature: any(row["feature_hits"][feature] for row in included) for feature in FEATURES}
        # Required all hints plus at least one actual exported mechanism row.  Current rows are hints/obstructions only.
        profiles.append({
            "sources": [row["source"] for row in included],
            "covered_features": union,
            "covered_feature_count": sum(union.values()),
            "feature_complete": all(union.values()),
            "has_exported_mechanism_row": any(row["exports_complete_selector_mechanism"] for row in included),
            "accepted_as_selector_mechanism": all(union.values()) and any(row["exports_complete_selector_mechanism"] for row in included),
        })
    return profiles


def build_matrix() -> dict[str, Any]:
    rows = build_hint_rows()
    profiles = glue_profiles(rows)
    feature_coverage = {feature: sum(1 for row in rows if row["feature_hits"][feature]) for feature in FEATURES}
    obligations = [
        {"obligation": "content_hints_not_human_schema_assumption", "satisfied": True, "detail": "the object is a typed hint sheaf, not a preselected familiar selector formula"},
        {"obligation": "finite_hint_sources_scanned", "satisfied": True, "detail": "V1, H42, P3035, and P3036 are scanned as current hint/obstruction sources"},
        {"obligation": "all_required_hint_features_covered_by_union", "satisfied": any(p["feature_complete"] for p in profiles), "detail": "the union covers all six hint features"},
        {"obligation": "one_row_exports_nonpremise_selector_mechanism", "satisfied": any(row["exports_complete_selector_mechanism"] for row in rows), "detail": "all rows are hints or no-go obstructions, not exported mechanisms"},
        {"obligation": "glued_profile_accepted", "satisfied": any(p["accepted_as_selector_mechanism"] for p in profiles), "detail": "feature coverage without an exported mechanism row is not selector closure"},
        {"obligation": "readout_and_unit_coupling_present", "satisfied": feature_coverage["readout_unit_coupling_needed"] > 0, "detail": "P3036 names the coupling need but exports no unit theorem"},
    ]
    return {
        "object": "UnknownSelectorMechanism_HintSheafAcceptanceMatrix",
        "assumption": "human schemas may miss the actual selector mechanism; existing artifacts are treated as typed hints rather than closure claims",
        "required_features": FEATURES,
        "hint_rows": rows,
        "glue_profiles": profiles,
        "feature_coverage": feature_coverage,
        "proof_obligations": obligations,
        "finite_certificate": {
            "hint_source_rows": len(rows),
            "required_features": len(FEATURES),
            "features_with_some_coverage": sum(1 for count in feature_coverage.values() if count > 0),
            "glue_profiles": len(profiles),
            "feature_complete_profiles": sum(1 for p in profiles if p["feature_complete"]),
            "accepted_selector_mechanism_profiles": sum(1 for p in profiles if p["accepted_as_selector_mechanism"]),
            "exported_mechanism_rows": sum(1 for row in rows if row["exports_complete_selector_mechanism"]),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "proof_obligations": len(obligations),
            "strict_selector_mechanism_exported": any(p["accepted_as_selector_mechanism"] for p in profiles),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3035)
    read_json(P3036)
    matrix = build_matrix()
    return {
        "status": "P3037_SELECTOR_MECHANISM_HINT_SHEAF_ACCEPTANCE_MATRIX_NO_SELECTOR_EXPORT",
        "input_hashes": {
            "V1": hashlib.sha256(V1.read_bytes()).hexdigest() if V1.exists() else None,
            "H42": hashlib.sha256(H42.read_bytes()).hexdigest() if H42.exists() else None,
            "P3035": hashlib.sha256(P3035.read_bytes()).hexdigest() if P3035.exists() else None,
            "P3036": hashlib.sha256(P3036.read_bytes()).hexdigest() if P3036.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The repository contains real selector hints outside a single familiar schema: viscosity/damping/memory, c-retardation with anisotropic anchor dependence, inversion-odd/chiral branch-source needs, and external unit/readout-coupling needs. The finite hint sheaf covers all six required feature atoms in some glue profiles, but zero source rows export a non-premise complete selector mechanism; therefore feature coverage is not selector closure.",
            "negative_export_flags": {k: False for k in ["strict_selector_mechanism_exported", "qw2191_discharged", "strict_selector_branch_source_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not force the selector into a familiar prechosen schema, but also do not promote hints to closure. The next proof-grade move should construct one new integrated candidate operator combining memory/viscosity, retardation anisotropy, an inversion-odd signed value, and explicit unit/readout coupling, then run this sheaf acceptance test on that single candidate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3037/S1987 selector-mechanism hint sheaf acceptance matrix", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- hint source rows: `{c['hint_source_rows']}`",
        f"- required features: `{c['required_features']}`",
        f"- features with some coverage: `{c['features_with_some_coverage']}`",
        f"- glue profiles: `{c['glue_profiles']}`",
        f"- feature-complete profiles: `{c['feature_complete_profiles']}`",
        f"- accepted selector-mechanism profiles: `{c['accepted_selector_mechanism_profiles']}`",
        f"- exported mechanism rows: `{c['exported_mechanism_rows']}`",
        f"- strict selector mechanism exported: `{c['strict_selector_mechanism_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3037/S1987 selector-mechanism hint sheaf acceptance matrix", "## P3037/S1987 selector-mechanism hint sheaf acceptance matrix\n\n`P3037/S1987` assumes the real selector mechanism may not match a familiar human schema and treats existing research as typed hints.  It constructs an `UnknownSelectorMechanism_HintSheafAcceptanceMatrix` from V1 informational-viscosity hints, H42 retardation/anisotropy hints, P3035 inversion/chiral branch-source needs, and P3036 unit/readout-coupling needs.  The finite sheaf has four hint rows, six required features, and fifteen glue profiles; some profiles cover all features, but zero rows export a non-premise complete selector mechanism and zero profiles are accepted.  Therefore no `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3037/S1987 selector-hint sheaf `L_total` guard", "## P3037/S1987 selector-hint sheaf `L_total` guard\n\n`P3037/S1987` adds no physical `L_total` term.  The hint sheaf covers viscosity/memory, retardation anisotropy, inversion-odd/chiral sign needs, localizer/coupling needs, and unit/readout needs, but it installs no actual non-premise selector operator or variational source term.\n")
    append_once(AGENTS, "Current selector-mechanism hint sheaf guardrail (P3037/S1987, 2026-06-23)", "## Current selector-mechanism hint sheaf guardrail (P3037/S1987, 2026-06-23)\n\n- P3037 assumes the true selector mechanism may not fit an existing human schema and organizes current research as typed hints rather than closure claims.\n- The finite hint sheaf covers viscosity/damping/memory, c-retardation anisotropy, inversion-odd/chiral branch-source needs, source-localizer/coupling-polarity needs, and external unit/readout coupling needs, but zero rows export a complete non-premise selector mechanism.\n- Do not promote hint coverage, viscosity/retardation motifs, inversion-odd representation type, branch-source needs, or unit/readout coupling needs to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move may construct one integrated candidate operator combining memory/viscosity, retardation anisotropy, an inversion-odd signed value, and explicit unit/readout coupling, then rerun this acceptance test on that single candidate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
