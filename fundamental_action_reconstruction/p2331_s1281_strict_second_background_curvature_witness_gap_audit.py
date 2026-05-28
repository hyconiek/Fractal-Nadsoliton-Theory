#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.json"
MD = GEN / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.md"

SOURCE_FILES = {
    "P1985_ADM_BIANCHI_NON_GB_LAPSE": GEN / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json",
    "P2030_TENSOR_RESOLVED_SOURCE_AUDIT": GEN / "p2030_s980_strict_tensor_resolved_projection_source_audit.json",
    "P2031_TENSOR_COMPONENT_TABLE_SCAFFOLD": GEN / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.json",
    "P2033_CURVED_B1_ANSATZ_NONAVAILABILITY": GEN / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.json",
    "P2034_QUOTIENT_ONLY_RENORMALIZATION": GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json",
    "P2036_BACKGROUND_TRANSPORT_CONTRACT": GEN / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.json",
    "P2296_GLOBAL_TASK1_REPLAY": GEN / "p2296_s1246_strict_global_task1_renormalization_replay_and_7task_reclassification_probe.json",
    "P2297_NON_GB_SPATIAL_EOM_PROVIDER": GEN / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json",
    "P2330_B1_GB_DEPENDENCE": GEN / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.json",
}

GREP_PATTERNS = (
    "second background",
    "tensor-resolved curvature",
    "tensor_component_profile_table",
    "background-global renormalization",
    "Gauss-Bonnet",
    "GB null",
    "quotient-only",
    "P1985",
    "P2297",
    "globalization",
)

REQUIRED_LIFT_FIELDS = (
    "independent_background_family",
    "full_tensor_component_profiles",
    "gb_channel_in_same_basis",
    "component_gram_rule",
    "same_basis_divergence_target",
    "transport_or_globalization_theorem",
)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": path.relative_to(REPO).as_posix()}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def grep_hits() -> list[dict[str, Any]]:
    candidates = [
        ROOT / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.py",
        ROOT / "p2030_s980_strict_tensor_resolved_projection_source_audit.py",
        ROOT / "p2031_s981_strict_b1_tensor_component_profile_table_scaffold.py",
        ROOT / "p2033_s983_strict_curved_b1_metric_ansatz_nonavailability_theorem.py",
        ROOT / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.py",
        ROOT / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.py",
        ROOT / "p2296_s1246_strict_global_task1_renormalization_replay_and_7task_reclassification_probe.py",
        ROOT / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.py",
        ROOT / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.py",
    ]
    candidates.extend(SOURCE_FILES.values())
    rows: list[dict[str, Any]] = []
    for path in candidates:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        excerpt = ""
        for index, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = index
                excerpt = line.strip()[:220]
                break
        rows.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": excerpt,
        })
    rows.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return rows


def source_capability_rows(artifacts: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    p1985 = artifacts["P1985_ADM_BIANCHI_NON_GB_LAPSE"]
    p2030 = artifacts["P2030_TENSOR_RESOLVED_SOURCE_AUDIT"]
    p2031 = artifacts["P2031_TENSOR_COMPONENT_TABLE_SCAFFOLD"]
    p2033 = artifacts["P2033_CURVED_B1_ANSATZ_NONAVAILABILITY"]
    p2036 = artifacts["P2036_BACKGROUND_TRANSPORT_CONTRACT"]
    p2297 = artifacts["P2297_NON_GB_SPATIAL_EOM_PROVIDER"]
    return [
        {
            "source_id": "P1985_ADM_BIANCHI_NON_GB_LAPSE",
            "status": p1985.get("status"),
            "independent_background_family": True,
            "full_tensor_component_profiles": False,
            "gb_channel_in_same_basis": False,
            "component_gram_rule": False,
            "same_basis_divergence_target": False,
            "transport_or_globalization_theorem": False,
            "positive_evidence": "ADM/Bianchi-I non-GB lapse residual is nonzero and strict-only.",
            "blocking_reason": "lapse-only and non-GB; does not export a full tensor component basis or GB row in the same basis",
            "numeric_indicator": {"numeric_replay_l2_norm": p1985.get("numeric_replay_l2_norm")},
        },
        {
            "source_id": "P2297_NON_GB_SPATIAL_EOM_PROVIDER",
            "status": p2297.get("status"),
            "independent_background_family": True,
            "full_tensor_component_profiles": False,
            "gb_channel_in_same_basis": False,
            "component_gram_rule": False,
            "same_basis_divergence_target": False,
            "transport_or_globalization_theorem": False,
            "positive_evidence": "tracefree spatial Bianchi-I provider matrix obstruction is exported.",
            "blocking_reason": "strict-core provider matrix is inconsistent and formal full-basis reconstruction is marked non-admissible",
            "numeric_indicator": p2297.get("gatekeeper_checks", {}),
        },
        {
            "source_id": "P2030_TENSOR_RESOLVED_SOURCE_AUDIT",
            "status": p2030.get("status"),
            "independent_background_family": False,
            "full_tensor_component_profiles": False,
            "gb_channel_in_same_basis": False,
            "component_gram_rule": False,
            "same_basis_divergence_target": False,
            "transport_or_globalization_theorem": False,
            "positive_evidence": ", ".join(p2030.get("known_positive_evidence", [])),
            "blocking_reason": "tensor_projection_ready is false; tensor component table, Gram rule, and divergence tensor target are missing",
            "numeric_indicator": {"tensor_projection_ready": p2030.get("tensor_projection_ready")},
        },
        {
            "source_id": "P2031_TENSOR_COMPONENT_TABLE_SCAFFOLD",
            "status": p2031.get("status"),
            "independent_background_family": False,
            "full_tensor_component_profiles": False,
            "gb_channel_in_same_basis": False,
            "component_gram_rule": False,
            "same_basis_divergence_target": False,
            "transport_or_globalization_theorem": False,
            "positive_evidence": "4x4 tensor table scaffold exists with conditional GB identity rows.",
            "blocking_reason": "all 16 required tensor entries are still marked missing and no component profile is derived",
            "numeric_indicator": p2031.get("table_summary", {}),
        },
        {
            "source_id": "P2033_CURVED_B1_ANSATZ_NONAVAILABILITY",
            "status": p2033.get("status"),
            "independent_background_family": False,
            "full_tensor_component_profiles": False,
            "gb_channel_in_same_basis": False,
            "component_gram_rule": False,
            "same_basis_divergence_target": False,
            "transport_or_globalization_theorem": False,
            "positive_evidence": "formal nonavailability theorem for the currently missing curved B1 ansatz/projection rule.",
            "blocking_reason": "minimal curved B1 ansatz and component projection rule are not exported",
            "numeric_indicator": p2033.get("gatekeeper_checks", {}),
        },
        {
            "source_id": "P2036_BACKGROUND_TRANSPORT_CONTRACT",
            "status": p2036.get("status"),
            "independent_background_family": False,
            "full_tensor_component_profiles": False,
            "gb_channel_in_same_basis": False,
            "component_gram_rule": False,
            "same_basis_divergence_target": False,
            "transport_or_globalization_theorem": False,
            "positive_evidence": "candidate contract for background transport exists.",
            "blocking_reason": "contract is a new-premise scaffold; global transport theorem and finite-part scheme transport are not exported",
            "numeric_indicator": p2036.get("gatekeeper_checks", {}),
        },
    ]


def coverage_matrix(rows: list[dict[str, Any]]) -> dict[str, Any]:
    mat = np.array([[1.0 if row[field] else 0.0 for field in REQUIRED_LIFT_FIELDS] for row in rows], dtype=float)
    coverage = {field: bool(np.any(mat[:, idx] > 0.5)) for idx, field in enumerate(REQUIRED_LIFT_FIELDS)}
    all_fields_met_by_union = all(coverage.values())
    admissible_rows = [
        row["source_id"] for row in rows
        if all(bool(row[field]) for field in REQUIRED_LIFT_FIELDS)
    ]
    return {
        "required_fields": list(REQUIRED_LIFT_FIELDS),
        "matrix_rows": [row["source_id"] for row in rows],
        "binary_matrix": mat.tolist(),
        "field_coverage_by_current_union": coverage,
        "covered_field_count": int(sum(coverage.values())),
        "required_field_count": len(REQUIRED_LIFT_FIELDS),
        "all_fields_met_by_current_union": all_fields_met_by_union,
        "admissible_single_source_rows": admissible_rows,
        "numeric_rank_of_capability_matrix": int(np.linalg.matrix_rank(mat)) if mat.size else 0,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    rows = source_capability_rows(artifacts)
    cap = coverage_matrix(rows)
    p2031_summary = artifacts["P2031_TENSOR_COMPONENT_TABLE_SCAFFOLD"].get("table_summary", {})
    p2330_cert = artifacts["P2330_B1_GB_DEPENDENCE"].get(
        "strict_b1_renormalization_gb_dependence_globalization_obstruction_probe", {}
    ).get("gb_dependence_certificate", {})

    # Formal current-export implication: if the tensor table is empty and the
    # capability union does not cover all lift fields, then no current export
    # can lift the B1 GB quotient dependence to a global 4-channel theorem.
    missing_tensor_entries = int(p2031_summary.get("missing_entry_count", 16))
    total_tensor_entries = int(p2031_summary.get("total_required_entries", 16))
    quotient_only_expr = sp.And(
        sp.Eq(sp.Symbol("missing_tensor_entries"), total_tensor_entries),
        sp.Not(sp.Symbol("all_lift_fields_met")),
    )
    quotient_only_current_exports = missing_tensor_entries == total_tensor_entries and not cap["all_fields_met_by_current_union"]

    current_backend_scope_theorem = {
        "theorem_name": "P2331 current-export quotient-scope renormalization boundary",
        "claim": "The repo already contains independent ADM/Bianchi-I obstruction evidence, but no current export supplies the full tensor-component profiles, GB row in the same basis, component Gram rule, same-basis divergence target, and transport theorem required to lift the B1 Gauss-Bonnet quotient dependence. Therefore, on current exports, strict Task-1 renormalization remains quotient-scope only.",
        "proof_witnesses": [
            "P1985/P2297 provide independent ADM/Bianchi-I non-GB lapse/spatial obstruction evidence, not a full tensor-resolved curvature witness.",
            "P2030 says tensor_projection_ready is false and names the missing tensor component table, Gram rule, and divergence target.",
            "P2031 reports 16/16 required tensor component entries missing and no derived component profile.",
            "P2033 proves current nonavailability of the curved B1 metric ansatz/component projection rule.",
            "P2330 shows the current B1 GB channel has rank-3/nullity-1 quotient dependence.",
        ],
        "formal_current_export_condition": str(quotient_only_expr),
        "formal_current_export_condition_evaluates_true": bool(quotient_only_current_exports),
        "scope_limits": [
            "This is a current-export boundary theorem, not an eternal no-go theorem.",
            "A future tensor-resolved background-family witness may reopen the lift test.",
            "No full/global renormalization, QW-2191 discharge, G1/G3 update, or ToE closure is claimed.",
        ],
        "strict_guardrails": {
            "strict_kernel_only": True,
            "no_legacy_kernel_role_transfer": True,
            "no_selector_premise_added": True,
            "no_qw2191_discharge_claimed": True,
            "no_g1_g3_update_claimed": True,
            "no_toe_closure_claimed": True,
        },
    }

    hits = grep_hits()
    probe = {
        "probe_id": "P2331_S1281_STRICT_SECOND_BACKGROUND_CURVATURE_WITNESS_GAP_AUDIT",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_patterns": list(GREP_PATTERNS),
            "hit_count": len(hits),
            "top_hits": hits[:20],
        },
        "second_background_candidate_rows": rows,
        "lift_capability_matrix": cap,
        "p2031_tensor_table_summary": p2031_summary,
        "p2330_gb_dependence_summary": {
            "numeric_rank": p2330_cert.get("numeric_rank"),
            "numeric_nullity": p2330_cert.get("numeric_nullity"),
            "gram_null_residual_l2": p2330_cert.get("gram_null_residual_l2"),
        },
        "current_backend_scope_theorem": current_backend_scope_theorem,
        "theorem_fingerprint_sha256": sha256_json(current_backend_scope_theorem),
    }

    gatekeeper_checks = {
        "grep_hits_found": len(hits) >= 5,
        "p1985_loaded": artifacts["P1985_ADM_BIANCHI_NON_GB_LAPSE"].get("packet_id") == "P1985",
        "p2297_loaded": artifacts["P2297_NON_GB_SPATIAL_EOM_PROVIDER"].get("packet_id") == "P2297",
        "p2030_tensor_projection_not_ready": artifacts["P2030_TENSOR_RESOLVED_SOURCE_AUDIT"].get("tensor_projection_ready") is False,
        "p2031_all_tensor_entries_missing": missing_tensor_entries == total_tensor_entries,
        "p2033_ansatz_nonavailable": artifacts["P2033_CURVED_B1_ANSATZ_NONAVAILABILITY"].get("gatekeeper_checks", {}).get("nonavailability_theorem_passed") is True,
        "p2330_rank3_nullity1_preserved": p2330_cert.get("numeric_rank") == 3 and p2330_cert.get("numeric_nullity") == 1,
        "no_admissible_second_background_lift_source": len(cap["admissible_single_source_rows"]) == 0,
        "current_union_does_not_cover_lift_fields": cap["all_fields_met_by_current_union"] is False,
        "current_exports_only_quotient_scope": quotient_only_current_exports,
        "no_full_global_renormalization_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2331_s1281_v1",
        "packet_id": "P2331",
        "stage_id": "S1281",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_CURRENT_EXPORTS_ONLY_QUOTIENT_SCOPE_RENORMALIZATION_BOUNDARY",
        "result_kind": "STRICT_SECOND_BACKGROUND_CURVATURE_WITNESS_GAP_AUDIT_NO_GLOBAL_RENORMALIZATION_CLAIM",
        "strict_second_background_curvature_witness_gap_audit": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": "Build the missing tensor-resolved curvature component table and same-basis divergence target for at least one independent background-family witness, then rerun the GB-lift rank test. Until then, current strict exports license only quotient-scope renormalization.",
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2331 second-background curvature witness gap audit\n\n"
        "Status: current strict exports still license only quotient-scope renormalization.\n\n"
        f"- Candidate source rows audited: `{len(rows)}`.\n"
        f"- Required lift fields covered by union: `{cap['covered_field_count']}/{cap['required_field_count']}`.\n"
        f"- P2031 missing tensor entries: `{missing_tensor_entries}/{total_tensor_entries}`.\n"
        f"- P2330 GB rank/nullity: `{p2330_cert.get('numeric_rank')}/{p2330_cert.get('numeric_nullity')}`.\n"
        "- No full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
