#!/usr/bin/env python3
"""P2030 S980: strict tensor-resolved projection source audit.

Professor-level decision after P2029:

Before attempting a tensor-resolved B1 projection, audit whether the repo
already exports the minimal tensor data needed for that projection.  The answer
is intentionally allowed to be negative.  A negative audit is useful because it
prevents upgrading scalar B1 quotient evidence or ADM lapse-only witnesses into
full tensor-channel renormalization closure.
"""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2030_s980_strict_tensor_resolved_projection_source_audit.json"
MD = GEN / "p2030_s980_strict_tensor_resolved_projection_source_audit.md"

SCHEMA_VERSION = "p2030_s980_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
CHANNELS = ("R2", "Ric2", "Riem2", "GB")
COMPONENTS_REQUIRED = ("00", "11", "22", "33")


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def present(payload: dict[str, Any], key: str) -> bool:
    return key in payload and payload.get("_missing") is None


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1848 = load("p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json")
    p1981 = load("p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json")
    p1982 = load("p1982_s932_strict_adm_bianchi_ricci2_lapse_variation_obligation.json")
    p1983 = load("p1983_s933_strict_adm_bianchi_riemann2_lapse_variation_obligation.json")
    p1984 = load("p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json")
    p1985 = load("p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json")
    p1988 = load("p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json")
    p1991 = load("p1991_s941_strict_augmented_provider_channel_matrix_witness.json")
    p2029 = load("p2029_s979_strict_task1_renormalization_quotient_ledger_update.json")

    h_pack = p1848.get("gravity_componentwise_variation_pack") or {}
    scalar_profiles = (p1848.get("gravity_operator_profiles_B1") or {}).get("profiles") or {}
    h_templates_present = all(key in h_pack for key in ("H_R2_munu", "H_Ric2_munu", "H_Riem2_munu", "H_GB_munu"))
    scalar_profiles_present = all(key in scalar_profiles for key in CHANNELS)

    adm_lapse_sources = [
        {
            "packet": "P1981",
            "channel": "R2",
            "result_kind": p1981.get("result_kind"),
            "has_lapse_operator": "r2_lapse_euler_operator" in p1981,
            "scope": "ADM/Bianchi-I lapse Euler operator",
        },
        {
            "packet": "P1982",
            "channel": "Ric2",
            "result_kind": p1982.get("result_kind"),
            "has_lapse_operator": "ricci2_lapse_euler_operator" in p1982,
            "scope": "ADM/Bianchi-I lapse Euler operator",
        },
        {
            "packet": "P1983",
            "channel": "Riem2",
            "result_kind": p1983.get("result_kind"),
            "has_lapse_operator": "riemann2_lapse_euler_operator" in p1983,
            "scope": "ADM/Bianchi-I lapse Euler operator",
        },
        {
            "packet": "P1984",
            "channel": "GB",
            "result_kind": p1984.get("result_kind"),
            "has_lapse_operator": "gb_lapse_euler_operator" in p1984,
            "scope": "ADM/Bianchi-I lapse Euler cancellation",
        },
    ]
    adm_lapse_chain_present = all(row["has_lapse_operator"] for row in adm_lapse_sources)
    gb_lapse_cancels = bool((p1984.get("gb_lapse_euler_operator") or {}).get("EL_N_GB_difference_is_zero") is True)
    non_gb_lapse_obstruction_present = (
        p1985.get("result_kind") == "PASS_NON_GB_LAPSE_OBSTRUCTION_WITNESS"
        and bool((p1985.get("weighted_non_gb_lapse_operator") or {}).get("is_identically_zero") is False)
    )
    spatial_family_obstruction_present = p1988.get("result_kind") == "PASS_NON_GB_SPATIAL_PROJECTION_OBSTRUCTION_WITNESS"
    provider_matrix_obstruction_present = p1991.get("result_kind") == "PASS_CHANNEL_MATRIX_OBSTRUCTION_WITNESS"

    # Strict requirement for a tensor-resolved B1 projection.  The current repo
    # exports H_{mu nu} templates and several reduced/background witnesses, but
    # not this full channel x component profile table.
    required_tensor_object = {
        "object_id": "strict_B1_tensor_component_profile_table_v1",
        "required_channels": list(CHANNELS),
        "required_components": list(COMPONENTS_REQUIRED),
        "required_entries": [
            f"{channel}.H_{component}_profile(d)"
            for channel in CHANNELS
            for component in COMPONENTS_REQUIRED
        ],
        "required_metadata": [
            "background metric and gauge convention",
            "component inner-product/Gram rule",
            "UV endpoint quadrature rule for every component",
            "GB identity row and topological/null classification",
            "one-loop divergence tensor target vector in the same component basis",
        ],
    }

    exported_tensor_component_table_present = False
    exported_component_gram_present = False
    exported_divergence_tensor_target_present = False
    p1848_has_only_templates_not_components = h_templates_present and not exported_tensor_component_table_present
    adm_is_lapse_not_full_tensor = adm_lapse_chain_present and not exported_tensor_component_table_present

    source_rows = [
        {
            "source": "P1848",
            "available": h_templates_present,
            "level": "covariant_H_munu_template",
            "useful_for": "symbolic operator definitions and strict index convention",
            "insufficient_for": "component-profile Gram matrix; no H_00/H_11/H_22/H_33 B1 profiles exported",
        },
        {
            "source": "P1848.gravity_operator_profiles_B1",
            "available": scalar_profiles_present,
            "level": "scalar_B1_profile",
            "useful_for": "rank-3 quotient scalar projection from P2027/P2028",
            "insufficient_for": "tensor-resolved separation of GB across components",
        },
        {
            "source": "P1981/P1982/P1983/P1984",
            "available": adm_lapse_chain_present,
            "level": "ADM_BianchiI_lapse_Euler_chain",
            "useful_for": "lapse-level GB cancellation and non-GB residual diagnosis",
            "insufficient_for": "full H_munu component table and tensor inner product",
        },
        {
            "source": "P1985",
            "available": non_gb_lapse_obstruction_present,
            "level": "weighted_non_GB_lapse_obstruction",
            "useful_for": "proving strict non-GB anisotropic lapse residual persists",
            "insufficient_for": "four-component tensor projection closure",
        },
        {
            "source": "P1988",
            "available": spatial_family_obstruction_present,
            "level": "spatial_family_projection_obstruction",
            "useful_for": "showing non-GB families exceed EH anisotropic capacity",
            "insufficient_for": "channel-resolved R2/Ric2/Riem2/GB spatial components",
        },
        {
            "source": "P1991",
            "available": provider_matrix_obstruction_present,
            "level": "augmented_provider_channel_matrix_obstruction",
            "useful_for": "testing non-strict augmented provider insufficiency",
            "insufficient_for": "strict tensor source or selector-safe closure",
        },
        {
            "source": "P2029",
            "available": p2029.get("result_kind") == "TASK1_LEDGER_UPDATED_TO_RANK3_QUOTIENT_PASS_TENSOR_EXTENSION_OPEN",
            "level": "task1_quotient_ledger",
            "useful_for": "preventing scalar quotient pass from becoming four-channel uniqueness",
            "insufficient_for": "new tensor data",
        },
    ]

    minimal_requirements = [
        {
            "requirement": "covariant_H_munu_templates",
            "present": h_templates_present,
            "source": "P1848",
        },
        {
            "requirement": "scalar_B1_profiles",
            "present": scalar_profiles_present,
            "source": "P1848/P2027",
        },
        {
            "requirement": "ADM_lapse_channel_chain",
            "present": adm_lapse_chain_present,
            "source": "P1981-P1984",
        },
        {
            "requirement": "tensor_component_profile_table",
            "present": exported_tensor_component_table_present,
            "source": "MISSING",
        },
        {
            "requirement": "tensor_component_Gram_rule",
            "present": exported_component_gram_present,
            "source": "MISSING",
        },
        {
            "requirement": "same_basis_divergence_tensor_target",
            "present": exported_divergence_tensor_target_present,
            "source": "MISSING",
        },
    ]
    blocking_requirements = [row["requirement"] for row in minimal_requirements if not row["present"]]
    tensor_projection_ready = not blocking_requirements

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2030",
        "stage_id": "S980",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "OPEN_TENSOR_RESOLVED_SOURCE_GAP_WITH_TRACE",
        "route": "strict_only",
        "strict_lane_assumptions": [
            "strict_kernel_only",
            "no_legacy_transfer",
            "scalar_B1_quotient_pass_not_tensor_projection",
            "ADM_lapse_witness_not_full_H_munu_component_table",
        ],
        "depends_on": {
            "p1848_present": p1848.get("_missing") is None,
            "p1981_present": p1981.get("_missing") is None,
            "p1982_present": p1982.get("_missing") is None,
            "p1983_present": p1983.get("_missing") is None,
            "p1984_present": p1984.get("_missing") is None,
            "p1985_present": p1985.get("_missing") is None,
            "p1988_present": p1988.get("_missing") is None,
            "p1991_present": p1991.get("_missing") is None,
            "p2029_present": p2029.get("_missing") is None,
        },
        "input_hashes": {
            "p1848_json_sha256": file_sha256(GEN / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json"),
            "p1984_json_sha256": file_sha256(GEN / "p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json"),
            "p1985_json_sha256": file_sha256(GEN / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json"),
            "p1988_json_sha256": file_sha256(GEN / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json"),
            "p1991_json_sha256": file_sha256(GEN / "p1991_s941_strict_augmented_provider_channel_matrix_witness.json"),
            "p2029_json_sha256": file_sha256(GEN / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.json"),
        },
        "source_audit_rows": source_rows,
        "adm_lapse_sources": adm_lapse_sources,
        "known_positive_evidence": {
            "covariant_H_munu_templates_present": h_templates_present,
            "scalar_B1_profiles_present": scalar_profiles_present,
            "adm_lapse_chain_present": adm_lapse_chain_present,
            "gb_lapse_cancels": gb_lapse_cancels,
            "non_gb_lapse_obstruction_present": non_gb_lapse_obstruction_present,
            "spatial_family_obstruction_present": spatial_family_obstruction_present,
            "provider_matrix_obstruction_present": provider_matrix_obstruction_present,
        },
        "required_tensor_object": required_tensor_object,
        "minimal_requirements": minimal_requirements,
        "blocking_requirements": blocking_requirements,
        "tensor_projection_ready": tensor_projection_ready,
        "professor_decision": {
            "decision": "DO_NOT_ATTEMPT_TENSOR_PROJECTION_UNTIL_COMPONENT_TABLE_EXISTS",
            "rationale": [
                "P1848 exports covariant templates and scalar B1 profiles, not componentwise tensor profiles.",
                "P1981-P1984 prove valuable ADM/Bianchi-I lapse facts, including GB lapse cancellation, but lapse-only is not full H_munu.",
                "P1988/P1991 expose spatial/provider obstructions rather than a strict tensor source table.",
            ],
            "next_route_preference": "construct_minimal_strict_B1_tensor_component_profile_table_v1",
        },
        "gatekeeper_checks": {
            "h_templates_present": h_templates_present,
            "scalar_profiles_present": scalar_profiles_present,
            "adm_lapse_chain_present": adm_lapse_chain_present,
            "p1848_has_only_templates_not_components": p1848_has_only_templates_not_components,
            "adm_is_lapse_not_full_tensor": adm_is_lapse_not_full_tensor,
            "tensor_component_profile_table_missing": not exported_tensor_component_table_present,
            "component_gram_missing": not exported_component_gram_present,
            "divergence_tensor_target_missing": not exported_divergence_tensor_target_present,
            "no_tensor_projection_claimed": not tensor_projection_ready,
            "no_toe_closure_claimed": True,
        },
        "theorem_scope": {
            "licensed_statement": "The repo currently has scalar B1 quotient evidence and ADM/Bianchi-I lapse/operator obstruction evidence, but it does not yet export the minimal tensor-component profile table needed for a tensor-resolved B1 projection.",
            "not_licensed": [
                "tensor-resolved Gram rank or coefficient identification",
                "independent a_GB identification",
                "four-channel counterterm uniqueness",
                "background-global renormalization",
                "BRST/Cutkosky unitarity closure",
                "QW-2191 selector closure",
                "ToE closure",
            ],
        },
        "false_pass_guard": "ADM lapse cancellation of GB and scalar B1 quotient identifiability are real, but neither is a full tensor-resolved H_munu projection.",
        "next_honest_step": "Build P2031: minimal strict_B1_tensor_component_profile_table_v1 scaffold with explicit missing entries for each channel/component, then fill only entries that can be derived from existing strict sources.",
        "lay_explanation": "Mamy wzory tensorowe, skalarne profile B1 oraz mocne rachunki lapse dla Bianchi-I. Brakuje jednak tabeli komponentow tensorowych dla kazdego kanalu. Bez niej nie wolno twierdzic, ze GB zostal rozdzielony tensorowo.",
        "environment": {
            "python": platform.python_version(),
        },
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = f"""# P2030 S980 Strict Tensor-Resolved Projection Source Audit

Status: `OPEN_OBSTRUCTION_WITH_TRACE`

Result kind: `OPEN_TENSOR_RESOLVED_SOURCE_GAP_WITH_TRACE`

## Professor Decision

`DO_NOT_ATTEMPT_TENSOR_PROJECTION_UNTIL_COMPONENT_TABLE_EXISTS`

The repo has:
- covariant `H_munu` templates from P1848;
- scalar B1 profiles from P1848/P2027;
- ADM/Bianchi-I lapse witnesses from P1981-P1984;
- non-GB lapse/spatial/provider obstruction witnesses from P1985/P1988/P1991.

But the repo does **not** yet export:
- a `channel x component` tensor profile table for `R2/Ric2/Riem2/GB`;
- a tensor-component Gram rule;
- a same-basis divergence tensor target vector.

Therefore tensor projection ready: `{tensor_projection_ready}`.

Blocking requirements:
{chr(10).join(f"- {item}" for item in blocking_requirements)}

## Minimal Required Object

`strict_B1_tensor_component_profile_table_v1`

Required components: `{', '.join(COMPONENTS_REQUIRED)}`.

Required channels: `{', '.join(CHANNELS)}`.

## Honest Interpretation

P1984's GB lapse cancellation is valuable evidence for topological behavior in
ADM/Bianchi-I minisuperspace, and P2028's scalar quotient theorem is valuable
Task-1 progress.  Neither supplies the full tensor-resolved projection needed
to identify `a_GB` or claim four-channel counterterm uniqueness.
"""
    MD.write_text(md, encoding="utf-8")
    print(OUT)
    print(MD)


if __name__ == "__main__":
    main()
