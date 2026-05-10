#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F704 = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_P705 = GENERATED / "p705_current_actual_strict_projective_operational_toe_os_support_with_invariant_mass_observable_probe_summary.json"
IN_Q2068 = REPO / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json" / "report_qw2068_sm_gr_parameter_registry.json"
IN_Q2069 = REPO / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json" / "report_qw2069_full_sm_gr_derivation_package.json"
IN_Q2084 = REPO / "report_qw2084_t1_nonanchor_strict_gate.json"
IN_Q2087 = REPO / "report_qw2087_alpha_s_nonanchor_boundary_gate.json"
IN_Q2093 = REPO / "report_qw2093_kernel_derived_nonanchor_inputs_plan_executor.json"
IN_Q2103 = REPO / "report_qw2103_gnewton_dimensionless_provenance_gate.json"
IN_Q2113 = REPO / "report_qw2113_gnewton_direct_dimensionless_pack_gate.json"
IN_GNEWTON_READY = REPO / "external_gnewton_bridge_qw2101.direct_dimensionless_ready.json"
IN_A8 = GENERATED / "a8_gravity_bridge_summary.json"

OUT = GENERATED / "f784_current_minimal_strict_sm_gr_bridge_registry_packet.json"
OUT_SUMMARY = GENERATED / "f784_current_minimal_strict_sm_gr_bridge_registry_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def first_group_entry_by_id(registry: dict[str, Any], group_name: str, entry_id: str) -> dict[str, Any]:
    groups = registry.get("groups") or {}
    group = groups.get(group_name) or []
    if not isinstance(group, list):
        return {}
    for item in group:
        if isinstance(item, dict) and item.get("id") == entry_id:
            return item
    return {}


def first_entry_by_id(entries: list[Any], entry_id: str) -> dict[str, Any]:
    for item in entries:
        if isinstance(item, dict) and item.get("id") == entry_id:
            return item
    return {}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_F704,
        IN_P705,
        IN_Q2068,
        IN_Q2069,
        IN_Q2084,
        IN_Q2087,
        IN_Q2093,
        IN_Q2103,
        IN_Q2113,
        IN_GNEWTON_READY,
        IN_A8,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F784",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f704 = load_json(IN_F704)
    p705 = load_json(IN_P705)
    q2068 = load_json(IN_Q2068)
    q2069 = load_json(IN_Q2069)
    q2084 = load_json(IN_Q2084)
    q2087 = load_json(IN_Q2087)
    q2093 = load_json(IN_Q2093)
    q2103 = load_json(IN_Q2103)
    q2113 = load_json(IN_Q2113)
    gnewton_ready = load_json(IN_GNEWTON_READY)
    a8 = load_json(IN_A8)

    m_proxy = ((f704.get("outputs") or {}).get("m_proxy_ascending")) or []
    if not (
        isinstance(m_proxy, list)
        and len(m_proxy) >= 2
        and all(isinstance(x, (int, float)) and math.isfinite(float(x)) and float(x) > 0.0 for x in m_proxy)
    ):
        artifact = {
            "stage": "F784",
            "status": "NOT_COMPUTABLE_INVALID_F704_MASS_PROXY_OBJECT",
            "as_of": AS_OF,
            "input": rel(IN_F704),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    m_proxy_f = [float(x) for x in m_proxy]
    adjacent_ratios = [m_proxy_f[i + 1] / m_proxy_f[i] for i in range(len(m_proxy_f) - 1)]

    q2069_entries = q2069.get("entries") or []
    if not isinstance(q2069_entries, list):
        q2069_entries = []

    q2068_sin2 = first_group_entry_by_id(q2068, "gauge_and_electroweak", "sin2_theta_w_mz")
    q2068_alpha_s = first_group_entry_by_id(q2068, "gauge_and_electroweak", "alpha_s_mz")

    q2069_sin2 = first_entry_by_id(q2069_entries, "sin2_theta_w_mz")
    q2069_alpha_s = first_entry_by_id(q2069_entries, "alpha_s_mz")
    q2069_gnewton = first_entry_by_id(q2069_entries, "G_newton")

    q2093_values = q2093.get("derived_values") or {}
    q2093_formulas = q2093.get("frozen_formulas") or {}
    q2113_contract = q2113.get("qw2108_contract") or {}
    q2113_range = q2113_contract.get("g_range") or {}

    targets = {
        "mass_ratio_ordering_layer": {
            "status_label": "strict-derived",
            "bridge_ready": True,
            "scope": "strict_internal_dimensionless_proxy_only",
            "external_observable_origin": False,
            "source_artifacts": [
                rel(IN_F704),
                rel(IN_P705),
                "fundamental_action_reconstruction/N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md",
                "fundamental_action_reconstruction/generated/p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json",
            ],
            "strict_input_chain": [
                "strict_projective_selector_closure_support_v3",
                "basis_invariant_mass_observable_F704",
                "N703_quadratic_mass_proxy_meaning_discipline",
            ],
            "observable_snapshot": {
                "n_modes": len(m_proxy_f),
                "m_proxy_ascending": m_proxy_f,
                "adjacent_ratio_ascending": adjacent_ratios,
                "span_ratio_max_over_min": m_proxy_f[-1] / m_proxy_f[0],
            },
            "promotion_blockers": [
                "No particle assignment is exported on this route.",
                "No physical-unit calibration belongs to this strict object.",
                "No host-matching or Standard Model identification claim follows from this layer alone.",
            ],
            "hard_limits": f704.get("hard_limits") or [],
        },
        "sin2_theta_w_eff": {
            "status_label": "strict-derived",
            "bridge_ready": True,
            "scope": "strict_internal_gate_candidate_observable_with_unresolved_legacy_semantics",
            "external_observable_origin": False,
            "source_artifacts": [
                rel(IN_Q2068),
                rel(IN_Q2069),
                rel(IN_Q2093),
                "fundamental_action_reconstruction/F9_LEGACY_WEINBERG_SEMANTIC_TRANSFER_REFINEMENT_PACKET.md",
            ],
            "strict_input_chain": [
                "QW-2063_shared_flavor_basis_reconstruction",
                "QW-2064_micro_derived_renormalization_constants",
                "QW-2093_frozen_nonanchor_formula_set",
                "QW-2098_nonanchor_ew_pole_chain",
                "QW-2069_full_sm_gr_derivation_package",
            ],
            "observable_snapshot": {
                "upstream_lineage_alias": "sin2_theta_w_eff",
                "promoted_object_id": "sin2_theta_w_mz",
                "lineage_value_from_q2093": float(q2093_values.get("sin2_theta_w_eff")),
                "promoted_value_from_q2069": q2069_sin2.get("predicted_value"),
                "reference_value_from_q2068_registry": q2068_sin2.get("value"),
                "lineage_formula": q2093_formulas.get("sin2_eff_kernel"),
                "promotion_method": q2069_sin2.get("method"),
                "promotion_notes": q2069_sin2.get("notes"),
            },
            "promotion_blockers": [
                "The strict-side observable is promoted in the registry/package lane, but the legacy Weinberg semantic-transfer verdict remains unresolved on current repo state.",
                "The upstream QW-2093 lineage still contains an explicit alpha_geo/12 touchpoint and therefore cannot be silently treated as a strict replacement of the legacy Weinberg-role semantics.",
                "No retained/replaced verdict for the legacy Weinberg role follows from this bridge entry alone.",
            ],
            "hard_limits": [
                "Do not treat this observable as discharge of the legacy Weinberg role.",
                "Do not promote the QW-2093 lineage touchpoint into retained semantics without an explicit verdict.",
                "Do not infer complete electroweak closure from this single promoted observable.",
            ],
        },
        "alpha_s_boundary_mu0_alpha0": {
            "status_label": "strict-derived",
            "bridge_ready": True,
            "scope": "strict_internal_gate_boundary_entry_with_frozen_ansatz_origin",
            "external_observable_origin": False,
            "source_artifacts": [
                rel(IN_Q2068),
                rel(IN_Q2069),
                rel(IN_Q2093),
                rel(IN_Q2087),
                rel(IN_Q2084),
            ],
            "strict_input_chain": [
                "QW-2063_shared_flavor_basis_reconstruction",
                "QW-2064_micro_derived_renormalization_constants",
                "QW-2093_frozen_nonanchor_formula_set",
                "QW-2087_nonanchor_boundary_gate",
                "QW-2069_full_sm_gr_derivation_package",
            ],
            "observable_snapshot": {
                "mu0_gev": float(q2093_values.get("alpha_s_boundary_mu0_gev")),
                "alpha0": float(q2093_values.get("alpha_s_boundary_alpha0")),
                "formula_mu0": q2093_formulas.get("alpha_s_boundary_mu0"),
                "formula_alpha0": q2093_formulas.get("alpha_s_boundary_alpha0"),
                "downstream_alpha_s_mz": q2069_alpha_s.get("predicted_value"),
                "reference_alpha_s_mz": q2068_alpha_s.get("value"),
                "downstream_method": q2069_alpha_s.get("method"),
                "boundary_gate_verdict": q2087.get("verdict"),
                "t1_gate_verdict": q2084.get("verdict"),
            },
            "promotion_blockers": [
                "The boundary source remains a frozen ansatz inside the strict pipeline, not a kernel-invariants-only bridge law.",
                "The route still depends on QW-2063 regime data that are not yet reduced to a kernel-invariants-only bridge registry layer.",
                "The strict promotion currently lives through downstream nonanchor gate support rather than through one standalone alpha_s boundary theorem.",
            ],
            "hard_limits": [
                "Do not treat this boundary object as full QCD closure.",
                "Do not treat downstream agreement as proof of kernel-alone derivation.",
                "Do not collapse alpha_s_mz downstream success into a claim that the boundary ansatz is uniquely forced.",
            ],
        },
        "g_dimensionless_mu_ref": {
            "status_label": "strict-derived-with-external-observable-origin",
            "bridge_ready": True,
            "scope": "strict_side_endpoint_with_explicit_external_origin_only",
            "external_observable_origin": True,
            "source_artifacts": [
                rel(IN_GNEWTON_READY),
                rel(IN_Q2103),
                rel(IN_Q2113),
                rel(IN_Q2069),
                rel(IN_A8),
            ],
            "strict_input_chain": [
                "QW-2113_direct_dimensionless_pack",
                "QW-2101_ready_external_dimensionless_payload",
                "QW-2103_provenance_gate",
                "QW-2092_dimensionless_to_si_bridge",
                "A8_partial_gravity_bridge_boundary",
            ],
            "observable_snapshot": {
                "mu_ref_gev": float(gnewton_ready.get("mu_ref_gev")),
                "g_dimensionless_mu_ref": float(gnewton_ready.get("g_dimensionless_mu_ref")),
                "bridge_observable_origin": gnewton_ready.get("bridge_observable_origin"),
                "acceptance_range_from_qw2108": {
                    "min": q2113_range.get("min"),
                    "max": q2113_range.get("max"),
                },
                "provenance_gate_verdict": q2103.get("verdict"),
                "pack_gate_verdict": q2113.get("verdict"),
                "downstream_gnewton_method": q2069_gnewton.get("method"),
            },
            "promotion_blockers": [
                "The bridge observable origin is explicitly external_dimensionless_observable, so this is not an internal strict derivation of G.",
                "A8 keeps the internal origin of the dimensionless G bridge observable explicitly open.",
                "No Einstein-Hilbert derivation, equivalence-principle derivation, or full GR closure follows from this endpoint alone.",
            ],
            "hard_limits": [
                "Do not relabel this external-origin endpoint as internally derived gravity closure.",
                "Do not infer full SM+GR reduction from this object.",
            ],
        },
    }

    counts = {
        "strict-derived": 0,
        "strict-derived-with-external-observable-origin": 0,
        "anchor-dependent": 0,
        "si-definition": 0,
        "non-strict": 0,
        "open": 0,
    }
    for item in targets.values():
        label = item["status_label"]
        counts[label] = counts.get(label, 0) + 1

    packet = {
        "stage": "F784",
        "packet_name": "MinimalStrictSMGRBridgeRegistry_v1",
        "status": "F784_EXECUTED_CURRENT_MINIMAL_STRICT_SM_GR_BRIDGE_REGISTRY_PACKET_NO_FALSE_PASS",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "minimal_strict_side_bridge_to_sm_gr_observables_only",
        "source_policy": {
            "strict_inputs_only": [
                "K_strict_gate",
                "strict-side nadsoliton ontology lane",
                "strict-side alpha_geo object when explicitly exported",
                "declared strict premises only",
            ],
            "forbidden_shortcuts": [
                "silent legacy-to-strict role transfer",
                "hidden GeV calibration inside strict bridge status",
                "host matching treated as strict derivation",
                "external gravity endpoint relabeled as internal origin of G",
            ],
        },
        "prerequisite_checks": {
            "basis_invariant_mass_observable_exported": bool(
                f704.get("status") == "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT"
            ),
            "strict_projective_os_support_v3_present": bool(
                p705.get("strict_projective_operational_toe_os_support_packet_v3_exported") is True
            ),
            "qw2069_full_sm_gr_package_present": bool(q2069.get("verdict") == "FULL_SM_GR_DERIVATION_PACKAGE_PASS"),
            "qw2093_frozen_candidate_values_present": bool(
                q2093.get("verdict") == "KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN"
            ),
            "alpha_s_nonanchor_boundary_gate_present": bool(
                q2087.get("verdict") == "ALPHA_S_NONANCHOR_BOUNDARY_GATE_PASS"
            ),
            "t1_nonanchor_gate_present": bool(
                q2084.get("verdict") == "T1_NONANCHOR_STRICT_GATE_PASS"
            ),
            "gnewton_dimensionless_provenance_ready": bool(
                q2103.get("verdict") == "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PASS_STRICT_READY"
            ),
            "gnewton_direct_pack_ready": bool(
                q2113.get("verdict") == "GNEWTON_DIRECT_DIMENSIONLESS_PACK_READY"
            ),
            "a8_partial_gravity_bridge_present": bool(
                (a8.get("a8") or {}).get("verdict")
                == "strict-scope partial gravity bridge integrated; foundational gravity closure remains open"
            ),
        },
        "targets": targets,
        "status_counts": counts,
        "current_honest_reading": [
            "The repo already exports one strict-derived internal mass-ordering layer in dimensionless proxy scope.",
            "The electroweak and alpha_s entries do have strict internal promotion downstream in the repo, but each still carries explicit blockers that prevent overclaim.",
            "The gravity entry is present as a strict-side endpoint with explicit external observable origin, not as an internal derivation of G.",
        ],
        "hard_limits": [
            "This packet does not claim Standard Model identification.",
            "This packet does not claim proxy-to-GeV calibration is strict.",
            "This packet does not claim internal origin of G is closed.",
            "This packet does not claim full SM/GR closure.",
            "This packet does not claim ToE closure.",
        ],
        "recommended_next_move": "Build one narrow promotion audit that tests whether the QW-2093 boundary/formula layer for alpha_s and sin2 can be re-expressed through a kernel-invariants-only bridge registry without hidden legacy transfer or hidden regime constants.",
        "no_false_pass": True,
    }

    summary = {
        "stage": "F784",
        "packet_name": packet["packet_name"],
        "status": packet["status"],
        "as_of": AS_OF,
        "targets": {
            key: {
                "status_label": value["status_label"],
                "bridge_ready": value["bridge_ready"],
                "external_observable_origin": value["external_observable_origin"],
            }
            for key, value in targets.items()
        },
        "status_counts": counts,
        "recommended_next_move": packet["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(packet, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
