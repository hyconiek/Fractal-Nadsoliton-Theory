#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
REPORTS_JSON = REPO / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json"

IN_F704 = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_P694 = GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
IN_P705 = GENERATED / "p705_current_actual_strict_projective_operational_toe_os_support_with_invariant_mass_observable_probe_summary.json"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_QW2063 = REPORTS_JSON / "report_qw2063_derivational_reconstruction_shared_flavor_basis.json"
IN_QW2070 = REPORTS_JSON / "report_qw2070_full_radiative_program_baseline.json"
IN_QW2087 = REPORTS_JSON / "report_qw2087_alpha_s_nonanchor_boundary_gate.json"
IN_QW2093 = REPORTS_JSON / "report_qw2093_kernel_derived_nonanchor_inputs_plan_executor.json"
IN_T1_ALPHA_S = REPO / "t1_nonanchor_alpha_s_input_qw2087.json"
IN_F787 = GENERATED / "f787_current_minimal_strict_sm_gr_bridge_export_refinement_packet.json"
IN_F801 = GENERATED / "f801_current_strict_sm_gr_minimal_bridge_registry_packet.json"

OUT = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
OUT_SUMMARY = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return (
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("_", " ")
        .replace("-", " ")
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def check_flag(summary: dict[str, Any], check_id: str) -> bool:
    for item in summary.get("checks") or []:
        if item.get("id") == check_id:
            return bool(item.get("pass"))
    return False


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_F704,
        IN_P694,
        IN_P705,
        IN_N703,
        IN_QW2063,
        IN_QW2070,
        IN_QW2087,
        IN_QW2093,
        IN_T1_ALPHA_S,
        IN_F787,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P788",
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
    p694 = load_json(IN_P694)
    p705 = load_json(IN_P705)
    n703_text = IN_N703.read_text(encoding="utf-8")
    qw2063 = load_json(IN_QW2063)
    qw2070 = load_json(IN_QW2070)
    qw2087 = load_json(IN_QW2087)
    qw2093 = load_json(IN_QW2093)
    t1_alpha_s = load_json(IN_T1_ALPHA_S)
    f787 = load_json(IN_F787)
    f801 = load_json(IN_F801) if IN_F801.exists() else None

    alpha_s_boundary = t1_alpha_s.get("alpha_s_boundary") or {}
    validation_points = t1_alpha_s.get("validation_points") or []
    qcd_baseline = qw2070.get("qcd_one_loop_baseline") or {}
    q2093_formulas = qw2093.get("frozen_formulas") or {}

    adapter_search_files = [
        IN_T1_ALPHA_S,
        IN_QW2070,
        IN_QW2087,
        IN_QW2093,
        IN_F787,
    ]
    if f801 is not None:
        adapter_search_files.append(IN_F801)

    adapter_tokens = [
        "mu0_dimensionless",
        "mu0_normalized",
        "alpha_s_boundary_ratio",
        "validation_points_normalized",
        "normalized_alpha_s_boundary",
        "dimensionless_alpha_s_boundary",
    ]
    adapter_hits = {}
    for path in adapter_search_files:
        text = path.read_text(encoding="utf-8")
        hits = [token for token in adapter_tokens if token in text]
        if hits:
            adapter_hits[rel(path)] = hits

    checks = [
        {
            "id": "strict_mass_proxy_layer_exists",
            "pass": (
                f704.get("status") == "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT"
                and check_flag(p705, "basis_invariant_mass_observable_exported")
                and contains_all(n703_text, ["dimensionless quadratic coefficients", "not yet physical masses in GeV"])
            ),
            "details": "F704/P705/N703 together must export a strict dimensionless mass proxy layer.",
        },
        {
            "id": "current_alpha_s_input_is_gev_interface",
            "pass": (
                alpha_s_boundary.get("mu0_gev") is not None
                and any(point.get("mu_gev") is not None for point in validation_points)
            ),
            "details": "Current alpha_s input schema uses GeV-level boundary and validation scales.",
        },
        {
            "id": "q2093_alpha_s_formula_uses_mass_chain_symbols",
            "pass": (
                q2093_formulas.get("alpha_s_boundary_mu0") == "m_bottom"
                and "m_top/m_bottom" in str(q2093_formulas.get("alpha_s_boundary_alpha0"))
            ),
            "details": "QW-2093 alpha_s formula still routes through mass-chain symbols rather than a dimensionless strict proxy boundary.",
        },
        {
            "id": "q2063_foundational_constants_not_strict_first_principles",
            "pass": ((qw2063.get("flags") or {}).get("strict_first_principles_foundational_constants_derived") is False),
            "details": "QW-2063 still keeps foundational constants below strict first-principles promotion.",
        },
        {
            "id": "qcd_baseline_is_gev_only_interface",
            "pass": (
                qcd_baseline.get("mu0_gev") is not None
                and any(
                    isinstance(sample, dict) and sample.get("mu_gev") is not None
                    for sample in (qcd_baseline.get("samples") or [])
                )
            ),
            "details": "The current QCD running baseline accepts GeV-scale inputs only.",
        },
        {
            "id": "no_exported_dimensionless_alpha_s_adapter_detected",
            "pass": len(adapter_hits) == 0,
            "details": "No exported file currently advertises a dimensionless/normalized alpha_s boundary adapter interface.",
        },
        {
            "id": "f787_already_demotes_alpha_s_from_minimal_export_set",
            "pass": "alpha_s_boundary_mu0_alpha0" in ((f787.get("export_sets") or {}).get("support_only_nonexport_ids") or []),
            "details": "F787 should already keep alpha_s outside the minimal export set on current repo state.",
        },
    ]

    failed = [item["id"] for item in checks if not item["pass"]]

    route_ready = (
        checks[0]["pass"]
        and not checks[1]["pass"]
        and not checks[2]["pass"]
        and not checks[3]["pass"]
        and not checks[4]["pass"]
        and not checks[5]["pass"]
    )

    if route_ready:
        status = "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_PRESENT"
    else:
        status = "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE"

    artifact = {
        "stage": "P788",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f704_mass_observable": rel(IN_F704),
            "p694_mass_proxy_probe": rel(IN_P694),
            "p705_os_support_v3_summary": rel(IN_P705),
            "n703_mass_proxy_meaning": rel(IN_N703),
            "qw2063_mass_chain": rel(IN_QW2063),
            "qw2070_qcd_baseline": rel(IN_QW2070),
            "qw2087_alpha_s_gate": rel(IN_QW2087),
            "qw2093_frozen_formula_set": rel(IN_QW2093),
            "t1_alpha_s_input": rel(IN_T1_ALPHA_S),
            "f787_export_refinement": rel(IN_F787),
            "f801_minimal_bridge_registry_optional": rel(IN_F801) if f801 is not None else None,
        },
        "existing_objects": {
            "strict_mass_proxy_layer": {
                "f704_status": f704.get("status"),
                "p694_mass_proxy_m2_sorted_ascending": p694.get("mass_proxy_m2_sorted_ascending"),
                "n703_meaning_boundary_present": contains_all(
                    n703_text,
                    ["dimensionless quadratic coefficients", "not yet physical masses in GeV"],
                ),
            },
            "current_alpha_s_boundary_lane": {
                "t1_alpha_s_boundary": alpha_s_boundary,
                "qw2087_boundary_source": ((qw2087.get("checks") or {}).get("boundary_source")),
                "qw2093_frozen_formula_mu0": q2093_formulas.get("alpha_s_boundary_mu0"),
                "qw2093_frozen_formula_alpha0": q2093_formulas.get("alpha_s_boundary_alpha0"),
                "qcd_baseline_mu0_gev": qcd_baseline.get("mu0_gev"),
            },
            "optional_overlap_note": {
                "f801_present": f801 is not None,
                "f801_alpha_s_scope": None if f801 is None else ((f801.get("registry_entries") or [])[2].get("scope") if len(f801.get("registry_entries") or []) > 2 else None),
            },
        },
        "checks": checks,
        "failed_checks": failed,
        "adapter_token_hits": adapter_hits,
        "missing_interfaces": [
            "A dimensionless or explicitly normalized alpha_s boundary-scale object replacing mu0_gev.",
            "A normalized validation-point grid replacing mu_gev samples for the alpha_s lane.",
            "An exported adapter from strict mass proxies / normalized strict observables into the alpha_s boundary schema.",
            "A QCD-running interface that can consume that normalized boundary without silently reintroducing GeV-level calibration.",
        ],
        "current_honest_reading": [
            "F704/N703/P705 do export a real strict dimensionless mass-proxy layer, but that layer is still only a proxy/order object and not an alpha_s boundary adapter.",
            "The active alpha_s lane still starts from mu0_gev plus mu_gev validation points, so the current interface is GeV-level by construction.",
            "QW-2093 still defines alpha_s through m_bottom and ln(m_top/m_bottom)+delta_eta, so the formula layer is not yet a dimensionless or explicitly normalized strict replacement route.",
            "QW-2063 still reports strict_first_principles_foundational_constants_derived=false, which blocks silent promotion of the mass chain into a strict boundary interface.",
            "The current repo therefore does not yet export a real dimensionless or explicitly normalized alpha_s boundary replacement route after F787.",
        ],
        "recommended_next_packet": {
            "id": "F788_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_TARGET_PACKET",
            "goal": "Export one explicit normalized alpha_s boundary interface object with mu0_tilde, alpha_s_mu0, normalized_validation_points, and normalization_rule_ref, without claiming QCD closure or GeV identification.",
            "minimum_fields": [
                "mu0_tilde",
                "alpha_s_mu0",
                "n_f_active_at_mu0",
                "normalized_validation_points",
                "normalization_rule_ref",
                "strict_input_chain",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P788",
        "status": status,
        "as_of": AS_OF,
        "replacement_route_ready": route_ready,
        "missing_interfaces": artifact["missing_interfaces"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
