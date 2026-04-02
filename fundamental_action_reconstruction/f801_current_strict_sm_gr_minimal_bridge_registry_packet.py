#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
REPORTS_JSON = REPO / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json"

IN_P694 = GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
IN_F704 = GENERATED / "f704_current_strict_invariant_mass_observable_from_diagonal_local_psi_hessian_eigensystem_export_packet_summary.json"
IN_P704 = GENERATED / "p704_current_nonstrict_standard_model_host_matching_from_f704_h_psi_eigenvalue_proxy_probe_summary.json"
IN_P710 = GENERATED / "p710_current_nonstrict_proxy_to_gev_calibration_map_from_f704_eigenspectrum_probe_summary.json"

IN_T1_EW = REPO / "t1_nonanchor_observables_input_qw2085_2086.json"
IN_T1_ALPHA_S = REPO / "t1_nonanchor_alpha_s_input_qw2087.json"
IN_QW2086 = REPORTS_JSON / "report_qw2086_mz_nonanchor_ew_pole_gate.json"
IN_QW2087 = REPORTS_JSON / "report_qw2087_alpha_s_nonanchor_boundary_gate.json"
IN_QW2093 = REPORTS_JSON / "report_qw2093_kernel_derived_nonanchor_inputs_plan_executor.json"

IN_G_DIMLESS = REPO / "external_gnewton_bridge_qw2101.direct_dimensionless_ready.json"
IN_QW2101 = REPORTS_JSON / "report_qw2101_gnewton_bridge_external_autocollector.json"
IN_QW2103 = REPORTS_JSON / "report_qw2103_gnewton_dimensionless_provenance_gate.json"
IN_QW2113 = REPORTS_JSON / "report_qw2113_gnewton_direct_dimensionless_pack_gate.json"
IN_QW2207 = REPORTS_JSON / "report_qw2207_planck_internalization_obstruction_gate.json"
IN_A8 = GENERATED / "a8_gravity_bridge_summary.json"

OUT = GENERATED / "f801_current_strict_sm_gr_minimal_bridge_registry_packet.json"
OUT_SUMMARY = GENERATED / "f801_current_strict_sm_gr_minimal_bridge_registry_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_P694,
        IN_F704,
        IN_P704,
        IN_P710,
        IN_T1_EW,
        IN_T1_ALPHA_S,
        IN_QW2086,
        IN_QW2087,
        IN_QW2093,
        IN_G_DIMLESS,
        IN_QW2101,
        IN_QW2103,
        IN_QW2113,
        IN_QW2207,
        IN_A8,
    ]
    missing = [rel(p) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "F801",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p694 = load_json(IN_P694)
    f704 = load_json(IN_F704)
    p704 = load_json(IN_P704)
    p710 = load_json(IN_P710)
    t1_ew = load_json(IN_T1_EW)
    t1_alpha_s = load_json(IN_T1_ALPHA_S)
    qw2086 = load_json(IN_QW2086)
    qw2087 = load_json(IN_QW2087)
    qw2093 = load_json(IN_QW2093)
    g_dimless = load_json(IN_G_DIMLESS)
    qw2101 = load_json(IN_QW2101)
    qw2103 = load_json(IN_QW2103)
    qw2113 = load_json(IN_QW2113)
    qw2207 = load_json(IN_QW2207)
    a8 = load_json(IN_A8)

    entries: list[dict[str, Any]] = [
        {
            "id": "mass_ratio_ordering_layer",
            "status_label": "strict-derived",
            "scope": "strict_side_dimensionless_proxy_ordering_only",
            "value_summary": {
                "mass_proxy_m2_sorted_ascending": p694.get("mass_proxy_m2_sorted_ascending"),
                "mass_observable_object_ref": f704.get("exported_object_ref"),
            },
            "provenance_chain": [
                rel(IN_P694),
                rel(IN_F704),
            ],
            "external_observable_origin": None,
            "strict_admissibility": (
                "Uses exported strict projective mass-proxy computation and exported strict invariant mass observable only."
            ),
            "current_blockers": [
                "No proxy-to-GeV calibration inside strict scope.",
                "No Standard Model particle identification inside the minimal bridge.",
            ],
            "hard_limits": [
                "P704 host matching is excluded as non-strict.",
                "P710 proxy-to-GeV calibration is excluded as non-strict.",
            ],
        },
        {
            "id": "sin2_theta_w_eff",
            "status_label": "strict-derived",
            "scope": "strict_side_candidate_electroweak_observable",
            "value_summary": {
                "sin2_theta_w_eff": ((t1_ew.get("mz_pole_chain") or {}).get("sin2_theta_w_eff")),
                "origin": ((t1_ew.get("mz_pole_chain") or {}).get("sin2_theta_w_eff_origin")),
                "source": ((t1_ew.get("mz_pole_chain") or {}).get("source_sin2_theta_w_eff")),
                "qw2086_verdict": qw2086.get("verdict"),
                "qw2086_strict_nonanchor_pass": ((qw2086.get("flags") or {}).get("strict_nonanchor_pass")),
            },
            "provenance_chain": [
                rel(IN_T1_EW),
                rel(IN_QW2086),
                rel(IN_QW2093),
            ],
            "external_observable_origin": None,
            "strict_admissibility": (
                "Kernel-derived EW-pole chain is present and the non-anchor EW pole gate passes without retune or scan."
            ),
            "current_blockers": [
                "No explicit retained/replaced semantic-transfer verdict for the legacy Weinberg role.",
            ],
            "hard_limits": [
                "Do not reinterpret this as legacy sin^2(theta_W)=alpha_geo/12 transfer.",
                "This is a strict-side candidate observable, not a legacy-role discharge theorem.",
            ],
        },
        {
            "id": "alpha_s_boundary_mu0_alpha0",
            "status_label": "strict-derived",
            "scope": "strict_side_qcd_boundary_observable",
            "value_summary": {
                "mu0_gev": ((t1_alpha_s.get("alpha_s_boundary") or {}).get("mu0_gev")),
                "alpha_s_mu0": ((t1_alpha_s.get("alpha_s_boundary") or {}).get("alpha_s_mu0")),
                "origin": ((t1_alpha_s.get("alpha_s_boundary") or {}).get("origin")),
                "source": ((t1_alpha_s.get("alpha_s_boundary") or {}).get("source")),
                "qw2087_verdict": qw2087.get("verdict"),
                "qw2087_strict_nonanchor_pass": ((qw2087.get("flags") or {}).get("strict_nonanchor_pass")),
            },
            "provenance_chain": [
                rel(IN_T1_ALPHA_S),
                rel(IN_QW2087),
                rel(IN_QW2093),
            ],
            "external_observable_origin": None,
            "strict_admissibility": (
                "Kernel-derived alpha_s boundary and validation points are exported and pass the strict non-anchor gate."
            ),
            "current_blockers": [
                "Boundary source remains an explicit frozen hierarchy-log ansatz inside the strict pipeline.",
            ],
            "hard_limits": [
                "Do not upgrade the boundary object into a full QCD closure theorem.",
                "This entry packages the boundary observable only, not the entire radiative program.",
            ],
        },
        {
            "id": "g_dimensionless_mu_ref",
            "status_label": "strict-derived-with-external-observable-origin",
            "scope": "strict_side_gravity_bridge_observable_with_explicit_external_origin",
            "value_summary": {
                "mu_ref_gev": g_dimless.get("mu_ref_gev"),
                "g_dimensionless_mu_ref": g_dimless.get("g_dimensionless_mu_ref"),
                "bridge_observable_origin": g_dimless.get("bridge_observable_origin"),
                "qw2101_verdict": qw2101.get("verdict"),
                "qw2103_verdict": qw2103.get("verdict"),
                "qw2113_verdict": qw2113.get("verdict"),
                "qw2207_verdict": qw2207.get("verdict"),
                "a8_gravity_bridge_status": a8.get("status"),
            },
            "provenance_chain": [
                rel(IN_G_DIMLESS),
                rel(IN_QW2101),
                rel(IN_QW2103),
                rel(IN_QW2113),
                rel(IN_QW2207),
                rel(IN_A8),
            ],
            "external_observable_origin": "external_dimensionless_observable",
            "strict_admissibility": (
                "The dimensionless gravity bridge observable is explicit, provenance-clean, and integrated into strict gravity bridge scope."
            ),
            "current_blockers": [
                "Internal origin of the dimensionless G-bridge observable remains open.",
                "Full Einstein-Hilbert / GR foundational closure remains open.",
            ],
            "hard_limits": [
                "Do not promote this entry to full internal derivation of G.",
                "Do not use SI conversion or host-side constants as evidence of full strict gravity closure.",
            ],
        },
    ]

    excluded_interfaces = [
        {
            "id": "p704_nonstrict_host_matching",
            "status_label": "non-strict",
            "summary_ref": rel(IN_P704),
            "status": p704.get("status"),
            "reason": "Depends on external dataset plus host-matching policy.",
        },
        {
            "id": "p710_nonstrict_proxy_to_gev_calibration",
            "status_label": "non-strict",
            "summary_ref": rel(IN_P710),
            "status": p710.get("status"),
            "reason": "Proxy-to-GeV calibration is explicitly outside strict scope.",
        },
    ]

    payload = {
        "stage": "F801",
        "status": "PASS_EXPORTED_STRICT_SM_GR_MINIMAL_BRIDGE_REGISTRY_PACKET",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "Export one minimal strict-side bridge registry for four initial SM/GR-facing observables "
            "with explicit provenance and no false-pass promotion."
        ),
        "policy": {
            "working_note_ref": "WORKING_NOTE_LEGACY_KEEP_CUT_AND_MINIMAL_STRICT_SM_GR_BRIDGE.md",
            "allowed_status_labels": [
                "strict-derived",
                "strict-derived-with-external-observable-origin",
                "non-strict",
                "open",
            ],
            "bridge_stops_before_units": True,
            "legacy_role_transfer_disallowed": True,
            "no_false_pass": True,
        },
        "registry_entries": entries,
        "explicitly_excluded_interfaces": excluded_interfaces,
        "hard_limits": [
            "No Standard Model host matching claim.",
            "No strict proxy-to-GeV claim.",
            "No legacy Weinberg-role transfer claim.",
            "No legacy gravity-hierarchy transfer claim.",
            "No full internal origin claim for the dimensionless G bridge observable.",
            "No Einstein-Hilbert derivation claim.",
            "No full GR closure claim.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    counts: dict[str, int] = {}
    for entry in entries + excluded_interfaces:
        label = str(entry.get("status_label"))
        counts[label] = counts.get(label, 0) + 1

    summary = {
        "stage": "F801",
        "status": payload["status"],
        "as_of": AS_OF,
        "registry_ids": [entry["id"] for entry in entries],
        "status_counts": counts,
        "excluded_interface_ids": [item["id"] for item in excluded_interfaces],
        "no_false_pass": True,
    }

    write_json(OUT, payload)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
