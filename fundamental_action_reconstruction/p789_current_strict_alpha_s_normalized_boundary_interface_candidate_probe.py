#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F704 = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_P694 = GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_F789 = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet.json"

OUT = GENERATED / "p789_current_strict_alpha_s_normalized_boundary_interface_candidate_probe.json"
OUT_SUMMARY = GENERATED / "p789_current_strict_alpha_s_normalized_boundary_interface_candidate_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def round_list(values: list[float], digits: int = 12) -> list[float]:
    return [round(value, digits) for value in values]


def geometric_mean(values: list[float]) -> float:
    return math.exp(sum(math.log(v) for v in values) / len(values))


def family_checks(*, canonical_anchor_exported: bool, normalization_rule_exported: bool, nf_mapping_exported: bool, normalized_running_consumer_exported: bool) -> dict[str, bool]:
    return {
        "dimensionless_points_present": True,
        "derived_from_exported_strict_objects": True,
        "no_gev_fields_used": True,
        "canonical_alpha_s_boundary_anchor_exported": canonical_anchor_exported,
        "normalization_rule_exported": normalization_rule_exported,
        "nf_mapping_exported": nf_mapping_exported,
        "normalized_running_consumer_exported": normalized_running_consumer_exported,
    }


def family_ready(checks: dict[str, bool]) -> bool:
    return all(checks.values())


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F704, IN_P694, IN_N703, IN_F789]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P789",
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
    f789 = load_json(IN_F789)

    m_proxy = [float(v) for v in ((f704.get("outputs") or {}).get("m_proxy_ascending") or [])]
    pair_m2_sorted = (p694.get("mass_proxy_m2_sorted_ascending") or [])
    pair_m = [math.sqrt(float(item["m2"])) for item in pair_m2_sorted]

    max_mode = max(m_proxy)
    gm_mode = geometric_mean(m_proxy)
    max_pair = max(pair_m)

    families = {
        "f704_max_mode_anchor_family": {
            "source_object": "MassObservable_diagonal_local_strict_derived_v1",
            "normalization_rule": "mu_tilde_i := m_proxy_i / max_j(m_proxy_j)",
            "mu0_tilde_candidate": 1.0,
            "normalized_validation_points_candidate": round_list([value / max_mode for value in m_proxy]),
            "normalization_rule_ref_candidate": "probe_local_rule:f704_max_mode_anchor_family",
            "strict_input_chain": [
                rel(IN_F704),
                rel(IN_N703),
            ],
            "checks": family_checks(
                canonical_anchor_exported=False,
                normalization_rule_exported=False,
                nf_mapping_exported=False,
                normalized_running_consumer_exported=False,
            ),
            "why_not_export_ready": [
                "F704 exports a dimensionless spectrum, but no current source upgrades the max mode into a canonical alpha_s boundary anchor.",
                "The normalization rule exists only probe-locally; no exported authority object currently names it.",
                "No strict-side map exports n_f_active_at_mu0 from this normalized boundary.",
                "No exported QCD-running interface consumes this normalized grid directly.",
            ],
        },
        "f704_geometric_mean_anchor_family": {
            "source_object": "MassObservable_diagonal_local_strict_derived_v1",
            "normalization_rule": "mu_tilde_i := m_proxy_i / geometric_mean(m_proxy)",
            "mu0_tilde_candidate": 1.0,
            "normalized_validation_points_candidate": round_list([value / gm_mode for value in m_proxy]),
            "normalization_rule_ref_candidate": "probe_local_rule:f704_geometric_mean_anchor_family",
            "strict_input_chain": [
                rel(IN_F704),
                rel(IN_N703),
            ],
            "checks": family_checks(
                canonical_anchor_exported=False,
                normalization_rule_exported=False,
                nf_mapping_exported=False,
                normalized_running_consumer_exported=False,
            ),
            "why_not_export_ready": [
                "The geometric-mean anchor is dimensionless and internal, but no current source exports it as the canonical alpha_s boundary anchor.",
                "This family introduces a useful normalized scale, yet still no strict-side n_f map or normalized running consumer is exported.",
            ],
        },
        "p694_pairmax_anchor_family": {
            "source_object": "P694 pair-rayleigh mass proxy summary",
            "normalization_rule": "mu_tilde_pair := m_pair_proxy / max_pair_proxy",
            "mu0_tilde_candidate": 1.0,
            "normalized_validation_points_candidate": round_list([value / max_pair for value in pair_m]),
            "normalization_rule_ref_candidate": "probe_local_rule:p694_pairmax_anchor_family",
            "strict_input_chain": [
                rel(IN_P694),
                rel(IN_N703),
            ],
            "checks": family_checks(
                canonical_anchor_exported=False,
                normalization_rule_exported=False,
                nf_mapping_exported=False,
                normalized_running_consumer_exported=False,
            ),
            "why_not_export_ready": [
                "This family is weaker than F704 because it depends on pair-rayleigh proxies rather than the full basis-invariant spectrum.",
                "No current source promotes a pair-max anchor into the canonical alpha_s boundary slot.",
                "No strict-side n_f map or normalized running consumer is exported here either.",
            ],
        },
    }

    for family in families.values():
        family["candidate_export_ready"] = family_ready(family["checks"])

    strongest_family_id = "f704_max_mode_anchor_family"
    candidate_family_constructible = True
    canonical_export_ready = any(family["candidate_export_ready"] for family in families.values())

    if candidate_family_constructible and not canonical_export_ready:
        status = "P789_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_CANDIDATE_FAMILY_PRESENT_CANONICAL_EXPORT_BLOCKED"
    elif canonical_export_ready:
        status = "P789_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_CANDIDATE_EXPORT_READY"
    else:
        status = "P789_REQUIRES_REVIEW"

    residual_blockers = [
        "No exported canonical alpha_s boundary anchor is selected from the strict dimensionless candidate families.",
        "No exported normalization_rule_ref currently upgrades any candidate family into an authoritative interface object.",
        "No exported n_f_active mapping is attached to the normalized candidate families.",
        "No exported normalized QCD-running consumer accepts mu_tilde-based validation grids.",
    ]

    artifact = {
        "stage": "P789",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f704_mass_observable": rel(IN_F704),
            "p694_mass_proxy_probe": rel(IN_P694),
            "n703_meaning_theorem": rel(IN_N703),
            "f789_interface_target": rel(IN_F789),
        },
        "candidate_family_constructible": candidate_family_constructible,
        "canonical_export_ready": canonical_export_ready,
        "strongest_family_id": strongest_family_id,
        "candidate_families": families,
        "current_honest_reading": [
            "A dimensionless candidate family for mu0_tilde plus normalized validation points can already be assembled from exported strict objects, strongest on the F704 basis-invariant mass spectrum.",
            "That is still weaker than a canonical interface export, because the repo does not yet choose one boundary anchor, one normalization rule, one n_f map, or one normalized QCD-running consumer.",
            "So the blocker has narrowed from 'no route at all' to 'candidate families exist but canonical export remains blocked'.",
        ],
        "residual_blockers_after_p789": residual_blockers,
        "recommended_next_move": "Export one narrow canonical-anchor probe that tests whether the strongest family f704_max_mode_anchor_family can be upgraded into an authoritative alpha_s boundary anchor without host matching or GeV reintroduction.",
        "no_false_pass": True,
    }

    summary = {
        "stage": "P789",
        "status": status,
        "as_of": AS_OF,
        "candidate_family_constructible": candidate_family_constructible,
        "canonical_export_ready": canonical_export_ready,
        "strongest_family_id": strongest_family_id,
        "residual_blocker_count": len(residual_blockers),
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
