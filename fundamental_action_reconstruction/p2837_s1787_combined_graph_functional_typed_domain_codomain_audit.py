#!/usr/bin/env python3
"""P2837/S1787: typed domain/codomain obligation audit.

P2836 closed the units/normalization attempt as a no-go.  P2837 attacks exactly
one remaining P2835/P2836 theorem obligation: whether the combined graph
functional has an exported typed map into K, L_total, a source density, or a
unit-bearing coefficient.  The audit is intentionally a typed obstruction table,
not another finite-separation refinement.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2835 = GEN / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.json"
P2836 = GEN / "p2836_s1786_combined_graph_functional_units_normalization_obligation_audit.json"
OUT = GEN / "p2837_s1787_combined_graph_functional_typed_domain_codomain_audit.json"
MD = GEN / "p2837_s1787_combined_graph_functional_typed_domain_codomain_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

EXPECTED_FULL_CARRIER_COUNT = 16828

REQUIRED_TYPED_PREMISES = (
    "source_domain_object",
    "target_codomain_object",
    "evaluation_or_pullback_map",
    "target_independent_units_or_scale",
    "locality_or_covariance_rule",
    "variational_variable_identification",
    "coupling_coefficient_rule",
)


def candidate_typed_maps() -> list[dict[str, Any]]:
    candidates = [
        {
            "candidate": "graph_to_kernel_scalar_control",
            "proposed_type": "F_combined: Graph_16,4 -> parameter/control for K(d)",
            "premises": {
                "source_domain_object": True,
                "target_codomain_object": True,
                "evaluation_or_pullback_map": False,
                "target_independent_units_or_scale": False,
                "locality_or_covariance_rule": False,
                "variational_variable_identification": False,
                "coupling_coefficient_rule": False,
            },
            "obstruction": "No exported map turns a finite graph invariant into a function of distance d or into a K(d) parameter with units.",
        },
        {
            "candidate": "graph_to_lagrangian_density_term",
            "proposed_type": "F_combined: Graph_16,4 -> local density contribution L_graph(x)",
            "premises": {
                "source_domain_object": True,
                "target_codomain_object": True,
                "evaluation_or_pullback_map": False,
                "target_independent_units_or_scale": False,
                "locality_or_covariance_rule": False,
                "variational_variable_identification": False,
                "coupling_coefficient_rule": False,
            },
            "obstruction": "No exported pullback embeds the finite graph carrier into spacetime/local field variables x, fields, or metric data.",
        },
        {
            "candidate": "graph_to_source_density",
            "proposed_type": "F_combined: Graph_16,4 -> J(x) or source density",
            "premises": {
                "source_domain_object": True,
                "target_codomain_object": True,
                "evaluation_or_pullback_map": False,
                "target_independent_units_or_scale": False,
                "locality_or_covariance_rule": False,
                "variational_variable_identification": False,
                "coupling_coefficient_rule": False,
            },
            "obstruction": "The carrier has no exported localization map assigning graph rows to x-dependent source support.",
        },
        {
            "candidate": "graph_to_dimensionless_coefficient",
            "proposed_type": "F_combined: Graph_16,4 -> dimensionless coefficient label",
            "premises": {
                "source_domain_object": True,
                "target_codomain_object": True,
                "evaluation_or_pullback_map": True,
                "target_independent_units_or_scale": False,
                "locality_or_covariance_rule": False,
                "variational_variable_identification": False,
                "coupling_coefficient_rule": False,
            },
            "obstruction": "A dimensionless coefficient label is typable as finite data, but P2836 leaves its scale free and no coupling law attaches it to K or L_total.",
        },
    ]
    for row in candidates:
        row["satisfied_premise_count"] = sum(bool(row["premises"][name]) for name in REQUIRED_TYPED_PREMISES)
        row["missing_premises"] = [name for name in REQUIRED_TYPED_PREMISES if not row["premises"][name]]
        row["accepted_as_typed_K_or_Ltotal_map"] = not row["missing_premises"]
    return candidates


def typed_obstruction_summary(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [row["candidate"] for row in candidates if row["accepted_as_typed_K_or_Ltotal_map"]]
    blocker_histogram: dict[str, int] = {name: 0 for name in REQUIRED_TYPED_PREMISES}
    for row in candidates:
        for missing in row["missing_premises"]:
            blocker_histogram[missing] += 1
    return {
        "candidate_count": len(candidates),
        "accepted_typed_map_count": len(accepted),
        "accepted_typed_maps": accepted,
        "blocker_histogram": blocker_histogram,
        "common_hard_blockers": [name for name, count in blocker_histogram.items() if count == len(candidates)],
    }


def build_audit(p2835: dict[str, Any], p2836: dict[str, Any]) -> dict[str, Any]:
    p2835_sep = p2835["combined_witness_source_law_theorem_obligation_audit"]["combined_separator"]
    candidates = candidate_typed_maps()
    summary = typed_obstruction_summary(candidates)
    return {
        "input_statuses_rechecked": {"P2835": p2835.get("status"), "P2836": p2836.get("status")},
        "combined_separator_rechecked": {
            "combined_class_count": p2835_sep["combined_class_count"],
            "combined_collision_class_count": p2835_sep["combined_collision_class_count"],
            "p2834_patch_graph_count": p2835_sep["p2834_patch_graph_count"],
        },
        "required_typed_premises": list(REQUIRED_TYPED_PREMISES),
        "candidate_typed_maps": candidates,
        "typed_obstruction_summary": summary,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    sep = audit["combined_separator_rechecked"]
    summary = audit["typed_obstruction_summary"]
    facts = {
        "p2835_combined_separator_rechecked": sep["combined_class_count"] == EXPECTED_FULL_CARRIER_COUNT and sep["combined_collision_class_count"] == 0,
        "at_least_one_candidate_typed_map_tested": summary["candidate_count"] > 0,
        "accepted_typed_K_or_Ltotal_map_exported": summary["accepted_typed_map_count"] > 0,
        "target_independent_units_available_from_p2836": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted_as_typed_map = all([
        facts["p2835_combined_separator_rechecked"],
        facts["at_least_one_candidate_typed_map_tested"],
        facts["accepted_typed_K_or_Ltotal_map_exported"],
        facts["target_independent_units_available_from_p2836"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_typed_domain_codomain_audit": facts["p2835_combined_separator_rechecked"] and facts["at_least_one_candidate_typed_map_tested"],
        "accepted_as_typed_K_or_Ltotal_map": accepted_as_typed_map,
        "accepted_as_typed_domain_codomain_no_go": not accepted_as_typed_map,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["combined_graph_functional_typed_domain_codomain_audit"]
    summary = audit["typed_obstruction_summary"]
    lines = [
        "# P2837/S1787 combined graph functional typed domain/codomain audit", "", f"Status: `{payload['status']}`", "",
        "## Rechecked finite separator",
        f"- combined_class_count={audit['combined_separator_rechecked']['combined_class_count']}",
        f"- combined_collision_class_count={audit['combined_separator_rechecked']['combined_collision_class_count']}", "",
        "## Typed obstruction result",
        f"- candidate_count={summary['candidate_count']}",
        f"- accepted_typed_map_count={summary['accepted_typed_map_count']}",
        f"- common_hard_blockers={summary['common_hard_blockers']}", "",
        "## Acceptance",
        f"- accepted_as_typed_domain_codomain_audit={payload['acceptance_matrix']['accepted_as_typed_domain_codomain_audit']}",
        f"- accepted_as_typed_K_or_Ltotal_map={payload['acceptance_matrix']['accepted_as_typed_K_or_Ltotal_map']}",
        f"- accepted_as_typed_domain_codomain_no_go={payload['acceptance_matrix']['accepted_as_typed_domain_codomain_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2835 = read_json(P2835)
    p2836 = read_json(P2836)
    audit = build_audit(p2835, p2836)
    payload: dict[str, Any] = {
        "status": "P2837_TYPED_DOMAIN_CODOMAIN_OBLIGATION_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2835": sha(P2835), "P2836": sha(P2836), "16_4_4.scd": sha(SCD)},
        "combined_graph_functional_typed_domain_codomain_audit": audit,
        "decision": {
            "negative_export_flags": {
                "typed_K_or_Ltotal_map_exported": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "proved_variational_derivative_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2837 attacks exactly one remaining theorem obligation: typed domain/codomain into K or L_total.  Four candidate typed maps are audited.  None satisfies all required premises.  The strongest finite option is only a dimensionless coefficient label; P2836 leaves its scale free, and no current artifact exports an evaluation/pullback, locality/covariance rule, variational variable identification, or unit-bearing coupling law into K or L_total.",
            "next_honest_step": "Do not replay graph separation, units averaging, or generic typed-map names.  The next admissible proof-grade move should attack exactly one remaining premise: a formal variational derivative theorem for the finite graph functional, or an explicit evaluation/pullback/localization map into field variables.  If neither is exported, preserve the P2831-P2837 finite-witness/no-units/no-typed-coupling boundary and pivot away from graph-source promotion.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2837/S1787 combined graph functional typed domain codomain audit", "## P2837/S1787 combined graph functional typed domain codomain audit\n\n`P2837/S1787` attacks exactly one P2835/P2836 missing theorem obligation: a typed domain/codomain map into `K` or `L_total`.  Four candidate maps are audited (`graph -> K(d)` control, `graph -> L_graph(x)` density, `graph -> J(x)` source density, and `graph -> dimensionless coefficient`).  None satisfies all premises: the finite graph domain is available, but current artifacts do not export the needed evaluation/pullback/localization map, target-independent units, locality/covariance rule, variational variable identification, or unit-bearing coupling coefficient.  No strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2837/S1787 typed domain codomain Ltotal guard", "## P2837/S1787 typed domain codomain Ltotal guard\n\n`P2837/S1787` adds no term to `L_total`.  Candidate maps from the combined graph functional into `K`, a local density, a source density, or a coefficient remain blocked by missing pullback/localization, units, locality/covariance, variational-variable, and coupling premises.\n")
    append_once(AGENTS, "Current typed domain-codomain obligation guardrail (P2837/S1787, 2026-06-17)", "## Current typed domain-codomain obligation guardrail (P2837/S1787, 2026-06-17)\n\n- P2837 attacks one remaining P2835/P2836 theorem obligation: a typed domain/codomain map from the combined P2833/P2834 graph functional into `K` or `L_total`.\n- Four candidate maps are audited and all fail at least one required premise; current artifacts do not export the needed evaluation/pullback/localization, target-independent units, locality/covariance, variational-variable identification, or unit-bearing coupling rule.\n- Do not promote P2837 to a strict graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must attack exactly one remaining premise, such as a formal variational derivative theorem or explicit localization map, or preserve the no-typed-coupling boundary.\n")
    return payload


if __name__ == "__main__":
    main()
