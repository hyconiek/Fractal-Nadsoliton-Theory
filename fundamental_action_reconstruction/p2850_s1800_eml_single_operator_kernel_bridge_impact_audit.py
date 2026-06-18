#!/usr/bin/env python3
"""P2850/S1800: EML single-operator paper impact audit for kernel bridge.

The user supplied arXiv:2603.21852v2 ("All elementary functions from a single
operator").  P2850 checks whether that information changes the current
P2849 damping/compression bridge boundary.

The paper's relevant claim is representational: the binary EML operator
EML(x, y) = exp(x) - log(y), with terminal 1, can generate the ordinary
scientific-calculator elementary basis.  Both legacy and strict kernels are
already elementary expressions, so EML may provide a uniform syntax/circuit
basis for writing them.  This audit tests whether such syntactic universality
also supplies any missing typed bridge/source premise.  It does not.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2849 = GEN / "p2849_s1799_damping_compression_kernel_bridge_atom_audit.json"
OUT = GEN / "p2850_s1800_eml_single_operator_kernel_bridge_impact_audit.json"
MD = GEN / "p2850_s1800_eml_single_operator_kernel_bridge_impact_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

LEGACY_ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY_OMEGA = math.pi / 4.0
LEGACY_PHI = math.pi / 6.0
LEGACY_BETA_TORS = 0.01
STRICT_OMEGA = 0.18575
STRICT_PHI = 0.16250
STRICT_BETA = 1.0
STRICT_ETA = 1.8
DISTANCES = tuple(range(1, 13))
ARXIV_HTML = "https://arxiv.org/html/2603.21852v2"
ARXIV_PDF = "https://arxiv.org/pdf/2603.21852"

REQUIRED_IMPACT_PREMISES = (
    "external_paper_checked",
    "eml_represents_elementary_syntax",
    "kernel_formulas_are_elementary",
    "eml_exports_parameter_source_law",
    "eml_exports_beta_eta_source",
    "eml_exports_amplitude_source",
    "eml_exports_phase_topological_selector",
    "eml_exports_completion_map_semantics",
    "eml_exports_role_transfer_theorem",
)


def eml(x: float, y: float) -> float:
    return math.exp(x) - math.log(y)


def legacy_kernel(d: int) -> float:
    return LEGACY_ALPHA_GEO * math.cos(LEGACY_OMEGA * d + LEGACY_PHI) / (1.0 + LEGACY_BETA_TORS * d)


def strict_kernel(d: int) -> float:
    return math.cos(STRICT_OMEGA * d + STRICT_PHI) / (1.0 + STRICT_BETA * (d ** STRICT_ETA))


def eml_sanity_rows() -> list[dict[str, Any]]:
    rows = []
    for x in (0.0, 0.5, 1.0, 2.0):
        via_eml = eml(x, 1.0)
        direct = math.exp(x)
        rows.append(
            {
                "x": x,
                "eml_x_1": via_eml,
                "exp_x": direct,
                "abs_error": abs(via_eml - direct),
            }
        )
    return rows


def kernel_value_rows() -> list[dict[str, Any]]:
    rows = []
    for d in DISTANCES:
        legacy = legacy_kernel(d)
        strict = strict_kernel(d)
        rows.append(
            {
                "d": d,
                "legacy_kernel": legacy,
                "strict_kernel": strict,
                "difference_strict_minus_legacy": strict - legacy,
                "absolute_difference": abs(strict - legacy),
            }
        )
    return rows


def elementary_syntax_classification() -> dict[str, Any]:
    components = {
        "legacy_kernel_components": ["multiplication", "cos", "addition", "division", "linear_polynomial"],
        "strict_kernel_components": ["cos", "addition", "division", "power", "nonlinear_polynomial_power"],
        "eml_paper_relevant_claim": "EML(x,y)=exp(x)-log(y) with terminal 1 generates ordinary elementary scientific-calculator functions; therefore elementary kernel formulas can in principle be re-expressed in that syntax.",
    }
    return {
        **components,
        "legacy_formula_elementary": True,
        "strict_formula_elementary": True,
        "eml_representation_relevance": "syntax_only_unless_typed_source_law_added",
    }


def premise_matrix() -> dict[str, Any]:
    premises = {
        "external_paper_checked": True,
        "eml_represents_elementary_syntax": True,
        "kernel_formulas_are_elementary": True,
        "eml_exports_parameter_source_law": False,
        "eml_exports_beta_eta_source": False,
        "eml_exports_amplitude_source": False,
        "eml_exports_phase_topological_selector": False,
        "eml_exports_completion_map_semantics": False,
        "eml_exports_role_transfer_theorem": False,
    }
    return {
        "premises": premises,
        "missing_premises": [key for key in REQUIRED_IMPACT_PREMISES if not premises[key]],
        "changes_p2849_bridge_boundary": False,
        "accepted_as_external_information_impact_audit": True,
    }


def build_payload(p2849: dict[str, Any]) -> dict[str, Any]:
    eml_rows = eml_sanity_rows()
    kernel_rows = kernel_value_rows()
    max_exp_error = max(row["abs_error"] for row in eml_rows)
    max_kernel_gap = max(row["absolute_difference"] for row in kernel_rows)
    matrix = premise_matrix()
    facts = {
        "p2849_rechecked": p2849.get("status") == "P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE",
        "eml_exp_identity_sanity_passes": max_exp_error < 1e-15,
        "kernel_formulas_classified_elementary": True,
        "kernel_values_not_identified_by_syntax": max_kernel_gap > 0.0,
        "no_parameter_source_law_exported": not matrix["premises"]["eml_exports_parameter_source_law"],
        "no_beta_eta_source_exported": not matrix["premises"]["eml_exports_beta_eta_source"],
        "p2849_boundary_unchanged": not matrix["changes_p2849_bridge_boundary"],
    }
    return {
        "status": "P2850_EML_SINGLE_OPERATOR_KERNEL_BRIDGE_IMPACT_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2849": sha(P2849)},
        "external_sources_checked": {
            "arxiv_html": ARXIV_HTML,
            "arxiv_pdf": ARXIV_PDF,
            "paper_title": "All elementary functions from a single operator",
            "arxiv_id_version": "2603.21852v2",
            "relevant_claim_summary": "EML(x,y)=Exp[x]-Log[y], with terminal 1, is claimed to reconstruct the elementary scientific-calculator basis; the paper frames this as syntactic/operator-basis universality and separates heuristic search from independent verification.",
        },
        "eml_impact_audit": {
            "input_statuses_rechecked": {"P2849": p2849.get("status")},
            "eml_sanity_rows": eml_rows,
            "max_exp_identity_abs_error": max_exp_error,
            "elementary_syntax_classification": elementary_syntax_classification(),
            "kernel_value_rows": kernel_rows,
            "max_legacy_strict_kernel_abs_gap_on_audited_distances": max_kernel_gap,
            "required_impact_premises": list(REQUIRED_IMPACT_PREMISES),
            "premise_matrix": matrix,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_eml_external_information_impact_audit": all(facts.values()),
            "exports_new_kernel_bridge_source": False,
        },
        "decision": {
            "negative_export_flags": {
                "p2849_boundary_changed": False,
                "kernel_bridge_exported": False,
                "beta_eta_source_exported": False,
                "amplitude_source_exported": False,
                "selector_source_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "The arXiv EML result changes the available expression syntax, not the typed source obligations.  It supports treating legacy and strict kernel formulas as elementary expressions that may be encoded in one operator basis, but it does not export beta/eta, amplitude, phase/topological selector, completion-map semantics, unit-bearing L_total coupling, or role-transfer source laws.  P2849's damping/compression obstruction is unchanged.",
            "next_honest_step": "Use the EML paper only as an optional syntax/compression lens.  The next proof-grade move remains one typed bridge-source atom: either a new strict source law for eta and target-independent beta, or a separate amplitude-normalization passage theorem.  If no such typed source premise is supplied, preserve the no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["eml_impact_audit"]
    lines = [
        "# P2850/S1800 EML single-operator kernel bridge impact audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## External sources checked",
        f"- HTML: {payload['external_sources_checked']['arxiv_html']}",
        f"- PDF: {payload['external_sources_checked']['arxiv_pdf']}",
        f"- relevant_claim_summary={payload['external_sources_checked']['relevant_claim_summary']}",
        "",
        "## EML sanity check",
        f"- max_exp_identity_abs_error={audit['max_exp_identity_abs_error']}",
        "",
        "## Elementary syntax classification",
        f"- legacy_formula_elementary={audit['elementary_syntax_classification']['legacy_formula_elementary']}",
        f"- strict_formula_elementary={audit['elementary_syntax_classification']['strict_formula_elementary']}",
        f"- eml_representation_relevance={audit['elementary_syntax_classification']['eml_representation_relevance']}",
        "",
        "## Kernel value separation",
        f"- max_legacy_strict_kernel_abs_gap_on_audited_distances={audit['max_legacy_strict_kernel_abs_gap_on_audited_distances']}",
        "",
        "## Premise matrix",
        f"- changes_p2849_bridge_boundary={audit['premise_matrix']['changes_p2849_bridge_boundary']}",
        f"- missing_premises={audit['premise_matrix']['missing_premises']}",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2849))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2850/S1800 EML single-operator kernel bridge impact audit", "## P2850/S1800 EML single-operator kernel bridge impact audit\n\n`P2850/S1800` checks the supplied arXiv `2603.21852v2` EML single-operator result against the current kernel-bridge boundary.  The result is relevant as syntax/operator-basis compression: `EML(x,y)=exp(x)-log(y)` with terminal `1` can represent elementary scientific-calculator functions, and both `K_legacy_ont` and `K_strict_gate` are elementary formulas.  However, this does not export a typed source law for `beta/eta`, amplitude, phase/topological selector data, completion-map semantics, unit-bearing `L_total`, role transfer, EOM, Hamiltonian, or ToE closure; P2849's damping/compression obstruction remains unchanged.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2850/S1800 EML syntax Ltotal guard", "## P2850/S1800 EML syntax `L_total` guard\n\n`P2850/S1800` adds no action term.  EML single-operator representability is a syntax/circuit-basis observation for elementary expressions, not a unit-bearing source density, coupling coefficient, variational chain rule, nonproxy `L_total` insertion, EOM, or Hamiltonian source.\n")
    append_once(AGENTS, "Current EML single-operator external impact guardrail (P2850/S1800, 2026-06-18)", "## Current EML single-operator external impact guardrail (P2850/S1800, 2026-06-18)\n\n- P2850 audits arXiv `2603.21852v2`: EML single-operator universality is relevant as an elementary-expression syntax/circuit compression lens.\n- It does not export typed source laws for strict `beta/eta`, amplitude passage, phase/topological selector data, completion-map semantics, unit-bearing `L_total`, role transfer, EOM, Hamiltonian, or ToE closure.\n- Do not promote EML representability of elementary formulas to kernel bridge closure or role transfer.\n- The next admissible proof-grade move remains one typed bridge-source atom: a strict source law for `eta` and target-independent `beta`, or an amplitude-normalization passage theorem; otherwise preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
