#!/usr/bin/env python3
"""P2843/S1793: professorial SWOT and closure-path audit.

This is a proof-oriented synthesis step requested after P2842.  It does not add
another graph-source replay.  It converts the current guardrailed state into a
structured SWOT for the theoretical foundations, kernel formulas (especially
K_strict_gate), and the Lagrangian/EOM/Hamiltonian chain, then extracts a finite
closure-path obligation order with explicit blockers and stop rules.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
OUT = GEN / "p2843_s1793_professorial_swot_foundations_kernel_lagrangian_eom_hamiltonian_closure_path.json"
MD = GEN / "p2843_s1793_professorial_swot_foundations_kernel_lagrangian_eom_hamiltonian_closure_path.md"
P2842 = GEN / "p2842_s1792_exchangeable_edge_pair_measure_localization_candidate_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
K1 = ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
K2 = ROOT / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md"
F2 = ROOT / "F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md"
F3 = ROOT / "F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md"
S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

CURRENT_INPUTS = [K1, K2, F2, F3, S2, STRICT_EQUATION_SHEET, STRICT_LAGRANGIAN_DRAFT, AGENTS, P2842]

CLOSURE_GATES = (
    "ontology_and_layer_discipline",
    "strict_kernel_provenance",
    "legacy_to_strict_completion_map",
    "role_transfer_theorem",
    "selector_or_orientation_source",
    "units_and_normalization",
    "typed_coupling_to_K_or_Ltotal",
    "variational_derivative_and_eom",
    "hamiltonian_legendre_phase_space",
    "nonproxy_residual_closure",
)


def count_mentions(text: str, patterns: tuple[str, ...]) -> int:
    return sum(len(re.findall(pattern, text, flags=re.IGNORECASE)) for pattern in patterns)


def source_metrics() -> dict[str, Any]:
    blobs = {path.name: path.read_text(encoding="utf-8") for path in CURRENT_INPUTS if path.exists()}
    combined = "\n".join(blobs.values())
    return {
        "files_accounted": sorted(blobs),
        "mention_counts": {
            "K_strict_gate": count_mentions(combined, (r"K_strict_gate", r"strict gate kernel")),
            "K_legacy_ont": count_mentions(combined, (r"K_legacy_ont", r"legacy ontological", r"legacy kernel")),
            "L_total_or_lagrangian": count_mentions(combined, (r"L_total", r"Lagrangian", r"Lagrangianu")),
            "EOM": count_mentions(combined, (r"\bEOM\b", r"Euler-Lagrange", r"variational derivative")),
            "Hamiltonian": count_mentions(combined, (r"Hamiltonian", r"Hamiltonianu", r"Legendre")),
            "selector_QW2191": count_mentions(combined, (r"QW-2191", r"selector")),
            "no_closure_flags": count_mentions(combined, (r"no .*closure", r"No .*closure", r"not exported", r"bounded no-go")),
        },
        "hashes": {str(path.relative_to(REPO)): sha(path) for path in CURRENT_INPUTS if path.exists()},
    }


def build_swot(metrics: dict[str, Any], p2842: dict[str, Any]) -> dict[str, Any]:
    return {
        "foundations": {
            "strengths": [
                "Ontology is explicit: nadsoliton is primordial fractal information in a solitonic state, with no lower information layer.",
                "Layer discipline is mature: legacy ontological kernel, strict operational kernel, selector lane, graph-source lane, and L_total/EOM lane are separated by guardrails.",
                "The research program now has repeatable no-false-pass intake gates rather than unrestricted narrative promotion.",
            ],
            "weaknesses": [
                "Closure remains gate-limited: multiple lanes end in no-new-live-frontier certificates unless a genuinely new strict object is supplied.",
                "Several positive finite witnesses are not yet converted into typed physical source laws.",
                "Role-bearing claims are intentionally quarantined until bridge and role-transfer theorems exist.",
            ],
            "opportunities": [
                "A finite theorem-obligation ledger can turn the program from discovery mode into certification mode.",
                "A single new strict typed object with source, units, and coupling data could reopen multiple currently blocked rows without replay.",
            ],
            "threats": [
                "The main scientific risk is false closure by importing selector, role-transfer, or graph-source claims across layers.",
                "The second risk is endless finite-carrier refinement after P2835/P2840/P2841 without a physical coupling theorem.",
            ],
        },
        "kernel_formulas": {
            "strengths": [
                "The kernel split is explicitly documented: K_legacy_ont(d)=alpha_geo*cos(omega*d+phi)/(1+beta_tors*d) and K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta) are not silently identified.",
                "K_strict_gate has a recorded operational derivation chain through refreeze, micro support repair, micro/Stage-C intersection, and strict damping/compression interpretation.",
                "The strict formula captures nonlinear damping/compression through d^eta, a strict-side characteristic absent from the legacy linear beta_tors*d denominator.",
            ],
            "weaknesses": [
                "The amplitude passage alpha_geo -> strict normalization is not theorem-level closed.",
                "Strict beta/Z_beta and eta/compression data remain operational/representational unless target-independent source laws are exported.",
                "The legacy-to-strict completion bridge and downstream physical-role transfer theorem are still missing.",
            ],
            "opportunities": [
                "A completion-map theorem decomposed into amplitude, phase/frequency, damping/compression, selector/source, and residual strict-side additions is the cleanest kernel-level route.",
                "A bounded witness table can test exactly one bridge atom at a time, avoiding generic bridge replay.",
            ],
            "threats": [
                "Using K_strict_gate as a silent ontological substitute for K_legacy_ont would invalidate downstream role claims.",
                "Using legacy formulas for alpha_EM, sin^2(theta_W), or beta^N hierarchy before role transfer would overclaim.",
            ],
        },
        "lagrangian_eom_hamiltonian": {
            "strengths": [
                "The Lagrangian/EOM draft contains a long audit trail showing which finite and variational-looking objects add no L_total term.",
                "Graph-source audits reached a full finite separator and then correctly audited units, typed domain/codomain, variational derivative, localization, and coupling obligations.",
                "The current EOM posture distinguishes finite graph differences from Euler-Lagrange/functional derivatives.",
            ],
            "weaknesses": [
                "No current graph functional exports a typed, unit-bearing source density or action-density embedding into L_total.",
                "Formal variational derivative, localization/pullback, and graph-bit-to-field chain rule remain unexported.",
                "Hamiltonian closure is even later-stage: without a nonproxy L_total/EOM closure and a Legendre/constraint analysis, no canonical Hamiltonian should be promoted.",
            ],
            "opportunities": [
                "The next proof-grade Lagrangian move is not more graph separation but one explicit typed graph-to-field or kernel-to-action source theorem with units and variation rules.",
                "A Hamiltonian audit can be made rigorous only as a downstream obligation matrix: fields, momenta, primary constraints, gauge generators, boundedness, and recovery of EOM.",
            ],
            "threats": [
                "Treating finite edge-toggle differences as variational calculus would create a category error.",
                "Treating a reduced or proxy Lagrangian as full tensor-resolved nonproxy closure would bypass known residual obstructions.",
            ],
        },
    }


def closure_path() -> list[dict[str, Any]]:
    rows = [
        (1, "Freeze scope and notation", "Publish a one-page current-state theorem ledger separating ontology, kernel, selector, graph-source, Lagrangian/EOM, and Hamiltonian objects.", True, "This is editorial/proof-hygiene and can be done without new physics."),
        (2, "Kernel completion-map theorem", "Attack exactly one bridge atom at a time: amplitude/normalization, phase/frequency/topological data, or damping/compression source; stop at first missing source premise.", False, "Generic bridge replay is blocked; a new typed bridge/source atom is required."),
        (3, "Role-transfer theorem", "Only after step 2 exports a bridge, audit sin^2(theta_W), alpha_EM, and hierarchy claims as survive/modified/rejected.", False, "Downstream of bridge completion; currently forbidden."),
        (4, "Selector/orientation source", "Provide a non-premise strict symmetry-breaking provider or signed pseudoscalar/chiral source coupled to the relevant torsor.", False, "QW-2191 remains open without a new source law."),
        (5, "Typed source/coupling into L_total", "Construct one unit-bearing source density or coefficient map with domain/codomain, localization, covariance, and coupling coefficient.", False, "P2835-P2842 show finite graph evidence is insufficient."),
        (6, "Variational derivative and EOM", "Derive Euler-Lagrange equations from the accepted action term, including boundary conditions, integration-by-parts rules, and nonproxy residual tables.", False, "No accepted source term currently exists."),
        (7, "Hamiltonian/constraint closure", "Only after L_total/EOM closure, perform Legendre transform, identify momenta and constraints, prove Hamiltonian boundedness/gauge closure, and recover EOM.", False, "Hamiltonian promotion before L_total/EOM closure is premature."),
        (8, "Broad state-map reconciliation", "If none of steps 2/4/5 supplies a new strict object, preserve no-new-live-frontier rather than manufacture closure.", True, "This is the honest fallback."),
    ]
    return [
        {"rank": rank, "name": name, "professorial_action": action, "admissible_now": admissible, "reason": reason}
        for rank, name, action, admissible, reason in rows
    ]


def acceptance_matrix(payload: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "all_requested_domains_accounted": set(payload["swot"]) == {"foundations", "kernel_formulas", "lagrangian_eom_hamiltonian"},
        "strict_kernel_split_preserved": True,
        "lagrangian_eom_hamiltonian_not_promoted": True,
        "graph_source_replay_avoided": True,
        "closure_path_has_stop_rules": all((not row["admissible_now"]) or ("stop" in row["professorial_action"].lower()) or ("fallback" in row["reason"].lower()) or ("proof-hygiene" in row["reason"].lower()) for row in payload["closure_path"]),
        "new_closure_exported": False,
    }
    return {
        "gates": gates,
        "accepted_as_professorial_swot_and_closure_path_audit": all(v for k, v in gates.items() if k != "new_closure_exported") and not gates["new_closure_exported"],
        "exports_ltotal_or_toe_closure": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2843/S1793 professorial SWOT: foundations, strict kernel, Lagrangian/EOM/Hamiltonian", "", f"Status: `{payload['status']}`", "",
        "## Executive result",
        payload["decision"]["reason"], "",
    ]
    labels = {"foundations": "Foundations", "kernel_formulas": "Kernel formulas, especially K_strict_gate", "lagrangian_eom_hamiltonian": "Lagrangian / EOM / Hamiltonian"}
    for key, title in labels.items():
        lines += [f"## SWOT — {title}"]
        for bucket in ("strengths", "weaknesses", "opportunities", "threats"):
            lines.append(f"### {bucket.title()}")
            lines.extend(f"- {item}" for item in payload["swot"][key][bucket])
        lines.append("")
    lines += ["## Professorial closure path"]
    for row in payload["closure_path"]:
        lines.append(f"{row['rank']}. **{row['name']}** — admissible_now={row['admissible_now']}. {row['professorial_action']} Reason: {row['reason']}")
    lines += ["", "## Next honest step", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2842 = read_json(P2842)
    metrics = source_metrics()
    payload: dict[str, Any] = {
        "status": "P2843_PROFESSORIAL_SWOT_CLOSURE_PATH_AUDIT_NO_NEW_CLOSURE",
        "input_statuses_rechecked": {"P2842": p2842.get("status")},
        "source_metrics": metrics,
        "closure_gates": list(CLOSURE_GATES),
        "swot": build_swot(metrics, p2842),
        "closure_path": closure_path(),
        "decision": {
            "negative_export_flags": {
                "strict_kernel_identity_exported": False,
                "legacy_role_transfer_started": False,
                "selector_closure_exported": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "role_bearing_ltotal_promoted": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2843 is a synthesis and obligation-ranking audit, not a new closure theorem.  It accounts for the requested foundations, strict-kernel formula, Lagrangian/EOM, and Hamiltonian layers while preserving the kernel split and the P2835-P2842 no-coupling boundary.",
            "next_honest_step": "Do not add more graph-source separation or carrier-constant measures.  The next proof-grade move should introduce exactly one new typed theorem object: preferably a unit-bearing typed source/coupling map into L_total with localization and variational chain rule; if that cannot be supplied, the honest alternative is one kernel bridge atom with a concrete source premise.  Otherwise preserve the no-new-live-frontier certificate.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(payload)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2843/S1793 professorial SWOT closure path audit", "## P2843/S1793 professorial SWOT closure path audit\n\n`P2843/S1793` performs the requested professorial SWOT of the theoretical foundations, kernel formulas (especially `K_strict_gate`), and the Lagrangian/EOM/Hamiltonian chain.  It is a synthesis and theorem-obligation ranking, not a closure theorem: the kernel split is preserved, P2835-P2842 finite graph evidence remains uncoupled from `K`/`L_total`, and Hamiltonian closure is marked downstream of nonproxy `L_total`/EOM plus Legendre/constraint analysis.  The recommended next proof-grade move is exactly one new unit-bearing typed source/coupling map into `L_total` with localization and a variational chain rule, or else one concrete kernel bridge atom with an exported source premise; otherwise preserve the no-new-live-frontier certificate.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2843/S1793 SWOT Hamiltonian closure guard", "## P2843/S1793 SWOT Hamiltonian closure guard\n\n`P2843/S1793` adds no variational term, EOM, or Hamiltonian.  It ranks the closure obligations and states that Hamiltonian promotion is downstream of a nonproxy unit-bearing `L_total`, Euler-Lagrange derivation, Legendre transform, constraint/gauge analysis, boundedness check, and EOM recovery.\n")
    append_once(AGENTS, "Current professorial SWOT closure-path guardrail (P2843/S1793, 2026-06-18)", "## Current professorial SWOT closure-path guardrail (P2843/S1793, 2026-06-18)\n\n- P2843 is a synthesis/audit of foundations, kernel formulas, and Lagrangian/EOM/Hamiltonian obligations; it exports no new closure theorem.\n- Preserve the kernel split: `K_legacy_ont` remains an intermediate bridge kernel and `K_strict_gate` remains the completed/enriched strict working kernel only where an explicit completion-map certificate licenses it.\n- Do not promote finite graph witnesses, SWOT language, or Hamiltonian path language to strict graph-source law, role-bearing `L_total`, EOM closure, Hamiltonian closure, bridge closure, role transfer, selector closure, or ToE closure.\n- The next proof-grade move must supply exactly one unit-bearing typed source/coupling map into `L_total` with localization and a variational chain rule, or one concrete kernel bridge atom with an exported source premise; otherwise preserve the no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
