#!/usr/bin/env python3
"""P3062/S2012: sigma-selector source-law candidate audit.

P3061 isolated the next missing object: a concrete strict source law computing
a nonzero non-premise sigma_selector value coupled to G_selector.  P3062 makes
that obligation proof-grade by auditing named current candidate source laws
against the full acceptance boundary instead of adding another formal sigma row.
"""
from __future__ import annotations

import hashlib, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3061_s2011_sign_odd_source_value_normalizer_acceptance_matrix import OUT as P3061

OUT = GEN / "p3062_s2012_sigma_selector_source_law_candidate_audit.json"
MD = GEN / "p3062_s2012_sigma_selector_source_law_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "concrete_sigma_selector_source_law": r"sigma_selector|sigma selector|signed source value|source law.*signed|strict source law.*sigma",
    "nonpremise_signed_source_coupling": r"non-premise.*signed source|nonpremise.*signed source|coupled to G_selector|coupling.*G_selector",
    "candidate_sign_carriers": r"boundary.*cocycle|chiral.*bispectrum|receiver.*winding|Pontryagin|eta.*asymmetry|Levi-Civita",
}

ACCEPTANCE_CRITERIA = (
    "strict_source_law_exported",
    "nonzero_signed_value_computed",
    "nonpremise_origin_or_localizer",
    "coupled_to_G_selector",
    "not_only_convention_or_readout",
)

CANDIDATES = [
    {
        "candidate": "boundary_cocycle_orientation_sign",
        "evidence_kind": "topological orientation clue",
        "criteria": {
            "strict_source_law_exported": False,
            "nonzero_signed_value_computed": True,
            "nonpremise_origin_or_localizer": False,
            "coupled_to_G_selector": False,
            "not_only_convention_or_readout": False,
        },
    },
    {
        "candidate": "chiral_bispectrum_Im_B_1_5_sign",
        "evidence_kind": "nonzero chiral readout marker",
        "criteria": {
            "strict_source_law_exported": False,
            "nonzero_signed_value_computed": True,
            "nonpremise_origin_or_localizer": False,
            "coupled_to_G_selector": False,
            "not_only_convention_or_readout": True,
        },
    },
    {
        "candidate": "receiver_winding_sign",
        "evidence_kind": "receiver/readout signed diagnostic",
        "criteria": {
            "strict_source_law_exported": False,
            "nonzero_signed_value_computed": True,
            "nonpremise_origin_or_localizer": False,
            "coupled_to_G_selector": False,
            "not_only_convention_or_readout": False,
        },
    },
    {
        "candidate": "levi_civita_orientation_density",
        "evidence_kind": "orientation-density template",
        "criteria": {
            "strict_source_law_exported": False,
            "nonzero_signed_value_computed": False,
            "nonpremise_origin_or_localizer": False,
            "coupled_to_G_selector": False,
            "not_only_convention_or_readout": False,
        },
    },
    {
        "candidate": "pontryagin_or_anomaly_density",
        "evidence_kind": "pseudoscalar density template",
        "criteria": {
            "strict_source_law_exported": False,
            "nonzero_signed_value_computed": False,
            "nonpremise_origin_or_localizer": False,
            "coupled_to_G_selector": False,
            "not_only_convention_or_readout": True,
        },
    },
    {
        "candidate": "eta_spectral_asymmetry",
        "evidence_kind": "spectral asymmetry template",
        "criteria": {
            "strict_source_law_exported": False,
            "nonzero_signed_value_computed": False,
            "nonpremise_origin_or_localizer": False,
            "coupled_to_G_selector": False,
            "not_only_convention_or_readout": True,
        },
    },
    {
        "candidate": "coefficient_sign_convention",
        "evidence_kind": "normalization convention",
        "criteria": {
            "strict_source_law_exported": False,
            "nonzero_signed_value_computed": True,
            "nonpremise_origin_or_localizer": False,
            "coupled_to_G_selector": False,
            "not_only_convention_or_readout": False,
        },
    },
]


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def evaluate_candidate(row: dict[str, Any]) -> dict[str, Any]:
    criteria = row["criteria"]
    missing = [name for name in ACCEPTANCE_CRITERIA if not criteria[name]]
    return {
        **row,
        "satisfied_criteria": sum(1 for name in ACCEPTANCE_CRITERIA if criteria[name]),
        "missing_criteria": missing,
        "accepted_as_sigma_selector_source_law": not missing,
        "blocker": "accepted" if not missing else "missing " + ", ".join(missing),
    }


def build_payload() -> dict[str, Any]:
    p3061 = read_json(P3061)
    grep_rows = content_grep()
    candidate_rows = [evaluate_candidate(row) for row in CANDIDATES]
    obligations = [
        {"obligation": "content_first_grep_before_sigma_source_audit", "satisfied": True, "detail": "searched by signed-source, G_selector coupling, and candidate carrier content"},
        {"obligation": "construct_sigma_source_law_acceptance_boundary", "satisfied": True, "detail": "five criteria separate source law, value, origin, coupling, and convention/readout status"},
        {"obligation": "audit_named_current_candidate_classes", "satisfied": True, "detail": "seven candidate source-law/carrier classes are checked"},
        {"obligation": "export_concrete_sigma_selector_source_law", "satisfied": False, "detail": "zero candidates satisfy all five acceptance criteria"},
        {"obligation": "export_selector_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE closure follows"},
    ]
    accepted = [row for row in candidate_rows if row["accepted_as_sigma_selector_source_law"]]
    return {
        "status": "P3062_SIGMA_SELECTOR_SOURCE_LAW_CANDIDATE_AUDIT_NO_EXPORT",
        "input_hashes": {"P3061": hashlib.sha256(P3061.read_bytes()).hexdigest() if P3061.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "sigma_selector_source_law_acceptance_boundary": {
                "object": "ConcreteSigmaSelectorSourceLawAcceptanceBoundary",
                "target_missing_object": "strict_source_law_computing_nonzero_nonpremise_sigma_selector_coupled_to_G_selector",
                "criteria": list(ACCEPTANCE_CRITERIA),
                "acceptance_rule": "accept exactly a candidate satisfying all five criteria; signed readouts or conventions alone are insufficient",
            },
            "candidate_rows": candidate_rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "acceptance_criteria": len(ACCEPTANCE_CRITERIA),
            "candidate_classes_audited": len(candidate_rows),
            "candidates_with_nonzero_signed_value": sum(1 for row in candidate_rows if row["criteria"]["nonzero_signed_value_computed"]),
            "candidates_coupled_to_G_selector": sum(1 for row in candidate_rows if row["criteria"]["coupled_to_G_selector"]),
            "accepted_sigma_selector_source_laws": len(accepted),
            "p3061_status_seen": p3061.get("status"),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3062 audits seven named current candidate classes against the concrete sigma_selector source-law boundary.  Four candidates carry a nonzero signed value, but zero are coupled to G_selector and zero satisfy all five criteria.  The result preserves the P3061 distinction: formal sign-odd normalization would work only after an exported strict non-premise sigma source law/value is supplied.  No G_selector, QW-2191 discharge, selector closure, L_total, bridge/role transfer, or ToE closure is exported.",
            "negative_export_flags": {k: False for k in ["concrete_sigma_selector_source_law_exported", "strict_nonpremise_sigma_selector_exported", "global_coefficient_sign_normalization_exported", "actual_g_selector_formula_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not recycle the seven audited carrier names as sigma_selector closure evidence.  The next proof-grade move must either derive an explicit coupling theorem from one concrete signed source law/value into G_selector, or pivot to another P3057 atom with content-first grep and bounded finite acceptance criteria.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3062/S2012 sigma-selector source-law candidate audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- acceptance criteria: `{c['acceptance_criteria']}`", f"- candidate classes audited: `{c['candidate_classes_audited']}`", f"- candidates with nonzero signed value: `{c['candidates_with_nonzero_signed_value']}`", f"- candidates coupled to G_selector: `{c['candidates_coupled_to_G_selector']}`", f"- accepted sigma-selector source laws: `{c['accepted_sigma_selector_source_laws']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3062/S2012 sigma-selector source-law candidate audit", "## P3062/S2012 sigma-selector source-law candidate audit\n\n`P3062/S2012` audits the exact missing object named by `P3061`: a strict source law computing a nonzero non-premise `sigma_selector` value coupled to `G_selector`.  The acceptance boundary has `5` criteria and checks `7` named current candidate classes: boundary-cocycle orientation sign, chiral-bispectrum sign, receiver-winding sign, Levi-Civita orientation density, Pontryagin/anomaly density, eta/spectral asymmetry, and coefficient-sign convention.  `4` candidates carry a nonzero signed value, but `0` are coupled to `G_selector` and `0` satisfy all criteria.  No concrete `sigma_selector` source law, actual `G_selector`, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3062/S2012 sigma-selector source-law audit `L_total` guard", "## P3062/S2012 sigma-selector source-law audit `L_total` guard\n\n`P3062/S2012` adds no physical `L_total` term.  It is a source-law acceptance audit; no audited candidate supplies a unit-bearing signed source/action/EOM carrier or nonproxy variational chain rule coupled to `G_selector`.\n")
    append_once(AGENTS, "Current sigma-selector source-law candidate audit guardrail (P3062/S2012, 2026-06-23)", "## Current sigma-selector source-law candidate audit guardrail (P3062/S2012, 2026-06-23)\n\n- P3062 audits the concrete missing object after P3061: a strict source law computing a nonzero non-premise `sigma_selector` value coupled to `G_selector`.\n- The audit checks `7` named candidate classes against `5` criteria; `4` candidates carry nonzero signed values, but `0` are coupled to `G_selector` and `0` are accepted as concrete sigma-selector source laws.\n- Do not recycle boundary-cocycle, chiral-bispectrum, receiver-winding, Levi-Civita, Pontryagin/anomaly, eta-asymmetry, or coefficient-convention rows as `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must derive an explicit coupling theorem from one concrete signed source law/value into `G_selector`, or pivot to another P3057 atom while preserving the P3048-P3062 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
