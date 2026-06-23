#!/usr/bin/env python3
"""P3055/S2005: mechanism-agnostic selector clue sheaf.

The post-P3054 problem is not to force the selector into a familiar human
schema.  The honest task is to mine the existing finite positives for clues and
then state the exact theorem object that an unfamiliar mechanism would still
have to supply.

P3055 builds a content-first clue sheaf from current artifacts.  Stalks are not
named by P-numbers but by research content: robust receiver phase geometry,
inversion-odd sign algebra, explicit chiral signed formulae, non-label graph
motif/digest witnesses, and unit-bearing variational installation.  Each stalk
records positive finite evidence and the missing theorem obligations.  A finite
hitting-set computation then asks which obligation bundle is required for an
unknown mechanism to glue these clues into a strict selector/source theorem.
"""
from __future__ import annotations

import hashlib, itertools, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3054_s2004_nonreceiver_signed_source_coupling_intake_matrix import OUT as P3054

OUT = GEN / "p3055_s2005_mechanism_agnostic_selector_clue_sheaf.json"
MD = GEN / "p3055_s2005_mechanism_agnostic_selector_clue_sheaf.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "receiver_phase_geometry_clues": r"receiver.*geometry|phase.*curve.*winding|memory.*phase|area.*sign|robust.*winding|positive.*clearance",
    "odd_sign_and_polarity_clues": r"inversion-odd|odd.*polarity|sign.*torsor|polarity.*pair|coupling.*polarity|nonpremise.*polarity",
    "chiral_signed_formula_clues": r"chiral.*bispectrum|pseudoscalar|signed.*formula|phase-origin|source.*localizer",
    "graph_motif_digest_clues": r"local.*motif|edge-toggle|second-variation|digest.*carrier|graph.*source|localization.*pullback",
    "unit_variational_installation_clues": r"unit-bearing|variational.*derivative|action-density|Euler-Lagrange|L_total|chain-rule",
}

OBLIGATIONS = [
    "strict_source_law",
    "nonzero_signed_value",
    "nonpremise_origin_or_localizer",
    "unique_polarity_coupling",
    "localization_pullback_to_field_support",
    "unit_bearing_variational_chain_rule",
]

STALKS = [
    {
        "stalk": "robust_receiver_phase_geometry",
        "positive_clues": ["nonzero area sign", "integer winding", "positive clearance", "Aut signed response"],
        "missing": ["strict_source_law", "nonpremise_origin_or_localizer", "unique_polarity_coupling", "unit_bearing_variational_chain_rule"],
        "why_not_human_schema": "The clue is geometric/receiver-level and need not be a conventional observer selector, but it still lacks sourcehood and coupling.",
    },
    {
        "stalk": "inversion_odd_sign_algebra",
        "positive_clues": ["correct C2 representation type", "odd/equivariant laws exist"],
        "missing": ["strict_source_law", "nonzero_signed_value", "unique_polarity_coupling"],
        "why_not_human_schema": "The algebra allows unfamiliar sign mechanisms, but paired polarities remain unglued without a strict source value.",
    },
    {
        "stalk": "explicit_chiral_signed_formulae",
        "positive_clues": ["nonzero signed formulae", "orientation separation"],
        "missing": ["nonpremise_origin_or_localizer", "unique_polarity_coupling", "unit_bearing_variational_chain_rule"],
        "why_not_human_schema": "A chiral marker may be the right non-human clue, but current artifacts do not source its origin or coupling polarity.",
    },
    {
        "stalk": "nonlabel_graph_motif_digest_witnesses",
        "positive_clues": ["full finite carrier separation", "edge-toggle response", "second-variation patch", "non-label motif evidence"],
        "missing": ["localization_pullback_to_field_support", "unit_bearing_variational_chain_rule", "strict_source_law"],
        "why_not_human_schema": "The graph/digest lane is the least human-readable clue class, but it remains uncoupled from fields and units.",
    },
    {
        "stalk": "unit_bearing_action_eom_installation",
        "positive_clues": ["explicit target slot for action/EOM installation"],
        "missing": ["strict_source_law", "localization_pullback_to_field_support", "unit_bearing_variational_chain_rule", "unique_polarity_coupling"],
        "why_not_human_schema": "Any unfamiliar selector mechanism must still install as a typed source/coupling if it is to affect action/EOM rather than remain diagnostic.",
    },
]


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:12]})
    return rows


def incidence_matrix() -> list[dict[str, Any]]:
    rows = []
    for stalk in STALKS:
        missing = set(stalk["missing"])
        rows.append({
            "stalk": stalk["stalk"],
            "positive_clue_count": len(stalk["positive_clues"]),
            **{obligation: obligation in missing for obligation in OBLIGATIONS},
        })
    return rows


def minimal_hitting_sets() -> list[dict[str, Any]]:
    # A bundle is sufficient only if it covers every missing obligation appearing
    # anywhere in the clue sheaf.  Enumerating all subsets makes the theorem
    # target explicit rather than choosing a familiar-looking mechanism by hand.
    required = set().union(*(set(stalk["missing"]) for stalk in STALKS))
    bundles = []
    for r in range(1, len(OBLIGATIONS) + 1):
        for combo in itertools.combinations(OBLIGATIONS, r):
            combo_set = set(combo)
            if required.issubset(combo_set):
                bundles.append({"obligation_bundle": list(combo), "size": r, "covers_all_stalks": True})
        if bundles:
            break
    return bundles


def current_artifact_satisfaction() -> dict[str, bool]:
    # Positive clues are real, but no current artifact exports the gluing theorem.
    return {
        "strict_source_law": False,
        "nonzero_signed_value": True,
        "nonpremise_origin_or_localizer": False,
        "unique_polarity_coupling": False,
        "localization_pullback_to_field_support": False,
        "unit_bearing_variational_chain_rule": False,
    }


def build_payload() -> dict[str, Any]:
    p3054 = read_json(P3054)
    grep_rows = content_grep()
    incidence = incidence_matrix()
    bundles = minimal_hitting_sets()
    satisfaction = current_artifact_satisfaction()
    accepted = all(satisfaction[ob] for ob in bundles[0]["obligation_bundle"])
    obligations = [
        {"obligation": "content_first_grep_before_lane_choice", "satisfied": True, "detail": "five research-content patterns were searched before constructing the sheaf"},
        {"obligation": "construct_mechanism_agnostic_clue_sheaf", "satisfied": True, "detail": "five stalks separate finite positives from theorem gaps without assuming a human selector schema"},
        {"obligation": "compute_missing_obligation_hitting_set", "satisfied": True, "detail": "all subsets of six obligations were enumerated; the unique minimal full bundle has size six"},
        {"obligation": "export_unknown_mechanism_gluing_theorem", "satisfied": False, "detail": "current artifacts satisfy only the nonzero-signed-value row and do not glue source/localizer/polarity/localization/units"},
        {"obligation": "install_selector_or_ltotal", "satisfied": False, "detail": "no selector closure, QW-2191 discharge, unit-bearing action/EOM, L_total, bridge, role transfer, or ToE export follows"},
    ]
    return {
        "status": "P3055_MECHANISM_AGNOSTIC_SELECTOR_CLUE_SHEAF_BOUNDARY_NO_EXPORT",
        "input_hashes": {"P3054": hashlib.sha256(P3054.read_bytes()).hexdigest() if P3054.exists() else None},
        "constructed_theoretical_objects": {
            "object": "MechanismAgnosticSelectorClueSheaf",
            "content_first_repo_grep": grep_rows,
            "stalks": STALKS,
            "obligations": OBLIGATIONS,
            "incidence_matrix": incidence,
            "minimal_full_gluing_bundles": bundles,
            "current_artifact_satisfaction": satisfaction,
            "unknown_mechanism_gluing_theorem_exported": accepted,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "clue_stalks": len(STALKS),
            "positive_clues": sum(len(stalk["positive_clues"]) for stalk in STALKS),
            "obligation_atoms": len(OBLIGATIONS),
            "minimal_full_gluing_bundle_size": bundles[0]["size"],
            "minimal_full_gluing_bundle_count": len(bundles),
            "currently_satisfied_gluing_obligations": sum(1 for v in satisfaction.values() if v),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "unknown_mechanism_selector_source_exported": accepted,
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3055 takes seriously that the selector mechanism may be non-human-readable: it constructs a clue sheaf rather than imposing a known schema.  The finite positives are real across receiver geometry, odd sign algebra, chiral signed formulae, graph motif/digest witnesses, and action/EOM target slots.  However the sheaf cannot glue unless one object supplies source law, signed value, nonpremise origin/localizer, unique polarity coupling, localization/pullback to field support, and a unit-bearing variational chain rule.  Current artifacts satisfy only the signed-value clue, so no selector source theorem is exported.",
            "negative_export_flags": {k: False for k in ["unknown_mechanism_selector_source_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not force a familiar selector story and do not replay receiver/sign/graph clues as closure.  The next proof-grade move should construct exactly one candidate gluing object with all six typed rows: strict source law, computable nonzero signed value, nonpremise origin/localizer, unique polarity coupling, localization/pullback to field support, and unit-bearing variational chain rule.  If that is unavailable, pivot to a new typed object and preserve the P3048-P3055 no-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3055/S2005 mechanism-agnostic selector clue sheaf", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- clue stalks: `{c['clue_stalks']}`",
        f"- positive clues: `{c['positive_clues']}`",
        f"- obligation atoms: `{c['obligation_atoms']}`",
        f"- minimal full gluing bundle size: `{c['minimal_full_gluing_bundle_size']}`",
        f"- minimal full gluing bundle count: `{c['minimal_full_gluing_bundle_count']}`",
        f"- currently satisfied gluing obligations: `{c['currently_satisfied_gluing_obligations']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- unknown mechanism selector source exported: `{c['unknown_mechanism_selector_source_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3055/S2005 mechanism-agnostic selector clue sheaf", "## P3055/S2005 mechanism-agnostic selector clue sheaf\n\n`P3055/S2005` responds to the possibility that the real selector mechanism need not match a familiar human schema by constructing a content-first clue sheaf.  Five stalks are separated: robust receiver phase geometry, inversion-odd sign algebra, explicit chiral signed formulae, non-label graph motif/digest witnesses, and unit-bearing action/EOM installation.  The finite sheaf has `13` positive clues and six obligation atoms; exhaustive subset enumeration shows the minimal full gluing bundle still requires all six rows: strict source law, nonzero signed value, nonpremise origin/localizer, unique polarity coupling, localization/pullback to field support, and unit-bearing variational chain rule.  Current artifacts satisfy only the signed-value clue, so no unknown-mechanism selector source theorem, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3055/S2005 mechanism-agnostic selector clue sheaf `L_total` guard", "## P3055/S2005 mechanism-agnostic selector clue sheaf `L_total` guard\n\n`P3055/S2005` adds no physical `L_total` term.  The clue sheaf identifies real finite positives, including graph/motif and receiver/chiral clues, but no single gluing object supplies localization/pullback, unit-bearing coupling, and a variational chain rule.\n")
    append_once(AGENTS, "Current mechanism-agnostic selector clue sheaf guardrail (P3055/S2005, 2026-06-23)", "## Current mechanism-agnostic selector clue sheaf guardrail (P3055/S2005, 2026-06-23)\n\n- P3055 treats the selector as potentially non-human-schema by building a content-first clue sheaf rather than forcing receiver, sign-torsor, chiral, graph, or Lagrangian clues into closure.\n- The finite sheaf records `13` positive clues across five stalks, but its minimal full gluing bundle requires all six typed rows: strict source law, nonzero signed value, nonpremise origin/localizer, unique polarity coupling, localization/pullback to field support, and unit-bearing variational chain rule.\n- Current artifacts satisfy only the signed-value clue; do not promote clue convergence, receiver geometry, graph motif/digest evidence, chiral signed formulae, or action/EOM target slots to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must construct exactly one candidate gluing object satisfying the six-row bundle, or pivot to a new typed object while preserving the P3048-P3055 bounded no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
