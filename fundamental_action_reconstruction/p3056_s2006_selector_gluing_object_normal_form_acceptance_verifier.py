#!/usr/bin/env python3
"""P3056/S2006: selector gluing-object normal-form acceptance verifier.

P3055 says an unfamiliar selector mechanism must still arrive as one gluing
object with six typed rows.  P3056 constructs that missing object as a precise
normal form and runs a finite acceptance verifier against current clue-carriers.
The point is proof-oriented: do not ask whether a clue has a familiar name; ask
whether one object carries all six morphism slots and the compatibility squares
that would make selector/readout installation non-premise.
"""
from __future__ import annotations

import hashlib, itertools, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3055_s2005_mechanism_agnostic_selector_clue_sheaf import OUT as P3055

OUT = GEN / "p3056_s2006_selector_gluing_object_normal_form_acceptance_verifier.json"
MD = GEN / "p3056_s2006_selector_gluing_object_normal_form_acceptance_verifier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "six_row_gluing_object": r"strict source law.*nonzero signed value|nonpremise origin.*unique polarity|localization/pullback.*variational chain rule|candidate gluing object",
    "single_object_source_to_readout": r"source.*localizer.*polarity|selector/readout.*coupling|one.*object.*source.*coupling|gluing theorem",
    "unit_field_variational_installation": r"field support.*unit-bearing|localized.*action density|variational chain rule|Euler-Lagrange.*coupling",
}

ROWS = [
    "strict_source_law",
    "nonzero_signed_value",
    "nonpremise_origin_or_localizer",
    "unique_polarity_coupling",
    "localization_pullback_to_field_support",
    "unit_bearing_variational_chain_rule",
]

# Current positive clue-carriers rewritten as candidate pieces of one normal-form
# object.  Each row records what is actually exported now, not what would be
# desired after adding a new theorem.
CARRIERS = [
    {
        "carrier": "receiver_phase_area_winding_package",
        "rows": {"nonzero_signed_value"},
        "obstruction": "receiver provenance and polarity remain imported/diagnostic, not strict source/localizer/coupling rows",
    },
    {
        "carrier": "chiral_bispectrum_signed_marker",
        "rows": {"nonzero_signed_value"},
        "obstruction": "nonzero chiral sign lacks source-origin localizer and unique torsor-coupling polarity",
    },
    {
        "carrier": "inversion_odd_sign_algebra",
        "rows": set(),
        "obstruction": "representation type exists only as paired polarity laws without a signed source value selected by artifacts",
    },
    {
        "carrier": "graph_motif_digest_witness",
        "rows": set(),
        "obstruction": "finite nonlabel evidence lacks field pullback, source law, and action/EOM chain rule",
    },
    {
        "carrier": "incidence_localization_readiness",
        "rows": {"localization_pullback_to_field_support"},
        "obstruction": "readiness/localizer skeleton is conditional and not tied to the signed selector source or units",
    },
    {
        "carrier": "formal_action_eom_target_slot",
        "rows": {"unit_bearing_variational_chain_rule"},
        "obstruction": "formal target slot is not sourced by the signed object and does not select polarity",
    },
]

COMPATIBILITY_SQUARES = [
    "source_law_to_signed_value",
    "source_law_to_origin_localizer",
    "signed_value_to_unique_polarity",
    "origin_localizer_to_field_pullback",
    "field_pullback_to_unit_variational_chain_rule",
    "unique_polarity_to_selector_readout",
]


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def carrier_rows(carrier: dict[str, Any]) -> dict[str, bool]:
    return {row: row in carrier["rows"] for row in ROWS}


def subset_pushout() -> list[dict[str, Any]]:
    # Exhaustively computes whether current clue-carriers can be glued into a
    # single object by mere bundling.  A subset may cover rows, but it is accepted
    # only if one carrier already has all rows and all compatibility squares are
    # exported; no current subset satisfies that stronger theorem condition.
    out = []
    for r in range(1, len(CARRIERS) + 1):
        for combo in itertools.combinations(CARRIERS, r):
            names = [c["carrier"] for c in combo]
            covered = set().union(*(c["rows"] for c in combo)) if combo else set()
            row_complete = all(row in covered for row in ROWS)
            single_object = any(all(row in c["rows"] for row in ROWS) for c in combo)
            compatibility_exported = False
            out.append({
                "carriers": names,
                "size": r,
                "covered_rows": sorted(covered),
                "missing_rows": [row for row in ROWS if row not in covered],
                "row_complete_by_bundling": row_complete,
                "single_carrier_has_all_rows": single_object,
                "compatibility_squares_exported": compatibility_exported,
                "accepted_gluing_object": row_complete and single_object and compatibility_exported,
            })
    return out


def normal_form() -> dict[str, Any]:
    return {
        "object": "SelectorGluingObjectNormalForm",
        "typed_rows": ROWS,
        "compatibility_squares": COMPATIBILITY_SQUARES,
        "acceptance_rule": "accepted iff one non-imported object supplies all six rows and exports every compatibility square from source law through selector/readout and unit-bearing action/EOM installation",
        "candidate_carriers": [{"carrier": c["carrier"], "row_flags": carrier_rows(c), "obstruction": c["obstruction"]} for c in CARRIERS],
    }


def build_payload() -> dict[str, Any]:
    p3055 = read_json(P3055)
    grep_rows = content_grep()
    subsets = subset_pushout()
    accepted = [row for row in subsets if row["accepted_gluing_object"]]
    row_complete = [row for row in subsets if row["row_complete_by_bundling"]]
    obligations = [
        {"obligation": "content_first_grep_before_verifier", "satisfied": True, "detail": "three research-content grep lanes were executed before constructing P3056"},
        {"obligation": "construct_six_row_normal_form", "satisfied": True, "detail": "strict source, signed value, localizer, polarity, field pullback, and variational chain-rule slots are explicit"},
        {"obligation": "exhaust_current_carrier_pushout", "satisfied": True, "detail": f"all {len(subsets)} nonempty carrier subsets were enumerated"},
        {"obligation": "find_row_complete_bundle", "satisfied": False, "detail": "no subset of current carriers covers strict_source_law or unique_polarity_coupling"},
        {"obligation": "export_single_object_compatibility_squares", "satisfied": False, "detail": "no carrier exports all six rows or any source-to-readout/action compatibility theorem"},
    ]
    return {
        "status": "P3056_SELECTOR_GLUING_OBJECT_NORMAL_FORM_ACCEPTANCE_VERIFIER_NO_EXPORT",
        "input_hashes": {"P3055": hashlib.sha256(P3055.read_bytes()).hexdigest() if P3055.exists() else None},
        "constructed_theoretical_objects": {
            **normal_form(),
            "content_first_repo_grep": grep_rows,
            "carrier_subset_pushout": subsets,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "typed_rows": len(ROWS),
            "compatibility_squares": len(COMPATIBILITY_SQUARES),
            "candidate_carriers": len(CARRIERS),
            "carrier_subsets_enumerated": len(subsets),
            "row_complete_bundles": len(row_complete),
            "accepted_gluing_objects": len(accepted),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "p3055_status_seen": p3055.get("status"),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3056 constructs the exact six-row selector gluing-object normal form requested by P3055 and exhaustively tests whether current clue-carriers can instantiate it.  The finite pushout over six carriers has 63 nonempty bundles; none covers the strict_source_law or unique_polarity_coupling rows, none is a single object with all rows, and no source-to-readout/action compatibility square is exported.  Thus the unknown selector mechanism remains possible only as a genuinely new object, not as a recombination of current clues.",
            "negative_export_flags": {k: False for k in ["selector_gluing_object_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not bundle current clue-carriers as if bundling were gluing.  The next proof-grade move must introduce one genuinely new source-law/polarity carrier that fills at least the two absent rows strict_source_law and unique_polarity_coupling and proves the compatibility squares to signed value, localizer, field pullback, selector/readout, and unit-bearing variational chain rule; otherwise pivot outside the selector clue lane and preserve the P3048-P3056 no-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3056/S2006 selector gluing-object normal-form acceptance verifier", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- typed rows: `{c['typed_rows']}`",
        f"- compatibility squares: `{c['compatibility_squares']}`",
        f"- candidate carriers: `{c['candidate_carriers']}`",
        f"- carrier subsets enumerated: `{c['carrier_subsets_enumerated']}`",
        f"- row-complete bundles: `{c['row_complete_bundles']}`",
        f"- accepted gluing objects: `{c['accepted_gluing_objects']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3056/S2006 selector gluing-object normal-form acceptance verifier", "## P3056/S2006 selector gluing-object normal-form acceptance verifier\n\n`P3056/S2006` constructs the six-row selector gluing-object normal form demanded after `P3055`: strict source law, nonzero signed value, nonpremise origin/localizer, unique polarity coupling, localization/pullback to field support, and unit-bearing variational chain rule, plus six compatibility squares from source through selector/readout and action/EOM installation.  The finite pushout over six current clue-carriers enumerates all `63` nonempty bundles; `0` bundles cover all rows, `0` carriers contain all rows, and `0` compatibility-square exports are present.  No selector gluing object, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3056/S2006 selector gluing-object normal-form `L_total` guard", "## P3056/S2006 selector gluing-object normal-form `L_total` guard\n\n`P3056/S2006` adds no physical `L_total` term.  It is a normal-form acceptance verifier: current carriers do not provide a single signed source with field pullback, unit-bearing coupling, and a nonproxy variational chain rule.\n")
    append_once(AGENTS, "Current selector gluing-object normal-form guardrail (P3056/S2006, 2026-06-23)", "## Current selector gluing-object normal-form guardrail (P3056/S2006, 2026-06-23)\n\n- P3056 constructs the exact six-row normal form for the unknown selector gluing object required after P3055 and tests current clue-carriers by finite pushout.\n- All `63` nonempty carrier bundles are enumerated; none covers the absent `strict_source_law` and `unique_polarity_coupling` rows, none is a single all-row object, and no compatibility squares from source to selector/readout/action installation are exported.\n- Do not promote clue bundling, receiver/chiral signs, graph motifs, incidence readiness, or formal action slots to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must introduce one genuinely new source-law/polarity carrier and prove the compatibility squares, or pivot outside the selector clue lane while preserving the P3048-P3056 bounded no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
