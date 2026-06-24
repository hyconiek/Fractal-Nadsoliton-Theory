#!/usr/bin/env python3
"""P3063/S2013: signed-source to G_selector coupling theorem matrix.

P3062 showed that several carriers have signed values but no exported coupling
into G_selector.  P3063 attacks exactly that gap by constructing the finite
normal form for an explicit coupling theorem from a concrete signed source
law/value into G_selector.  It also records a lay analogy boundary: the strict
kernel may be compared to a law/rulebook only operationally, while the selector
is closer to an initial branch/choice rule than to a proven "simulation start".
"""
from __future__ import annotations

import hashlib, json, subprocess
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3062_s2012_sigma_selector_source_law_candidate_audit import OUT as P3062

OUT = GEN / "p3063_s2013_signed_source_to_gselector_coupling_theorem_matrix.json"
MD = GEN / "p3063_s2013_signed_source_to_gselector_coupling_theorem_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "explicit_signed_source_to_gselector_coupling": r"explicit coupling theorem|coupling theorem.*G_selector|G_selector.*coupling theorem|signed source.*G_selector",
    "selector_start_lay_analogy_boundary": r"simulation start|start.*simulation|neural network universe|strict kernel.*laws|kernel.*laws of universe|selector.*start",
    "unique_polarity_and_nonpremise_source": r"unique polarity|non-premise.*source law|nonpremise.*source law|polarity.*coupling",
}

SOURCE_ROWS = [
    {"source": "boundary_cocycle_orientation_sign", "signed_value": True, "strict_source_law_exported": False, "nonpremise_origin": False},
    {"source": "chiral_bispectrum_Im_B_1_5_sign", "signed_value": True, "strict_source_law_exported": False, "nonpremise_origin": False},
    {"source": "receiver_winding_sign", "signed_value": True, "strict_source_law_exported": False, "nonpremise_origin": False},
    {"source": "coefficient_sign_convention", "signed_value": True, "strict_source_law_exported": False, "nonpremise_origin": False},
]

COUPLING_POLARITIES = ("source_plus_maps_to_G_plus", "source_plus_maps_to_G_minus")
THEOREM_CRITERIA = (
    "concrete_signed_value",
    "strict_source_law_exported",
    "nonpremise_origin_or_localizer",
    "coupling_map_exists",
    "unique_polarity_selected_by_theory",
    "installed_as_G_selector_source_row",
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def evaluate_pair(source: dict[str, Any], polarity: str) -> dict[str, Any]:
    criteria = {
        "concrete_signed_value": source["signed_value"],
        "strict_source_law_exported": source["strict_source_law_exported"],
        "nonpremise_origin_or_localizer": source["nonpremise_origin"],
        "coupling_map_exists": True,
        "unique_polarity_selected_by_theory": False,
        "installed_as_G_selector_source_row": False,
    }
    missing = [name for name in THEOREM_CRITERIA if not criteria[name]]
    return {
        "source": source["source"],
        "coupling_polarity": polarity,
        "criteria": criteria,
        "missing_criteria": missing,
        "accepted_coupling_theorem": not missing,
        "blocker": "accepted" if not missing else "missing " + ", ".join(missing),
    }


def build_payload() -> dict[str, Any]:
    p3062 = read_json(P3062)
    grep_rows = content_grep()
    theorem_rows = [evaluate_pair(source, polarity) for source, polarity in product(SOURCE_ROWS, COUPLING_POLARITIES)]
    accepted = [row for row in theorem_rows if row["accepted_coupling_theorem"]]
    obligations = [
        {"obligation": "content_first_grep_before_coupling_theorem_matrix", "satisfied": True, "detail": "searched by coupling theorem, lay analogy, and polarity/source-law content"},
        {"obligation": "construct_explicit_coupling_theorem_normal_form", "satisfied": True, "detail": "six criteria separate source value, source provenance, coupling existence, unique polarity, and G_selector installation"},
        {"obligation": "enumerate_current_signed_source_by_coupling_polarity_rows", "satisfied": True, "detail": "four signed current rows times two possible coupling polarities"},
        {"obligation": "export_unique_nonpremise_coupling_theorem", "satisfied": False, "detail": "both coupling polarities remain paired and no strict source provenance installs a row"},
        {"obligation": "export_selector_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, observed physics, L_total, bridge, role transfer, or ToE closure follows"},
    ]
    return {
        "status": "P3063_SIGNED_SOURCE_TO_GSELECTOR_COUPLING_THEOREM_MATRIX_NO_EXPORT",
        "input_hashes": {"P3062": hashlib.sha256(P3062.read_bytes()).hexdigest() if P3062.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "lay_analogy_boundary": {
                "kernel_as_laws": "useful but limited metaphor: K_strict_gate behaves like an operational rule/kernel in the strict pipeline, not yet a licensed complete lawbook of observed physics or L_total",
                "selector_as_start": "useful weak but limited metaphor: a selector would choose an initial branch/direction/representative, but current artifacts do not prove a non-premise start button for a universe simulation or neural network",
                "ontology_correction": "nadsoliton is the primordial fractal information in a solitonic state; there is no deeper neural substrate below it in current guardrails",
            },
            "coupling_theorem_normal_form": {
                "object": "ExplicitSignedSourceToGSelectorCouplingTheoremMatrix",
                "target_missing_theorem": "derive a non-premise unique-polarity coupling from one concrete signed source law/value into G_selector",
                "criteria": list(THEOREM_CRITERIA),
                "acceptance_rule": "accept only if exactly one polarity is selected by a strict source theorem and installed as a G_selector source row",
            },
            "theorem_rows": theorem_rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "signed_source_rows": len(SOURCE_ROWS),
            "coupling_polarities_per_source": len(COUPLING_POLARITIES),
            "coupling_theorem_rows": len(theorem_rows),
            "rows_with_coupling_map_exists": sum(1 for row in theorem_rows if row["criteria"]["coupling_map_exists"]),
            "rows_with_unique_polarity_selected_by_theory": sum(1 for row in theorem_rows if row["criteria"]["unique_polarity_selected_by_theory"]),
            "accepted_coupling_theorems": len(accepted),
            "p3062_status_seen": p3062.get("status"),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3063 constructs the explicit signed-source-to-G_selector coupling theorem normal form requested after P3062.  Four currently signed rows each admit two formal coupling polarities, so the finite matrix has 8 rows.  All 8 rows have a formal coupling map, but 0 rows have a strict source law, non-premise localizer, unique theory-selected polarity, or installed G_selector source row.  Thus the result is a coupling-theorem obstruction, not a selector export.",
            "negative_export_flags": {k: False for k in ["explicit_gselector_coupling_theorem_exported", "strict_nonpremise_sigma_selector_exported", "actual_g_selector_formula_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "simulation_start_exported", "neural_network_universe_start_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not repeat formal two-polarity coupling maps.  The next proof-grade move must either construct a strict polarity-source theorem that selects one coupling polarity and installs it as a G_selector row, or pivot to a different P3057 atom with content-first grep and a finite acceptance matrix.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3063/S2013 signed-source to G_selector coupling theorem matrix", "", f"Status: `{payload['status']}`", "", "## Lay analogy boundary", "- `K_strict_gate` as 'laws of the universe' is a useful **weak metaphor** only: it is a strict-pipeline operational kernel/rule, not yet a proven complete physical lawbook or `L_total`.", "- The selector as 'start of a simulation/neural network of the universe' is also a useful **weak metaphor** only: it points to choosing a branch/direction/representative, not to an exported cosmic start button.", "- Current ontology keeps the nadsoliton itself as primordial fractal information in a solitonic state; no lower neural substrate is added.", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- signed source rows: `{c['signed_source_rows']}`", f"- coupling polarities per source: `{c['coupling_polarities_per_source']}`", f"- coupling theorem rows: `{c['coupling_theorem_rows']}`", f"- rows with formal coupling map: `{c['rows_with_coupling_map_exists']}`", f"- rows with unique theory-selected polarity: `{c['rows_with_unique_polarity_selected_by_theory']}`", f"- accepted coupling theorems: `{c['accepted_coupling_theorems']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3063/S2013 signed-source to `G_selector` coupling theorem matrix", "## P3063/S2013 signed-source to `G_selector` coupling theorem matrix\n\n`P3063/S2013` attacks the exact post-`P3062` gap: an explicit coupling theorem from one concrete signed source law/value into `G_selector`.  It also records the lay analogy boundary: `K_strict_gate` may be compared to an operational rule/kernel, and the selector may be compared to a branch/start choice, but neither comparison licenses observed physics, a universe-simulation start, a neural-network substrate, `QW-2191` discharge, `L_total`, bridge closure, role transfer, or ToE closure.  The finite theorem matrix checks `4` signed source rows times `2` coupling polarities.  All `8` rows have a formal coupling map, but `0` have a strict source law, non-premise localizer, theory-selected unique polarity, or installed `G_selector` source row.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3063/S2013 signed-source to `G_selector` coupling theorem `L_total` guard", "## P3063/S2013 signed-source to `G_selector` coupling theorem `L_total` guard\n\n`P3063/S2013` adds no physical `L_total` term.  It is a finite coupling-theorem obstruction matrix; no row exports a unit-bearing signed action/EOM carrier or nonproxy variational chain rule.\n")
    append_once(AGENTS, "Current signed-source to `G_selector` coupling theorem guardrail (P3063/S2013, 2026-06-23)", "## Current signed-source to `G_selector` coupling theorem guardrail (P3063/S2013, 2026-06-23)\n\n- P3063 constructs the post-P3062 explicit coupling theorem matrix from currently signed source rows into `G_selector`.\n- The matrix has `4 x 2 = 8` rows: four signed source rows and two possible coupling polarities.  All `8` rows have a formal coupling map, but `0` have a strict source law, non-premise localizer, unique theory-selected polarity, or installed `G_selector` source row.\n- The lay analogy that the strict kernel is like laws and the selector is like a simulation/neural-network start is allowed only as weak intuition; it must not add a lower information layer under the nadsoliton or be promoted to observed physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must construct a strict polarity-source theorem selecting one coupling polarity and installing it as a `G_selector` row, or pivot to a different P3057 atom while preserving the P3048-P3063 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
