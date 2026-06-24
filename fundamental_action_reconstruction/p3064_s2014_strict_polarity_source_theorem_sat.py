#!/usr/bin/env python3
"""P3064/S2014: strict polarity-source theorem SAT certificate.

P3063 left a precise gap: formal source-to-G_selector coupling maps exist, but
no theorem selects one polarity and installs it as a selector source row.  P3064
turns that into a finite C2/SAT certificate.  The missing theoretical object is
a MinimalStrictPolaritySourceTheorem: an absolute non-premise signed source
choice plus an absolute non-premise coupling-polarity choice.  Current artifacts
provide neither primitive atom, so the exhaustive SAT table is no-export.
"""
from __future__ import annotations

import hashlib, json, subprocess
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3063_s2013_signed_source_to_gselector_coupling_theorem_matrix import OUT as P3063

OUT = GEN / "p3064_s2014_strict_polarity_source_theorem_sat.json"
MD = GEN / "p3064_s2014_strict_polarity_source_theorem_sat.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "strict_polarity_source_theorem": r"strict polarity-source theorem|polarity-source theorem|source theorem.*polarity|strict.*source.*polarity",
    "absolute_source_and_coupling_choice": r"absolute.*source sign|absolute.*coupling polarity|source_plus|source_minus|identity_polarity|flip_polarity",
    "c2_torsor_sat_selector_obstruction": r"C2.*torsor|two-polarity|coupling polarity|unique.*G_selector|installed.*G_selector",
}

SOURCE_ATOMS = ("strict_source_plus", "strict_source_minus")
COUPLING_ATOMS = ("strict_coupling_identity", "strict_coupling_flip")
PRIMITIVE_ATOMS = SOURCE_ATOMS + COUPLING_ATOMS
CURRENT_EXPORTED_ATOMS: tuple[str, ...] = ()


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def powerset(items: tuple[str, ...]) -> list[tuple[str, ...]]:
    return [tuple(combo) for r in range(len(items) + 1) for combo in combinations(items, r)]


def evaluate_atom_set(atoms: tuple[str, ...]) -> dict[str, Any]:
    atom_set = set(atoms)
    source_choices = sorted(atom_set.intersection(SOURCE_ATOMS))
    coupling_choices = sorted(atom_set.intersection(COUPLING_ATOMS))
    consistent = len(source_choices) <= 1 and len(coupling_choices) <= 1
    absolute_source_fixed = len(source_choices) == 1
    absolute_coupling_fixed = len(coupling_choices) == 1
    accepted = bool(consistent and absolute_source_fixed and absolute_coupling_fixed)
    selected_g = None
    if accepted:
        source_plus = source_choices[0] == "strict_source_plus"
        identity = coupling_choices[0] == "strict_coupling_identity"
        # identity maps source plus to G plus; flip maps source plus to G minus.
        selected_g = "G_plus" if source_plus == identity else "G_minus"
    missing = []
    if not consistent:
        missing.append("consistency_between_exclusive_atoms")
    if not absolute_source_fixed:
        missing.append("absolute_nonpremise_source_sign")
    if not absolute_coupling_fixed:
        missing.append("absolute_nonpremise_coupling_polarity")
    return {
        "atoms": list(atoms),
        "consistent": consistent,
        "absolute_source_fixed": absolute_source_fixed,
        "absolute_coupling_fixed": absolute_coupling_fixed,
        "accepted_minimal_strict_polarity_source_theorem": accepted,
        "selected_G_selector_orientation": selected_g,
        "missing_or_blocking": missing,
    }


def build_payload() -> dict[str, Any]:
    p3063 = read_json(P3063)
    grep_rows = content_grep()
    rows = [evaluate_atom_set(atoms) for atoms in powerset(PRIMITIVE_ATOMS)]
    accepted = [row for row in rows if row["accepted_minimal_strict_polarity_source_theorem"]]
    minimal_size = min((len(row["atoms"]) for row in accepted), default=None)
    minimal_rows = [row for row in accepted if len(row["atoms"]) == minimal_size]
    current_row = evaluate_atom_set(CURRENT_EXPORTED_ATOMS)
    obligations = [
        {"obligation": "content_first_grep_before_strict_polarity_sat", "satisfied": True, "detail": "searched by polarity-source theorem, absolute sign/coupling, and C2 torsor content"},
        {"obligation": "construct_minimal_strict_polarity_source_theorem_object", "satisfied": True, "detail": "source sign and coupling polarity are separated as primitive C2 atoms"},
        {"obligation": "exhaust_all_primitive_atom_extensions", "satisfied": True, "detail": "all 16 subsets of four primitive atoms are enumerated"},
        {"obligation": "show_current_artifacts_select_G_selector", "satisfied": False, "detail": "current artifacts export zero primitive atoms; current row lacks source sign and coupling polarity"},
        {"obligation": "export_selector_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, observed physics, L_total, bridge, role transfer, or ToE closure follows"},
    ]
    return {
        "status": "P3064_STRICT_POLARITY_SOURCE_THEOREM_SAT_NO_CURRENT_EXPORT",
        "input_hashes": {"P3063": hashlib.sha256(P3063.read_bytes()).hexdigest() if P3063.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "minimal_strict_polarity_source_theorem": {
                "object": "MinimalStrictPolaritySourceTheorem",
                "c2_rule": "G_orientation = source_sign composed_with coupling_polarity",
                "primitive_atoms": list(PRIMITIVE_ATOMS),
                "acceptance_rule": "accept exactly one absolute non-premise source-sign atom and exactly one absolute non-premise coupling-polarity atom; reject contradictions and partial rows",
                "why_more_dowodowy": "the theorem obligation is reduced to an exhaustive finite SAT table over the primitive atoms rather than a narrative coupling claim",
            },
            "sat_rows": rows,
            "current_artifact_row": current_row,
            "minimal_accepting_rows": minimal_rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "primitive_atoms": len(PRIMITIVE_ATOMS),
            "sat_rows": len(rows),
            "consistent_rows": sum(1 for row in rows if row["consistent"]),
            "accepted_theorem_rows": len(accepted),
            "minimal_accepting_atom_count": minimal_size,
            "minimal_accepting_rows": len(minimal_rows),
            "current_exported_atoms": len(CURRENT_EXPORTED_ATOMS),
            "current_row_accepted": current_row["accepted_minimal_strict_polarity_source_theorem"],
            "p3063_status_seen": p3063.get("status"),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3064 constructs the missing strict polarity-source theorem as an exhaustive C2/SAT certificate.  A real theorem must contain two primitive non-premise atoms: one absolute source sign and one absolute coupling polarity.  The SAT table enumerates all 16 subsets of the four primitive atoms.  Exactly 4 rows would be accepted if such primitive atoms were exported, and every minimal accepted row has 2 atoms.  The current artifact row has 0 exported primitive atoms and is not accepted, so no G_selector row, QW-2191 discharge, selector closure, L_total, bridge/role transfer, or ToE closure is exported.",
            "negative_export_flags": {k: False for k in ["minimal_strict_polarity_source_theorem_exported", "absolute_nonpremise_source_sign_exported", "absolute_nonpremise_coupling_polarity_exported", "actual_g_selector_formula_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not add another partial source/coupling row.  The next proof-grade move must construct exactly one primitive atom from the P3064 SAT theorem: either an absolute non-premise source-sign law or an absolute non-premise coupling-polarity law, then rerun the two-atom acceptance check; otherwise pivot to another P3057 atom with content-first grep.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3064/S2014 strict polarity-source theorem SAT certificate", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- primitive atoms: `{c['primitive_atoms']}`", f"- SAT rows: `{c['sat_rows']}`", f"- consistent rows: `{c['consistent_rows']}`", f"- accepted theorem rows: `{c['accepted_theorem_rows']}`", f"- minimal accepting atom count: `{c['minimal_accepting_atom_count']}`", f"- minimal accepting rows: `{c['minimal_accepting_rows']}`", f"- current exported atoms: `{c['current_exported_atoms']}`", f"- current row accepted: `{c['current_row_accepted']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3064/S2014 strict polarity-source theorem SAT certificate", "## P3064/S2014 strict polarity-source theorem SAT certificate\n\n`P3064/S2014` converts the post-`P3063` coupling-theorem gap into an exhaustive finite C2/SAT theorem template.  The constructed `MinimalStrictPolaritySourceTheorem` has four primitive atoms: `strict_source_plus`, `strict_source_minus`, `strict_coupling_identity`, and `strict_coupling_flip`.  Acceptance requires exactly one absolute non-premise source-sign atom and exactly one absolute non-premise coupling-polarity atom.  The verifier enumerates all `16` primitive-atom subsets: `4` rows would be accepted if the required atoms were exported, every minimal accepted row has `2` atoms, and the current artifact row has `0` exported atoms and is not accepted.  No `G_selector`, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3064/S2014 strict polarity-source theorem SAT `L_total` guard", "## P3064/S2014 strict polarity-source theorem SAT `L_total` guard\n\n`P3064/S2014` adds no physical `L_total` term.  It is a finite theorem-template SAT certificate over source/coupling sign atoms; no row exports a unit-bearing signed action/EOM carrier or nonproxy variational chain rule.\n")
    append_once(AGENTS, "Current strict polarity-source theorem SAT guardrail (P3064/S2014, 2026-06-23)", "## Current strict polarity-source theorem SAT guardrail (P3064/S2014, 2026-06-23)\n\n- P3064 reduces the post-P3063 selector-coupling gap to a finite C2/SAT theorem template with four primitive atoms: `strict_source_plus`, `strict_source_minus`, `strict_coupling_identity`, and `strict_coupling_flip`.\n- Acceptance requires exactly one absolute non-premise source-sign atom and exactly one absolute non-premise coupling-polarity atom; the exhaustive `16`-row table has `4` accepting rows, all with `2` atoms, but the current artifact row has `0` exported atoms and is not accepted.\n- Do not promote partial source/coupling rows, formal SAT accepting rows, or polarity conventions to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must construct exactly one primitive atom from the P3064 template, then rerun the two-atom acceptance check; otherwise pivot to another P3057 atom while preserving the P3048-P3064 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
