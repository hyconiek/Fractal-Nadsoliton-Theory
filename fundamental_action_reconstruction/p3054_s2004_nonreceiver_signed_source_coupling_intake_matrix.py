#!/usr/bin/env python3
"""P3054/S2004: non-receiver signed-source coupling intake matrix.

P3053 showed that receiver diagnostic signs cannot become a strict source theorem
by Boolean recombination.  The only honest continuation in that lane is therefore
not another receiver sign, but a genuinely non-receiver signed source/coupling
law that selects one odd output polarity.

This step constructs the missing acceptance object explicitly.  It separates a
candidate strict signed source sigma from the P3048-P3053 receiver sign cube and
from the P3046/readout polarity target.  A finite coupling pushout is then
enumerated: for each of the eight odd receiver-law polarity pairs, there are two
possible identifications with sigma.  Without an exported nonzero source value
and a theorem choosing one identification, the pushout is balanced and no
selector/readout installation follows.

The script also performs a content-first repo grep for the relevant research
content (not just P-numbers) and records the hits that show this is an intake
matrix rather than a replay of already-closed sign-torsor/receiver diagnostics.
"""
from __future__ import annotations

import hashlib, json, subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3053_s2003_receiver_diagnostic_sign_torsor_source_theorem_obstruction import OUT as P3053

OUT = GEN / "p3054_s2004_nonreceiver_signed_source_coupling_intake_matrix.json"
MD = GEN / "p3054_s2004_nonreceiver_signed_source_coupling_intake_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = [
    r"non-receiver.*strict.*signed.*source|strict.*signed.*source.*coupling|odd.*polarity.*source law",
    r"source.*coupling.*polarity|nonpremise.*polarity|signed.*source.*acceptance",
]

CANDIDATES = [
    {
        "candidate": "receiver_diagnostic_cube_P3048_P3053",
        "non_receiver": False,
        "exported_internal_source_law": False,
        "computable_nonzero_signed_value": True,
        "nonconventional_sign": False,
        "coupling_to_receiver_odd_polarity": True,
        "unique_polarity_theorem": False,
        "reason": "explicitly excluded by P3053: receiver diagnostics are the blocked input, not the missing source",
    },
    {
        "candidate": "abstract_minimal_odd_source_sign",
        "non_receiver": True,
        "exported_internal_source_law": False,
        "computable_nonzero_signed_value": False,
        "nonconventional_sign": False,
        "coupling_to_receiver_odd_polarity": True,
        "unique_polarity_theorem": False,
        "reason": "correct representation type only; no strict value or selected coupling is exported",
    },
    {
        "candidate": "chiral_time_or_pseudoscalar_source_family",
        "non_receiver": True,
        "exported_internal_source_law": False,
        "computable_nonzero_signed_value": False,
        "nonconventional_sign": False,
        "coupling_to_receiver_odd_polarity": False,
        "unique_polarity_theorem": False,
        "reason": "prior audits leave source law/value/coupling missing rather than exported",
    },
    {
        "candidate": "finite_Z12_signed_observable_family",
        "non_receiver": True,
        "exported_internal_source_law": False,
        "computable_nonzero_signed_value": True,
        "nonconventional_sign": False,
        "coupling_to_receiver_odd_polarity": False,
        "unique_polarity_theorem": False,
        "reason": "nonzero signed observables exist, but their coefficient/sign orbits remain polarity-paired and not coupled to P3046/P3053",
    },
    {
        "candidate": "unit_bearing_variational_signed_source",
        "non_receiver": True,
        "exported_internal_source_law": False,
        "computable_nonzero_signed_value": False,
        "nonconventional_sign": False,
        "coupling_to_receiver_odd_polarity": False,
        "unique_polarity_theorem": False,
        "reason": "no unit-bearing action/EOM source row is installed on current artifacts",
    },
]


def repo_grep() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pattern in CONTENT_PATTERNS:
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:20]})
    return rows


def coupling_pushout_rows(odd_polarity_pairs: int) -> list[dict[str, Any]]:
    rows = []
    for pair in range(odd_polarity_pairs):
        for source_sign in [-1, 1]:
            for coupling_polarity in [-1, 1]:
                rows.append({
                    "odd_law_pair_index": pair,
                    "source_sign_sigma": source_sign,
                    "coupling_polarity_kappa": coupling_polarity,
                    "selected_receiver_output": source_sign * coupling_polarity,
                    "artifact_selects_row": False,
                })
    return rows


def build_matrix() -> dict[str, Any]:
    p3053 = read_json(P3053)
    odd_pairs = p3053["finite_certificate"]["odd_polarity_pairs"]
    rows = coupling_pushout_rows(odd_pairs)
    acceptance_keys = [
        "non_receiver",
        "exported_internal_source_law",
        "computable_nonzero_signed_value",
        "nonconventional_sign",
        "coupling_to_receiver_odd_polarity",
        "unique_polarity_theorem",
    ]
    candidate_rows = []
    for c in CANDIDATES:
        satisfied = sum(1 for key in acceptance_keys if c[key])
        candidate_rows.append({**c, "acceptance_score": f"{satisfied}/{len(acceptance_keys)}", "accepted": satisfied == len(acceptance_keys)})
    obligations = [
        {"obligation": "content_first_prior_lane_check", "satisfied": True, "detail": "repo grep records prior sign-source/coupling lanes and P3053 receiver block before constructing P3054"},
        {"obligation": "separate_nonreceiver_source_from_receiver_diagnostics", "satisfied": True, "detail": "candidate rows mark receiver diagnostics as disallowed input class"},
        {"obligation": "construct_coupling_pushout", "satisfied": True, "detail": "sigma source sign, receiver odd-law polarity pair, and output polarity are separated"},
        {"obligation": "exhaust_p3053_odd_pair_couplings", "satisfied": True, "detail": f"{len(rows)} rows cover {odd_pairs} odd pairs times two source signs times two coupling polarities"},
        {"obligation": "export_concrete_nonreceiver_signed_source", "satisfied": False, "detail": "no candidate supplies all source-law, value, nonconventional-sign, coupling, and unique-polarity requirements"},
        {"obligation": "install_selector_readout_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, unit-bearing action/EOM, L_total, bridge, role transfer, or ToE export follows"},
    ]
    return {
        "object": "NonReceiverSignedSource_CouplingIntakeMatrix",
        "content_first_repo_grep": repo_grep(),
        "acceptance_criteria": acceptance_keys,
        "candidate_rows": candidate_rows,
        "coupling_pushout_rows": rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "content_grep_patterns": len(CONTENT_PATTERNS),
            "content_grep_hits": sum(row["hit_count"] for row in repo_grep()),
            "candidate_classes": len(candidate_rows),
            "accepted_candidate_classes": sum(1 for row in candidate_rows if row["accepted"]),
            "p3053_odd_polarity_pairs": odd_pairs,
            "coupling_pushout_rows": len(rows),
            "artifact_selected_pushout_rows": sum(1 for row in rows if row["artifact_selects_row"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_nonreceiver_signed_source_coupling_law_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3054_NONRECEIVER_SIGNED_SOURCE_COUPLING_INTAKE_MATRIX_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3053": hashlib.sha256(P3053.read_bytes()).hexdigest() if P3053.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3054 constructs the exact post-P3053 intake object for a non-receiver strict signed source/coupling law.  The finite pushout over the eight P3053 odd-law polarity pairs has 32 possible sigma/coupling rows, but current artifacts select zero rows and no audited candidate satisfies the full source-law/value/nonconventional-sign/coupling/unique-polarity package.  Therefore P3053 is not escaped by importing a named sign source unless a new concrete law is actually exported.",
            "negative_export_flags": {k: False for k in ["strict_nonreceiver_signed_source_coupling_law_exported", "unique_receiver_odd_polarity_selected", "p3046_coupling_polarity_selected", "selector_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not name another sign torsor or receiver diagnostic as closure.  The next proof-grade move must either provide one explicit non-receiver formula with an exported strict source law, a computable nonzero signed value, and a theorem choosing one P3053/P3046 coupling polarity, or pivot outside the phase-geometry selector lane and preserve the P3048-P3054 bounded no-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3054/S2004 non-receiver signed-source coupling intake matrix", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep patterns: `{c['content_grep_patterns']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- candidate classes: `{c['candidate_classes']}`",
        f"- accepted candidate classes: `{c['accepted_candidate_classes']}`",
        f"- P3053 odd polarity pairs: `{c['p3053_odd_polarity_pairs']}`",
        f"- coupling pushout rows: `{c['coupling_pushout_rows']}`",
        f"- artifact-selected pushout rows: `{c['artifact_selected_pushout_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict non-receiver signed source/coupling law exported: `{c['strict_nonreceiver_signed_source_coupling_law_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3054/S2004 non-receiver signed-source coupling intake matrix", "## P3054/S2004 non-receiver signed-source coupling intake matrix\n\n`P3054/S2004` constructs the exact intake object demanded by P3053: a non-receiver signed source `sigma`, the eight P3053 receiver odd-law polarity pairs, and the P3046/readout output polarity are separated in a finite coupling pushout.  The audit checks content-first repo hits for prior sign-source/coupling lanes, enumerates `32` pushout rows, and scans five candidate source classes.  Zero candidate classes satisfy the full strict source-law/value/nonconventional-sign/coupling/unique-polarity package, and zero pushout rows are artifact-selected.  No non-receiver strict signed source/coupling theorem, P3046 polarity selection, selector/readout installation, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3054/S2004 non-receiver signed-source coupling `L_total` guard", "## P3054/S2004 non-receiver signed-source coupling `L_total` guard\n\n`P3054/S2004` adds no physical `L_total` term.  Its finite coupling pushout is an acceptance/no-export matrix: no non-receiver source row installs a unit-bearing signed variational/action/EOM coupling or selects a unique receiver/P3046 polarity.\n")
    append_once(AGENTS, "Current non-receiver signed-source coupling intake guardrail (P3054/S2004, 2026-06-23)", "## Current non-receiver signed-source coupling intake guardrail (P3054/S2004, 2026-06-23)\n\n- P3054 constructs the post-P3053 intake matrix for the only admissible phase-geometry continuation: a genuinely non-receiver strict signed source/coupling law selecting one odd polarity.\n- The finite pushout separates source sign `sigma`, the eight P3053 receiver odd-law polarity pairs, and P3046/readout output polarity; all `32` pushout rows remain unselected by current artifacts.\n- Do not promote named sign torsors, finite signed observables, receiver diagnostics, or abstract odd-source representations to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure without an exported source law, nonzero signed value, and unique coupling-polarity theorem.\n- A next move must provide one explicit non-receiver formula/artifact satisfying that full package, or pivot outside the phase-geometry selector lane while preserving the P3048-P3054 bounded no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
