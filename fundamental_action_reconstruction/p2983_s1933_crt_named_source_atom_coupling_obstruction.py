#!/usr/bin/env python3
"""P2983/S1933: CRT named source-atom coupling obstruction.

P2982 left two concrete CRT routes. This audit attacks exactly one of them:
coupling the Z12 CRT idempotent projector split to named missing source atoms.
It does not replay CRT provenance, nonpremise factor semantics, nilradical lanes,
selector closure, generic bridge completion, role transfer, or formal L_total
promotion.

The finite calculation treats projectors 4 and 9 as exact algebraic factor
anchors and tests whether they supply a coupling witness for four named atoms:
selector/orientation sign, positive damping beta/Z_beta, legacy-to-strict bridge
source, and action-density/variational source. The projector split is exact, but
it remains orientation-blind, scale/measure-unsourced, bridge-law silent, and
action-measure silent.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2982_s1932_crt_factor_semantics_obstruction import OUT as P2982

OUT = GEN / "p2983_s1933_crt_named_source_atom_coupling_obstruction.json"
MD = GEN / "p2983_s1933_crt_named_source_atom_coupling_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MODULUS = 12
PROJECTORS = (4, 9)
NAMED_SOURCE_ATOMS = [
    "selector_orientation_sign_atom",
    "positive_beta_Z_beta_damping_atom",
    "legacy_to_strict_bridge_source_atom",
    "action_density_variational_source_atom",
]


def crt_projector_coupling_witness() -> dict[str, object]:
    projector_rows = []
    for projector in PROJECTORS:
        image = sorted({(projector * x) % MODULUS for x in range(MODULUS)})
        kernel = [x for x in range(MODULUS) if (projector * x) % MODULUS == 0]
        projector_rows.append({
            "projector": projector,
            "residue_signature": {"mod_3": projector % 3, "mod_4": projector % 4},
            "image": image,
            "image_size": len(image),
            "kernel": kernel,
            "kernel_size": len(kernel),
            "idempotent": (projector * projector - projector) % MODULUS == 0,
        })
    signed_orientation_scores = {
        "+omega": {str(p): len({(p * x) % MODULUS for x in range(MODULUS)}) for p in PROJECTORS},
        "-omega": {str(p): len({(p * x) % MODULUS for x in range(MODULUS)}) for p in PROJECTORS},
    }
    ratio_samples = {str(p): {"p_over_12": p / MODULUS, "image_size_over_12": len({(p * x) % MODULUS for x in range(MODULUS)}) / MODULUS} for p in PROJECTORS}
    return {
        "modulus": MODULUS,
        "projectors": list(PROJECTORS),
        "orthogonal_completion_pair": {"a": 4, "b": 9, "product_mod_12": (4 * 9) % MODULUS, "sum_mod_12": (4 + 9) % MODULUS},
        "projector_rows": projector_rows,
        "signed_orientation_scores": signed_orientation_scores,
        "orientation_score_gap": {str(p): signed_orientation_scores["+omega"][str(p)] - signed_orientation_scores["-omega"][str(p)] for p in PROJECTORS},
        "ratio_samples": ratio_samples,
        "exports_positive_beta_scale": False,
        "exports_bridge_completion_law": False,
        "exports_action_measure": False,
    }


def source_atom_rows(witness: dict[str, object]) -> list[dict[str, object]]:
    return [
        {
            "atom": "selector_orientation_sign_atom",
            "finite_crt_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "nonpremise_factor_semantics_available": False,
            "orientation_sensitive_coupling": any(gap != 0 for gap in witness["orientation_score_gap"].values()),
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "projector image-size scores are identical for +omega and -omega",
        },
        {
            "atom": "positive_beta_Z_beta_damping_atom",
            "finite_crt_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "nonpremise_factor_semantics_available": False,
            "orientation_sensitive_coupling": False,
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "projector ratios 4/12 and 9/12 are algebraic samples, not target-independent beta/Z_beta or measure sources",
        },
        {
            "atom": "legacy_to_strict_bridge_source_atom",
            "finite_crt_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "nonpremise_factor_semantics_available": False,
            "orientation_sensitive_coupling": False,
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "CRT factor split does not provide amplitude, phase/topological, or nonlinear damping completion map",
        },
        {
            "atom": "action_density_variational_source_atom",
            "finite_crt_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "nonpremise_factor_semantics_available": False,
            "orientation_sensitive_coupling": False,
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "no unit-bearing density, measure, field variable, or nonproxy variational chain is attached",
        },
        {
            "atom": "completed_strict_CRT_named_atom_coupling_schema",
            "finite_crt_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": True,
            "nonpremise_factor_semantics_available": True,
            "orientation_sensitive_coupling": True,
            "positive_scale_or_measure_source": True,
            "bridge_completion_law": True,
            "action_variational_installation": True,
            "accepted_current_coupling": False,
            "obstruction": "schema row only; no current artifact exports the required theorem package",
        },
    ]


def coupling_obligations(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    current = [r for r in rows if r["atom"] != "completed_strict_CRT_named_atom_coupling_schema"]
    return [
        {"obligation": "finite_CRT_projector_witness", "satisfied": all(r["finite_crt_witness"] for r in current), "evidence": "projectors 4 and 9 satisfy 4*9=0 and 4+9=1 mod 12"},
        {"obligation": "exactly_named_source_atom_targeted", "satisfied": all(r["named_atom_targeted"] for r in current) and len(current) == len(NAMED_SOURCE_ATOMS), "evidence": "selector/orientation, beta/Z_beta, bridge-source, and action-density atoms are enumerated explicitly"},
        {"obligation": "strict_provenance_available", "satisfied": any(r["strict_provenance_available"] for r in current), "evidence": "P2981 found no strict CRT provenance"},
        {"obligation": "nonpremise_factor_semantics_available", "satisfied": any(r["nonpremise_factor_semantics_available"] for r in current), "evidence": "P2982 found only algebraic factor labels"},
        {"obligation": "orientation_sensitive_coupling", "satisfied": any(r["orientation_sensitive_coupling"] for r in current), "evidence": "CRT projector scores do not distinguish +omega from -omega"},
        {"obligation": "positive_scale_or_measure_source", "satisfied": any(r["positive_scale_or_measure_source"] for r in current), "evidence": "formal projector ratios are not target-independent beta/Z_beta or measure sources"},
        {"obligation": "bridge_completion_law", "satisfied": any(r["bridge_completion_law"] for r in current), "evidence": "no amplitude/phase/damping completion map is supplied"},
        {"obligation": "action_variational_installation", "satisfied": any(r["action_variational_installation"] for r in current), "evidence": "no unit-bearing named density or variational chain is installed"},
        {"obligation": "accepted_current_coupling", "satisfied": any(r["accepted_current_coupling"] for r in current), "evidence": "all current rows fail at least one strict coupling premise"},
    ]


def acceptance_matrix() -> list[dict[str, object]]:
    names = ["finite_CRT_witness", "named_atom", "strict_provenance", "nonpremise_factor_semantics", "orientation_sensitive", "positive_scale_measure", "bridge_law", "action_variational"]
    return [{"present": dict(zip(names, bits)), "accepts_CRT_named_source_atom_coupling": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2982_path: object) -> dict[str, object]:
    witness = crt_projector_coupling_witness()
    rows = source_atom_rows(witness)
    obligations = coupling_obligations(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2983_CRT_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2982": hashlib.sha256(p2982_path.read_bytes()).hexdigest() if p2982_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "CRTNamedSourceAtomCoupling_ObstructionMatrix",
            "crt_projector_coupling_witness": witness,
            "named_source_atom_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "coupling_certificate": {
            "named_source_atom_count": len(NAMED_SOURCE_ATOMS),
            "candidate_row_count": len(rows),
            "projectors": list(PROJECTORS),
            "orthogonal_completion_pair": witness["orthogonal_completion_pair"],
            "orientation_score_gap": witness["orientation_score_gap"],
            "accepted_current_couplings": [r["atom"] for r in rows if r["accepted_current_coupling"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_CRT_named_source_atom_coupling"]),
        },
        "decision": {
            "positive_progress": "P2983 attacks exactly one remaining CRT route after P2982: named source-atom coupling for the CRT idempotent projector split.",
            "breakthrough": "Bounded no-go: projectors 4 and 9 form an exact finite CRT split, but they do not export strict provenance, nonpremise factor semantics, orientation-sensitive coupling, target-independent positive beta/Z_beta or measure, a bridge-completion law, or a unit-bearing action/variational installation.",
            "negative_export_flags": {k: False for k in ["named_source_atom_coupling_exported", "CRT_strict_provenance_exported", "nonpremise_factor_semantics_exported", "selector_or_orientation_exported", "damping_source_exported", "bridge_closure_exported", "action_density_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay CRT source-atom rows, projector ratios, residue signatures, nilradical lanes, selector replay, bridge maps, role transfer, or L_total placeholders.  The only remaining CRT-lane route is action/variational installation with a genuinely unit-bearing named density and nonproxy chain; otherwise introduce a genuinely new strict typed object/theorem/provider or preserve the P2929-P2983 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, object]) -> None:
    cert = payload["coupling_certificate"]
    lines = [
        "# P2983/S1933 CRT named source-atom coupling obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling certificate",
        f"- named source atoms: `{cert['named_source_atom_count']}`",
        f"- candidate rows: `{cert['candidate_row_count']}`",
        f"- projectors: `{cert['projectors']}`",
        f"- orthogonal completion pair: `{cert['orthogonal_completion_pair']}`",
        f"- orientation score gap: `{cert['orientation_score_gap']}`",
        f"- accepted current couplings: `{cert['accepted_current_couplings']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
        "",
        "## Lay summary",
        payload["decision"]["positive_progress"],
        payload["decision"]["breakthrough"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, object]:
    read_json(P2982)
    payload = build_payload(P2982)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2983/S1933 CRT named source-atom coupling obstruction", "## P2983/S1933 CRT named source-atom coupling obstruction\n\n`P2983/S1933` attacks exactly one P2982 remaining CRT route: named source-atom coupling for the CRT idempotent projector split.  The finite matrix audits selector/orientation sign, positive `beta/Z_beta` damping, legacy-to-strict bridge-source, and action-density/variational atoms.  The result is bounded no-go: projectors `4` and `9` satisfy `4*9=0` and `4+9=1 mod 12`, but their orientation-score gaps are `0`, their ratios are only algebraic samples, no amplitude/phase/damping completion law is supplied, and no unit-bearing named action density or nonproxy variational chain is installed.  No named source-atom coupling, damping source, selector/orientation source, bridge closure, action-density export, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2983/S1933 CRT source-atom coupling `L_total` guard", "## P2983/S1933 CRT source-atom coupling `L_total` guard\n\n`P2983/S1933` does not add a CRT source-atom coupling term to `L_total`.  The named source-atom audit finds no strict provenance, nonpremise factor semantics, positive scale/measure source, bridge law, unit-bearing action density, or nonproxy variational chain, so EOM, Hamiltonian, role transfer, bridge closure, and ToE remain unexported.\n")
    append_once(AGENTS, "Current CRT named source-atom coupling obstruction guardrail (P2983/S1933, 2026-06-20)", "## Current CRT named source-atom coupling obstruction guardrail (P2983/S1933, 2026-06-20)\n\n- P2983 attacks exactly one P2982 remaining CRT route: named source-atom coupling for the `Z/12Z` CRT idempotent projector split.\n- Finite positives remain real: projectors `4` and `9` satisfy `4*9=0` and `4+9=1 mod 12`; the four named atoms are explicitly audited; the acceptance matrix has `256` profiles with only the full profile accepting.\n- The current route is bounded no-go: no strict provenance, nonpremise factor semantics, orientation-sensitive coupling, target-independent positive `beta/Z_beta` or measure, bridge-completion law, or unit-bearing action/variational installation is exported.\n- Do not promote CRT source-atom rows, projector ratios, residue signatures, nilradical replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  The only remaining CRT-lane route is action/variational installation with a genuinely unit-bearing named density and nonproxy chain; otherwise introduce a genuinely new strict typed object/theorem/provider or preserve the P2929-P2983 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
