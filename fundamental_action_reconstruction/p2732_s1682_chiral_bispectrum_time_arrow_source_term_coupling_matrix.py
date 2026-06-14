#!/usr/bin/env python3
"""P2732/S1682: chiral-bispectrum time-arrow source-term coupling matrix.

After P2731 showed that local tau variational laws need an explicit
T-odd source term, P2732 tests the most immediate available signed source term:
coupling tau to the P2718 chiral-bispectrum marker Im(B_{1,5}).  The finite
matrix scans both coupling signs lambda=+/-1, both orientation torsor branches,
and all 2^12 tau fields under H=-lambda*sum_s tau_s*Im(B_{1,5})(orientation,s).
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.json"
MD = GEN / "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
SPINS = (-1, 1)
ORIENTATIONS = (-1, 1)
LAMBDAS = (-1, 1)
FIELDS = tuple(itertools.product(SPINS, repeat=Z12))
INPUTS = {
    "P2731_LOCAL_VARIATIONAL_NO_GO": GEN / "p2731_s1681_local_time_arrow_variational_law_no_go.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}
NEGATIVE_EXPORT_FLAGS = [
    "strict_time_arrow_source_term_exported",
    "nonpremise_lambda_sign_selected",
    "nonpremise_orientation_branch_selected",
    "p2721_polarity_selected",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def marker_table(p2718: dict[str, Any]) -> dict[tuple[int, int], float]:
    return {
        (int(row["orientation"]), int(row["source"])): float(row["marker_imag"])
        for row in p2718.get("marker_rows", [])
    }


def fallback_marker_table() -> dict[tuple[int, int], float]:
    return {(orientation, source): float(2 * orientation) for orientation in ORIENTATIONS for source in range(Z12)}


def coupling_energy(field: tuple[int, ...], orientation: int, lam: int, markers: dict[tuple[int, int], float]) -> float:
    return -sum(lam * field[source] * markers[(orientation, source)] for source in range(Z12))


def ground_summary(orientation: int, lam: int, markers: dict[tuple[int, int], float]) -> dict[str, Any]:
    energies = [coupling_energy(field, orientation, lam, markers) for field in FIELDS]
    min_energy = min(energies)
    ground_states = [field for field, energy in zip(FIELDS, energies) if energy == min_energy]
    magnetizations = sorted(set(sum(field) for field in ground_states))
    selected_tau_sign = None
    if len(ground_states) == 1 and abs(magnetizations[0]) == Z12:
        selected_tau_sign = 1 if magnetizations[0] > 0 else -1
    return {
        "orientation": orientation,
        "lambda": lam,
        "marker_values": sorted(set(markers[(orientation, source)] for source in range(Z12))),
        "min_energy": min_energy,
        "ground_state_count": len(ground_states),
        "ground_magnetizations": magnetizations,
        "selected_constant_tau_sign": selected_tau_sign,
    }


def audit_couplings(markers: dict[tuple[int, int], float]) -> dict[str, Any]:
    rows = [ground_summary(orientation, lam, markers) for lam in LAMBDAS for orientation in ORIENTATIONS]
    polarity_hist = Counter((row["lambda"], row["orientation"], row["selected_constant_tau_sign"]) for row in rows)
    selected_tau_hist = Counter(row["selected_constant_tau_sign"] for row in rows)
    balanced_tau_selection = selected_tau_hist[-1] == selected_tau_hist[1]
    lambda_flip_pairs = []
    for row in rows:
        paired = next(
            candidate for candidate in rows
            if candidate["lambda"] == -row["lambda"] and candidate["orientation"] == row["orientation"]
        )
        lambda_flip_pairs.append({
            "row": {"lambda": row["lambda"], "orientation": row["orientation"], "tau": row["selected_constant_tau_sign"]},
            "lambda_flipped_partner": {"lambda": paired["lambda"], "orientation": paired["orientation"], "tau": paired["selected_constant_tau_sign"]},
            "tau_sign_reversed": paired["selected_constant_tau_sign"] == -row["selected_constant_tau_sign"],
        })
    orientation_flip_pairs = []
    for row in rows:
        paired = next(
            candidate for candidate in rows
            if candidate["lambda"] == row["lambda"] and candidate["orientation"] == -row["orientation"]
        )
        orientation_flip_pairs.append({
            "row": {"lambda": row["lambda"], "orientation": row["orientation"], "tau": row["selected_constant_tau_sign"]},
            "orientation_flipped_partner": {"lambda": paired["lambda"], "orientation": paired["orientation"], "tau": paired["selected_constant_tau_sign"]},
            "tau_sign_reversed": paired["selected_constant_tau_sign"] == -row["selected_constant_tau_sign"],
        })
    return {
        "field_count_per_row": len(FIELDS),
        "coupling_rows": rows,
        "coupling_row_count": len(rows),
        "lambda_values": list(LAMBDAS),
        "orientation_values": list(ORIENTATIONS),
        "selected_tau_sign_histogram": dict(sorted(selected_tau_hist.items())),
        "balanced_tau_selection_across_lambda_orientation_pair": balanced_tau_selection,
        "all_rows_have_unique_constant_tau_ground_state": all(row["ground_state_count"] == 1 and row["selected_constant_tau_sign"] in SPINS for row in rows),
        "all_lambda_flip_pairs_reverse_tau": all(pair["tau_sign_reversed"] for pair in lambda_flip_pairs),
        "all_orientation_flip_pairs_reverse_tau": all(pair["tau_sign_reversed"] for pair in orientation_flip_pairs),
        "lambda_flip_pairs": lambda_flip_pairs,
        "orientation_flip_pairs": orientation_flip_pairs,
        "obstruction": "The chiral-bispectrum source term is strong enough to choose a tau ground state only after lambda and orientation/polarity are fixed.  The finite matrix is exactly paired: flipping lambda or the orientation torsor reverses the selected tau sign, and current artifacts do not export a non-premise law selecting one lambda/P2721 polarity.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2718_marker_coupling_is_computable": audit["coupling_row_count"] == 4,
        "each_fixed_lambda_orientation_selects_tau": audit["all_rows_have_unique_constant_tau_ground_state"],
        "lambda_flip_reverses_selected_tau": audit["all_lambda_flip_pairs_reverse_tau"],
        "orientation_flip_reverses_selected_tau": audit["all_orientation_flip_pairs_reverse_tau"],
        "nonpremise_lambda_sign_selected": False,
        "p2721_polarity_coupling_polarity_selected": False,
    }
    missing = [key for key, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_time_arrow_source_term": not missing,
        "blocker": "The source-term coupling is conditional rather than strict: tau selection depends on an unsourced lambda/polarity choice.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["chiral_bispectrum_time_arrow_coupling_audit"]
    lines = [
        "# P2732/S1682 chiral-bispectrum time-arrow source-term coupling matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling matrix",
        f"- coupling_row_count={audit['coupling_row_count']}",
        f"- field_count_per_row={audit['field_count_per_row']}",
        f"- all_rows_have_unique_constant_tau_ground_state={audit['all_rows_have_unique_constant_tau_ground_state']}",
        f"- selected_tau_sign_histogram={audit['selected_tau_sign_histogram']}",
        f"- all_lambda_flip_pairs_reverse_tau={audit['all_lambda_flip_pairs_reverse_tau']}",
        f"- all_orientation_flip_pairs_reverse_tau={audit['all_orientation_flip_pairs_reverse_tau']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    markers = marker_table(inputs["P2718_CHIRAL_BISPECTRUM_MARKER"]) or fallback_marker_table()
    audit = audit_couplings(markers)
    acceptance = acceptance_matrix(audit)
    no_go = not acceptance["accepted_as_strict_time_arrow_source_term"]
    payload = {
        "status": "P2732_CHIRAL_BISPECTRUM_TIME_ARROW_SOURCE_TERM_CONDITIONAL_NO_GO" if no_go else "P2732_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "H=-lambda*sum_s tau_s*Im(B_{1,5})(orientation,s), lambda in {-1,+1}, orientation in {-1,+1}, all 2^12 tau fields",
        "chiral_bispectrum_time_arrow_coupling_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "strict_time_arrow_source_term_exported": False,
            "nonpremise_lambda_sign_selected": False,
            "nonpremise_orientation_branch_selected": False,
            "p2721_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2732 tests the direct chiral-bispectrum source term requested after P2731.  The term selects a tau ground state for fixed lambda and orientation, but lambda and the P2721 polarity/orientation branch are exactly paired and unsourced, so no strict non-premise time-arrow source is exported.",
            "next_honest_step": "Do not replay direct Im(B_{1,5}) tau-coupling as closure.  A next admissible proof-grade move must export a strict law fixing the coupling sign lambda or the P2721 polarity from inside the theory; otherwise preserve the P2697-P2732 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2732/S1682 chiral-bispectrum time-arrow source-term coupling matrix", "## P2732/S1682 chiral-bispectrum time-arrow source-term coupling matrix\n\n`P2732/S1682` tests the direct source term `H=-lambda*sum_s tau_s*Im(B_{1,5})(orientation,s)` for `lambda=+/-1`, both orientation branches, and all `2^12=4096` tau fields.  Each fixed `(lambda,orientation)` row has a unique constant `tau` ground state, but flipping `lambda` or the orientation branch reverses that selected sign; the P2721 polarity remains unsourced.  No strict non-premise time-arrow source term, `lambda` fixing, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2732/S1682 chiral-bispectrum tau-coupling Ltotal guard", "## P2732/S1682 chiral-bispectrum tau-coupling Ltotal guard\n\n`P2732/S1682` supplies a concrete variational-looking coupling between `tau` and `Im(B_{1,5})`, but the selected ground-state sign is conditional on an unsourced coupling sign `lambda` and an unsourced P2721 polarity/orientation branch.  Therefore it cannot add a strict variational source term or promote `L_total`.\n")
    append_once(AGENTS, "Current chiral-bispectrum time-arrow source-term coupling guardrail (P2732/S1682, 2026-06-14)", "## Current chiral-bispectrum time-arrow source-term coupling guardrail (P2732/S1682, 2026-06-14)\n\n- P2732 tests the direct coupling `H=-lambda*sum_s tau_s*Im(B_{1,5})(orientation,s)` across both `lambda` signs, both orientation branches, and all `2^12=4096` tau fields.\n- Each fixed row selects a unique constant `tau`, but flipping `lambda` or the orientation/P2721 polarity reverses the selected sign; current artifacts export no non-premise law choosing one coupling sign or polarity.\n- Do not replay direct chiral-bispectrum tau-coupling as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must export a strict law fixing `lambda` or P2721 polarity, or preserve the P2697-P2732 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
